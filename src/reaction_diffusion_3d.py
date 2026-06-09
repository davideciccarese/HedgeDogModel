"""
reaction_diffusion_3d.py  --  Model E: a genuine 3D field coupled to the biomass
=================================================================================
Unlike `field_3d.py`, which lifted a 1D radial solution into 3D under an exact
symmetry assumption, this module solves the reaction-diffusion equations on the
full 3D voxel grid, with the source and sink terms supplied by *discrete,
heterogeneous biomass*.  The field is therefore no longer radially symmetric:
oxygen dips locally beneath dense corncob clusters, and the metabolite plumes
follow the actual distribution of fermenting cells.

Pipeline
--------
1. Build a 3D hedgehog: ~140 Corynebacterium filaments radiating into the upper
   hemisphere, corncob clusters of Streptococcus + Haemophilus at a random subset
   of filament tips, a sub-rim annulus of microaerophilic rods, and a sparse base.
2. Rasterise two biomass fields onto the grid and smooth them:
       W_resp  (O2 consumers: Strep, Haem, Neisseriaceae, Corynebacterium)
       W_ferm  (lactate / CO2 producers: Streptococcus)
3. March three coupled PDEs to quasi-steady state on the voxel grid:
       dO/dt = D_O lap(O) - W_resp * Qmax * O/(Km+O)
       dL/dt = D_L lap(L) + W_ferm * Lp - kL L
       dC/dt = D_C lap(C) + W_ferm * Cp - kC C
   Saliva (O2 = O2_sat, metabolites washed) bathes the exposed dome surface
   (a reservoir at r > R); the tooth plane at z = 0 is impermeable (no flux).

Outputs (../figures/):
  E_biomass_3d.gif          the discrete 3D biomass
  E_field3d_O2.gif          rotating transparent O2 field (heterogeneous)
  E_field3d_CO2.gif         rotating transparent CO2 field
  E_field3d_lactate.gif     rotating transparent lactate field
  E_slices.png              vertical + horizontal O2 slices with biomass overlaid
  E_heterogeneity.png       O2 around a ring: coupled (wavy) vs symmetric (flat)
  E_coupled_multipanel.gif  biomass + O2 cloud (rotating) + the two slices
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.colors import to_rgb
import scipy.ndimage as ndi
import params as P

FIG = os.path.join(os.path.dirname(__file__), "..", "figures")
os.makedirs(FIG, exist_ok=True)
RNG = np.random.default_rng(19)
R = P.R_HEDGEHOG
FRAMES = 72

# ---------------------------------------------------------------------
# 1. discrete 3D biomass
# ---------------------------------------------------------------------
def dome_dir(maxc=0.99):
    cz = RNG.uniform(0.05, maxc)
    s = np.sqrt(1 - cz*cz); th = RNG.uniform(0, 2*np.pi)
    return np.array([s*np.cos(th), s*np.sin(th), cz])


def build_biomass():
    """Return dict of taxon -> Nx3 positions, plus filament polylines."""
    cells = {k: [] for k in ("Streptococcus", "Haemophilus/Aggr.", "Neisseriaceae",
                              "Fusobacterium", "Leptotrichia", "Capnocytophaga",
                              "Actinomyces", "other", "Corynebacterium")}
    fils = []
    for _ in range(140):
        d = dome_dir()
        L = RNG.uniform(0.85, 1.0) * R
        pts = np.outer(np.linspace(0.05, 1.0, 26), d) * L
        # gentle bow
        perp = np.cross(d, [0, 0, 1]); 
        if np.linalg.norm(perp) > 1e-6:
            perp /= np.linalg.norm(perp)
            pts = pts + np.outer(np.sin(np.linspace(0, np.pi, 26)), perp) * RNG.uniform(-10, 10)
        fils.append(pts)
        cells["Corynebacterium"].extend(pts[::2])
        # corncob at the oxic tip on a random subset of filaments
        if pts[-1][2] >= 0 and np.linalg.norm(pts[-1]) > 0.85*R and RNG.random() < 0.7:
            tip = pts[-1]; axis = (tip - pts[-3]); axis /= (np.linalg.norm(axis)+1e-9)
            u = np.cross(axis, [0, 0, 1.0]); u /= (np.linalg.norm(u)+1e-9)
            v = np.cross(axis, u)
            m = RNG.integers(10, 26)                 # variable cluster size -> heterogeneity
            for _ in range(m):
                ang = RNG.uniform(0, 2*np.pi); rad = RNG.uniform(3.0, 7.5)
                along = RNG.uniform(-7, 2)
                p = tip + axis*along + (u*np.cos(ang)+v*np.sin(ang))*rad
                if p[2] >= 0:
                    cells["Streptococcus"].append(p)
            for _ in range(m//2):
                ang = RNG.uniform(0, 2*np.pi); rad = RNG.uniform(7.5, 11)
                p = tip + axis*RNG.uniform(-3, 4) + (u*np.cos(ang)+v*np.sin(ang))*rad
                if p[2] >= 0:
                    cells["Haemophilus/Aggr."].append(p)
            if RNG.random() < 0.5:
                cells["Neisseriaceae"].append(tip + (u*RNG.uniform(-6,6)+v*RNG.uniform(-6,6)))
    # sub-rim microaerophilic annulus
    for _ in range(420):
        d = dome_dir(); rr = RNG.uniform(0.5, 0.82)*R; p = d*rr
        if p[2] >= 0:
            cells[RNG.choice(["Fusobacterium", "Leptotrichia", "Capnocytophaga"])].append(p)
    # sparse base
    for _ in range(120):
        d = dome_dir(0.6); rr = RNG.uniform(2, 40); p = d*rr
        if p[2] >= 0:
            cells[RNG.choice(["Actinomyces", "other"])].append(p)
    return {k: np.array(v) for k, v in cells.items() if len(v)}, fils


# ---------------------------------------------------------------------
# 2. rasterise biomass -> grid source/sink fields
# ---------------------------------------------------------------------
NX = NY = 58
dx = 2*R / (NX - 1)
ZMAX = 1.12 * R
NZ = int(ZMAX / dx) + 1
xg = np.linspace(-R, R, NX)
yg = np.linspace(-R, R, NY)
zg = np.linspace(0, ZMAX, NZ)
GX, GY, GZ = np.meshgrid(xg, yg, zg, indexing="ij")
RR = np.sqrt(GX**2 + GY**2 + GZ**2)
INSIDE = RR <= R
EXT = ~INSIDE                                  # saliva reservoir surrounding the dome

RESP_W = {"Streptococcus": 1.0, "Haemophilus/Aggr.": 0.8, "Neisseriaceae": 0.7,
          "Corynebacterium": 0.6}
FERM_W = {"Streptococcus": 1.0}


def deposit(cells, weights):
    field = np.zeros((NX, NY, NZ))
    for taxon, w in weights.items():
        if taxon not in cells:
            continue
        pos = cells[taxon]
        ix = np.clip(((pos[:, 0] + R) / dx).round().astype(int), 0, NX-1)
        iy = np.clip(((pos[:, 1] + R) / dx).round().astype(int), 0, NY-1)
        iz = np.clip((pos[:, 2] / dx).round().astype(int), 0, NZ-1)
        np.add.at(field, (ix, iy, iz), w)
    field = ndi.gaussian_filter(field, sigma=1.1)
    return field


# ---------------------------------------------------------------------
# 3. solve the coupled 3D PDEs
# ---------------------------------------------------------------------
def solve(Wresp, Wferm):
    # normalise biomass weights to a sensible reaction regime
    Wresp = Wresp / np.percentile(Wresp[INSIDE][Wresp[INSIDE] > 0], 90)
    Wferm = Wferm / max(Wferm.max(), 1e-9)
    dx2 = dx*dx
    dt = 0.18 * dx2 / (6 * max(P.D_O2, P.D_CO2, P.D_LAC))
    steps = int(4.0 / dt)

    O = np.where(EXT, P.O2_SAT, P.O2_SAT).astype(float)
    Lf = np.zeros_like(O); Cf = np.zeros_like(O)
    Qmax = 6.0                                   # 3D: O2 enters from all sides, needs stronger sink
    for _ in range(steps):
        O += dt * (P.D_O2 * ndi.laplace(O, mode="nearest")/dx2
                   - Wresp * Qmax * O/(P.KM_O2 + O))
        Lf += dt * (P.D_LAC * ndi.laplace(Lf, mode="nearest")/dx2
                    + Wferm * P.LAC_PROD*60 - 1.0*Lf)
        Cf += dt * (P.D_CO2 * ndi.laplace(Cf, mode="nearest")/dx2
                    + Wferm * P.CO2_PROD*60 - 0.15*Cf)
        O[EXT] = P.O2_SAT; Lf[EXT] = 0.0; Cf[EXT] = 0.0
        O[:, :, 0] = O[:, :, 1]                  # impermeable enamel (no flux at base)
        Lf[:, :, 0] = Lf[:, :, 1]; Cf[:, :, 0] = Cf[:, :, 1]
        np.clip(O, 0, None, out=O)
    return O, Lf/max(Lf.max(), 1e-9), Cf/max(Cf.max(), 1e-9)


# ---------------------------------------------------------------------
# 4. rendering
# ---------------------------------------------------------------------
def cube_edges(ax):
    c = [-R, R]
    pts = [(x, y, z) for x in c for y in c for z in [0, R]]
    for i, a in enumerate(pts):
        for b in pts[i+1:]:
            if sum(1 for k in range(3) if a[k] != b[k]) == 1:
                ax.plot(*zip(a, b), color="#2a3344", lw=0.8, alpha=0.7)
    th = np.linspace(0, 2*np.pi, 60)
    ax.plot(R*np.cos(th), R*np.sin(th), np.zeros_like(th), color="#3f72b0", lw=1.2, alpha=0.8)


def style3d(ax, title):
    ax.set_facecolor("#0a0e14")
    ax.set_xlim(-R, R); ax.set_ylim(-R, R); ax.set_zlim(0, R)
    ax.set_box_aspect((2, 2, 1.2)); ax.grid(False)
    ax.set_title(title, color="#e7ecf3", fontsize=12, pad=2)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_pane_color((0.04, 0.055, 0.08, 1.0)); axis.line.set_color("#2a3344")
    ax.tick_params(colors="#5b6678", labelsize=7)
    ax.set_xticks([-R, 0, R]); ax.set_yticks([-R, 0, R]); ax.set_zticks([0, R])


def field_points(F, cmap, gamma=1.7, amax=0.5, thr=0.035):
    val = F[INSIDE]
    a = np.clip(val, 0, 1)**gamma * amax
    keep = a > thr
    xs, ys, zs = GX[INSIDE][keep], GY[INSIDE][keep], GZ[INSIDE][keep]
    rgba = plt.get_cmap(cmap)(np.clip(val[keep], 0, 1)); rgba[:, 3] = a[keep]
    return xs, ys, zs, rgba


def gif_field(F, cmap, title, fname, gamma=1.7):
    fig = plt.figure(figsize=(6.6, 6.0), facecolor="#0a0e14")
    ax = fig.add_subplot(111, projection="3d")
    xs, ys, zs, rgba = field_points(F, cmap, gamma=gamma)
    ax.scatter(xs, ys, zs, c=rgba, s=10, depthshade=False, edgecolors="none")
    cube_edges(ax); style3d(ax, title)
    sm = matplotlib.cm.ScalarMappable(cmap=cmap); sm.set_array([0, 1])
    cb = fig.colorbar(sm, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label("concentration (normalised)", color="#8b96a8")
    cb.ax.yaxis.set_tick_params(color="#5b6678"); plt.setp(plt.getp(cb.ax.axes, "yticklabels"), color="#8b96a8")
    fig.text(0.5, 0.04, "3D solve coupled to discrete biomass  ·  field is not radially symmetric",
             color="#8b96a8", ha="center", fontsize=8.5)
    anim = animation.FuncAnimation(fig, lambda f: ax.view_init(22, f*(360/FRAMES)) or [],
                                   frames=FRAMES, interval=80, blit=False)
    out = os.path.join(FIG, fname)
    anim.save(out, writer=animation.PillowWriter(fps=10)); plt.close(fig)
    return out


def gif_biomass(cells, fils, fname):
    fig = plt.figure(figsize=(6.4, 6.0), facecolor="#0a0e14")
    ax = fig.add_subplot(111, projection="3d")
    for pts in fils:
        ax.plot(pts[:, 0], pts[:, 1], pts[:, 2], color=P.TAXA["Corynebacterium"]["color"],
                lw=0.6, alpha=0.5)
    for taxon, pos in cells.items():
        if taxon == "Corynebacterium" or not len(pos):
            continue
        s = 14 if taxon in ("Streptococcus", "Haemophilus/Aggr.", "Neisseriaceae") else 6
        ax.scatter(pos[:, 0], pos[:, 1], pos[:, 2], c=P.TAXA[taxon]["color"], s=s,
                   depthshade=False, edgecolors="none")
    cube_edges(ax); style3d(ax, "discrete 3D biomass (sources & sinks)")
    anim = animation.FuncAnimation(fig, lambda f: ax.view_init(22, f*(360/FRAMES)) or [],
                                   frames=FRAMES, interval=80, blit=False)
    out = os.path.join(FIG, fname)
    anim.save(out, writer=animation.PillowWriter(fps=10)); plt.close(fig)
    return out


def slices_png(O, cells):
    iy0 = NY//2; izm = int(0.45*R/dx)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.8), facecolor="white")
    # vertical slice (x,z) at y=0  -- the classic wedge cross-section
    sl = O[:, iy0, :].T
    im1 = a1.imshow(sl, origin="lower", extent=[-R, R, 0, ZMAX], cmap="inferno",
                    vmin=0, vmax=P.O2_SAT, aspect="auto")
    a1.contour(np.linspace(-R, R, NX), zg, sl, levels=[P.ANOXIC_THR], colors="cyan", linewidths=1)
    th = np.linspace(0, np.pi, 80); a1.plot(R*np.cos(th), R*np.sin(th), color="#3f72b0", lw=1)
    a1.set_title("vertical slice (y=0): O$_2$ with anoxic contour"); a1.set_xlabel("x (\u03bcm)"); a1.set_ylabel("z (\u03bcm)")
    fig.colorbar(im1, ax=a1, fraction=0.046, pad=0.04)
    # horizontal slice (x,y) at mid height, with corncob clusters overlaid
    sl2 = O[:, :, izm].T
    im2 = a2.imshow(sl2, origin="lower", extent=[-R, R, -R, R], cmap="inferno",
                    vmin=0, vmax=P.O2_SAT)
    if "Streptococcus" in cells:
        sp = cells["Streptococcus"]
        near = np.abs(sp[:, 2] - izm*dx) < dx*2
        a2.scatter(sp[near, 0], sp[near, 1], s=6, c="#35d24b", alpha=.7, label="Streptococcus")
        a2.legend(loc="upper right", fontsize=7)
    a2.set_title(f"horizontal slice (z={izm*dx:.0f} \u03bcm): O$_2$ dips under clusters")
    a2.set_xlabel("x (\u03bcm)"); a2.set_ylabel("y (\u03bcm)")
    fig.colorbar(im2, ax=a2, fraction=0.046, pad=0.04)
    fig.suptitle("Heterogeneous O$_2$ field from the coupled 3D solve", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = os.path.join(FIG, "E_slices.png"); fig.savefig(out, dpi=130); plt.close(fig)
    return out


def heterogeneity_png(O):
    """O2 sampled around a horizontal ring vs azimuth: coupled is wavy, symmetric is flat."""
    izm = int(0.5*R/dx); r_ring = 0.78*R
    th = np.linspace(0, 2*np.pi, 240)
    xs = r_ring*np.cos(th); ys = r_ring*np.sin(th); zs = np.full_like(th, izm*dx)
    ix = np.clip(((xs+R)/dx).round().astype(int), 0, NX-1)
    iy = np.clip(((ys+R)/dx).round().astype(int), 0, NY-1)
    vals = O[ix, iy, izm]
    fig, ax = plt.subplots(figsize=(8, 4.2), facecolor="white")
    ax.plot(np.degrees(th), vals, color="#ff8a3d", lw=2, label="coupled 3D solve (this model)")
    ax.axhline(vals.mean(), color="#888", ls="--", lw=1.5,
               label="radially symmetric model (azimuth-independent)")
    ax.set_xlabel("azimuth around the ring  (degrees)")
    ax.set_ylabel("O$_2$  (fraction of saliva value)")
    ax.set_title(f"Azimuthal heterogeneity at r = {r_ring:.0f} \u03bcm\n"
                 f"coupling to discrete biomass produces real angular variation "
                 f"(CV = {vals.std()/vals.mean()*100:.0f}%)")
    ax.legend(); ax.set_xlim(0, 360)
    fig.tight_layout()
    out = os.path.join(FIG, "E_heterogeneity.png"); fig.savefig(out, dpi=130); plt.close(fig)
    return out


def multipanel(cells, fils, O):
    fig = plt.figure(figsize=(11.5, 10.5), facecolor="#0a0e14")
    ax_b = fig.add_subplot(2, 2, 1, projection="3d")
    for pts in fils:
        ax_b.plot(pts[:, 0], pts[:, 1], pts[:, 2], color=P.TAXA["Corynebacterium"]["color"], lw=0.5, alpha=0.4)
    for taxon, pos in cells.items():
        if taxon == "Corynebacterium" or not len(pos):
            continue
        s = 12 if taxon in ("Streptococcus", "Haemophilus/Aggr.", "Neisseriaceae") else 5
        ax_b.scatter(pos[:, 0], pos[:, 1], pos[:, 2], c=P.TAXA[taxon]["color"], s=s, depthshade=False, edgecolors="none")
    cube_edges(ax_b); style3d(ax_b, "discrete biomass")

    ax_o = fig.add_subplot(2, 2, 2, projection="3d")
    xs, ys, zs, rgba = field_points(O/max(O.max(), 1e-9), "turbo")
    ax_o.scatter(xs, ys, zs, c=rgba, s=9, depthshade=False, edgecolors="none")
    cube_edges(ax_o); style3d(ax_o, "O$_2$ field (3D, coupled)")

    iy0 = NY//2; izm = int(0.45*R/dx)
    ax1 = fig.add_subplot(2, 2, 3); ax1.set_facecolor("#0a0e14")
    sl = O[:, iy0, :].T
    ax1.imshow(sl, origin="lower", extent=[-R, R, 0, ZMAX], cmap="inferno", vmin=0, vmax=P.O2_SAT, aspect="auto")
    ax1.contour(np.linspace(-R, R, NX), zg, sl, levels=[P.ANOXIC_THR], colors="cyan", linewidths=1)
    ax1.set_title("vertical O$_2$ slice (y=0)", color="#e7ecf3", fontsize=11)
    ax1.set_xlabel("x (\u03bcm)", color="#8b96a8"); ax1.set_ylabel("z (\u03bcm)", color="#8b96a8")
    ax1.tick_params(colors="#5b6678", labelsize=7)

    ax2 = fig.add_subplot(2, 2, 4); ax2.set_facecolor("#0a0e14")
    ax2.imshow(O[:, :, izm].T, origin="lower", extent=[-R, R, -R, R], cmap="inferno", vmin=0, vmax=P.O2_SAT)
    if "Streptococcus" in cells:
        sp = cells["Streptococcus"]; near = np.abs(sp[:, 2]-izm*dx) < dx*2
        ax2.scatter(sp[near, 0], sp[near, 1], s=5, c="#35d24b", alpha=.7)
    ax2.set_title(f"horizontal O$_2$ slice (z={izm*dx:.0f} \u03bcm)", color="#e7ecf3", fontsize=11)
    ax2.set_xlabel("x (\u03bcm)", color="#8b96a8"); ax2.set_ylabel("y (\u03bcm)", color="#8b96a8")
    ax2.tick_params(colors="#5b6678", labelsize=7)

    fig.suptitle("Coupled 3D reaction-diffusion: the field follows the cells",
                 color="#e7ecf3", fontsize=14, y=0.97)

    def update(f):
        for ax in (ax_b, ax_o):
            ax.view_init(22, f*(360/FRAMES))
        return []
    anim = animation.FuncAnimation(fig, update, frames=FRAMES, interval=80, blit=False)
    out = os.path.join(FIG, "E_coupled_multipanel.gif")
    anim.save(out, writer=animation.PillowWriter(fps=10)); plt.close(fig)
    return out


def topview_png(cells, fils, O):
    """View from the top (looking down the z-axis): colonization + O2 column map."""
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11.4, 5.6), facecolor="#0a0e14")
    th = np.linspace(0, 2*np.pi, 80)
    for ax in (a1, a2):
        ax.set_facecolor("#0a0e14"); ax.set_aspect("equal")
        ax.set_xlim(-R*1.05, R*1.05); ax.set_ylim(-R*1.05, R*1.05); ax.axis("off")
        ax.plot(R*np.cos(th), R*np.sin(th), color="#3f72b0", lw=1.1, ls="--", alpha=0.8)

    # --- colonization, top view ---
    for pts in fils:
        a1.plot(pts[:, 0], pts[:, 1], color=P.TAXA["Corynebacterium"]["color"], lw=0.5, alpha=0.30)
    order = ["Actinomyces", "other", "Leptotrichia", "Fusobacterium", "Capnocytophaga",
             "Neisseriaceae", "Streptococcus", "Haemophilus/Aggr."]
    for taxon in order:
        pos = cells.get(taxon)
        if pos is None or not len(pos):
            continue
        s = 16 if taxon in ("Streptococcus", "Haemophilus/Aggr.", "Neisseriaceae") else 7
        a1.scatter(pos[:, 0], pos[:, 1], c=P.TAXA[taxon]["color"], s=s, edgecolors="none", label=taxon)
    a1.set_title("colonization, view from the top", color="#e7ecf3", fontsize=12)
    a1.legend(loc="upper right", fontsize=6.2, facecolor="#10151f", edgecolor="#2a3344",
              labelcolor="#cdd6e4", framealpha=.9, markerscale=1.3)

    # --- O2, top view (mean through the column, dome only) ---
    Omask = np.where(INSIDE, O, np.nan)
    top = np.nanmean(Omask, axis=2).T                 # rows = y, cols = x
    im = a2.imshow(top, origin="lower", extent=[-R, R, -R, R], cmap="inferno",
                   vmin=0, vmax=P.O2_SAT)
    a2.set_title("O$_2$, view from the top (column mean)", color="#e7ecf3", fontsize=12)
    cb = fig.colorbar(im, ax=a2, fraction=0.046, pad=0.04)
    cb.set_label("O$_2$ (fraction of saliva value)", color="#8b96a8")
    cb.ax.yaxis.set_tick_params(color="#5b6678")
    plt.setp(plt.getp(cb.ax.axes, "yticklabels"), color="#8b96a8")

    fig.suptitle("Top-down view of the coupled 3D hedgehog", color="#e7ecf3", fontsize=13, y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = os.path.join(FIG, "E_topview.png")
    fig.savefig(out, dpi=130, facecolor="#0a0e14"); plt.close(fig)
    return out


if __name__ == "__main__":
    print("[E] building discrete 3D biomass ...")
    cells, fils = build_biomass()
    nctot = sum(len(v) for v in cells.values())
    print(f"    biomass: {nctot} cells across {len(cells)} taxa, {len(fils)} filaments")
    print("    rasterising sources/sinks and solving 3D PDEs ...")
    Wresp = deposit(cells, RESP_W); Wferm = deposit(cells, FERM_W)
    O, Lf, Cf = solve(Wresp, Wferm)
    anoxic_frac = (O[INSIDE] < P.ANOXIC_THR).mean()
    print(f"    anoxic fraction of dome volume: {anoxic_frac*100:.0f}%")
    for fn in (
        gif_biomass(cells, fils, "E_biomass_3d.gif"),
        gif_field(O/max(O.max(), 1e-9), "turbo", "O$_2$ field (3D, coupled)", "E_field3d_O2.gif", gamma=1.9),
        gif_field(Cf, "viridis", "CO$_2$ field (3D, coupled)", "E_field3d_CO2.gif", gamma=1.5),
        gif_field(Lf, "YlGn", "lactate field (3D, coupled)", "E_field3d_lactate.gif", gamma=1.3),
        slices_png(O, cells), heterogeneity_png(O), multipanel(cells, fils, O),
        topview_png(cells, fils, O),
    ):
        print(f"    wrote {os.path.basename(fn)}")
