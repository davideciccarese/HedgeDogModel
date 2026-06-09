"""
succession_3d.py  --  Model F: dynamic succession with evolving coupled fields
================================================================================
A time-resolved version of the coupled model.  The consortium colonises in
succession waves (early base colonisers -> radial Corynebacterium filaments ->
peripheral corncobs -> micro-aerophilic annulus), and at every frame the 3D
reaction-diffusion fields for O2, CO2 and lactate are re-solved against the
*current* biomass (warm-started from the previous frame).  Oxygen therefore
changes dynamically: the dome starts fully oxygenated and the anoxic core
emerges only as the consuming biomass accumulates.

The 3D fields are colour-coded by concentration (a colormap per species), not
by a single hue.  Outputs animate the succession and the field evolution:

  F_succession_topview.gif      top-down: colonisation waves + O2 + CO2 maps
  F_succession_multipanel.gif   3D: species succession + O2 + CO2 + lactate
  F_succession_O2_3d.gif        3D O2 field alone, colour-coded, evolving
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import animation
import scipy.ndimage as ndi
import params as P

FIG = os.path.join(os.path.dirname(__file__), "..", "figures")
os.makedirs(FIG, exist_ok=True)
RNG = np.random.default_rng(23)
R = P.R_HEDGEHOG
T = 48                                   # animation frames (= time steps of growth)
STEPS_PER_FRAME = 240                    # PDE relaxation steps per growth increment

# ---- grid ----------------------------------------------------------
NX = NY = 48
dx = 2*R/(NX-1)
ZMAX = 1.12*R
NZ = int(ZMAX/dx) + 1
xg = np.linspace(-R, R, NX); yg = np.linspace(-R, R, NY); zg = np.linspace(0, ZMAX, NZ)
GX, GY, GZ = np.meshgrid(xg, yg, zg, indexing="ij")
RR = np.sqrt(GX**2 + GY**2 + GZ**2)
INSIDE = RR <= R
EXT = ~INSIDE
IX, IY, IZ = GX[INSIDE], GY[INSIDE], GZ[INSIDE]
SUB = RNG.choice(IX.size, size=min(3200, IX.size), replace=False)   # fixed cloud subsample

RESP_W = {"Streptococcus": 1.0, "Haemophilus/Aggr.": 0.8, "Neisseriaceae": 0.7, "Corynebacterium": 0.6}
FERM_W = {"Streptococcus": 1.0}
CMAP = {"O2": "turbo", "CO2": "viridis", "lactate": "YlGn"}


def dome_dir(maxc=0.99):
    cz = RNG.uniform(0.05, maxc); s = np.sqrt(1-cz*cz); th = RNG.uniform(0, 2*np.pi)
    return np.array([s*np.cos(th), s*np.sin(th), cz])


# ---- timed biomass: every cell carries an appearance time t_app in [0,1] ----
def build_timed():
    pos, tax, tapp = [], [], []
    fils = []
    for _ in range(120):
        d = dome_dir(); L = RNG.uniform(0.85, 1.0)*R
        npts = 24; ss = np.linspace(0.05, 1.0, npts)
        pts = np.outer(ss, d)*L
        perp = np.cross(d, [0, 0, 1.0])
        if np.linalg.norm(perp) > 1e-6:
            perp /= np.linalg.norm(perp)
            pts = pts + np.outer(np.sin(np.linspace(0, np.pi, npts)), perp)*RNG.uniform(-9, 9)
        tp = 0.04 + 0.50*ss                       # filaments grow outward over time
        fils.append((pts, tp))
        for i in range(0, npts, 2):
            pos.append(pts[i]); tax.append("Corynebacterium"); tapp.append(tp[i])
        tip = pts[-1]
        if np.linalg.norm(tip) > 0.80*R and tip[2] >= 0 and RNG.random() < 0.7:
            axis = tip - pts[-3]; axis /= (np.linalg.norm(axis)+1e-9)
            u = np.cross(axis, [0, 0, 1.0]); u /= (np.linalg.norm(u)+1e-9); v = np.cross(axis, u)
            t_tip = tp[-1]
            m = RNG.integers(10, 24)
            for _ in range(m):
                ang = RNG.uniform(0, 2*np.pi); rad = RNG.uniform(3, 7.5)
                p = tip + axis*RNG.uniform(-7, 2) + (u*np.cos(ang)+v*np.sin(ang))*rad
                if p[2] >= 0:
                    pos.append(p); tax.append("Streptococcus"); tapp.append(min(1, t_tip+RNG.uniform(0.02, 0.10)))
            for _ in range(m//2):
                ang = RNG.uniform(0, 2*np.pi); rad = RNG.uniform(7.5, 11)
                p = tip + axis*RNG.uniform(-3, 4) + (u*np.cos(ang)+v*np.sin(ang))*rad
                if p[2] >= 0:
                    pos.append(p); tax.append("Haemophilus/Aggr."); tapp.append(min(1, t_tip+RNG.uniform(0.06, 0.16)))
            if RNG.random() < 0.5:
                pos.append(tip + u*RNG.uniform(-6, 6)+v*RNG.uniform(-6, 6))
                tax.append("Neisseriaceae"); tapp.append(min(1, t_tip+RNG.uniform(0.05, 0.15)))
    for _ in range(380):                          # annulus wave
        d = dome_dir(); p = d*RNG.uniform(0.5, 0.82)*R
        if p[2] >= 0:
            pos.append(p); tax.append(RNG.choice(["Fusobacterium", "Leptotrichia", "Capnocytophaga"]))
            tapp.append(RNG.uniform(0.45, 0.80))
    for _ in range(110):                          # early base colonisers
        d = dome_dir(0.6); p = d*RNG.uniform(2, 40)
        if p[2] >= 0:
            pos.append(p); tax.append(RNG.choice(["Actinomyces", "other"])); tapp.append(RNG.uniform(0.0, 0.12))
    return (np.array(pos), np.array(tax, dtype=object), np.array(tapp)), fils


def deposit(pos, tax, weights):
    field = np.zeros((NX, NY, NZ))
    for taxon, w in weights.items():
        m = tax == taxon
        if not m.any():
            continue
        p = pos[m]
        ix = np.clip(((p[:, 0]+R)/dx).round().astype(int), 0, NX-1)
        iy = np.clip(((p[:, 1]+R)/dx).round().astype(int), 0, NY-1)
        iz = np.clip((p[:, 2]/dx).round().astype(int), 0, NZ-1)
        np.add.at(field, (ix, iy, iz), w)
    return ndi.gaussian_filter(field, sigma=1.1)


def run_succession():
    (pos, tax, tapp), fils = build_timed()
    dx2 = dx*dx
    dt = 0.18*dx2/(6*max(P.D_O2, P.D_CO2, P.D_LAC))
    O = np.full((NX, NY, NZ), P.O2_SAT); Lf = np.zeros_like(O); Cf = np.zeros_like(O)
    # normalisation reference from the FINAL biomass so colours are comparable across time
    Wresp_full = deposit(pos, tax, RESP_W)
    resp_ref = np.percentile(Wresp_full[INSIDE][Wresp_full[INSIDE] > 0], 90)
    Wferm_full = deposit(pos, tax, FERM_W); ferm_ref = max(Wferm_full.max(), 1e-9)
    snaps = []
    metrics = dict(tf=[], anox=[], meanO2=[], meanCO2=[], meanLac=[], ncells=[])
    print("    marching succession + dynamic fields ...")
    for f in range(T):
        tf = (f+1)/T
        vis = tapp <= tf
        Wr = deposit(pos[vis], tax[vis], RESP_W)/resp_ref
        Wf = deposit(pos[vis], tax[vis], FERM_W)/ferm_ref
        for _ in range(STEPS_PER_FRAME):
            O += dt*(P.D_O2*ndi.laplace(O, mode="nearest")/dx2 - Wr*6.0*O/(P.KM_O2+O))
            Lf += dt*(P.D_LAC*ndi.laplace(Lf, mode="nearest")/dx2 + Wf*P.LAC_PROD*60 - 1.0*Lf)
            Cf += dt*(P.D_CO2*ndi.laplace(Cf, mode="nearest")/dx2 + Wf*P.CO2_PROD*60 - 0.15*Cf)
            O[EXT] = P.O2_SAT; Lf[EXT] = 0; Cf[EXT] = 0
            O[:, :, 0] = O[:, :, 1]; Lf[:, :, 0] = Lf[:, :, 1]; Cf[:, :, 0] = Cf[:, :, 1]
            np.clip(O, 0, None, out=O)
        metrics["tf"].append(tf)
        metrics["anox"].append(float((O[INSIDE] < P.ANOXIC_THR).mean()))
        metrics["meanO2"].append(float(O[INSIDE].mean()))
        metrics["meanCO2"].append(float(Cf[INSIDE].mean()))
        metrics["meanLac"].append(float(Lf[INSIDE].mean()))
        metrics["ncells"].append(int(vis.sum()))
        Ln = Lf/max(Lf.max(), 1e-9); Cn = Cf/max(Cf.max(), 1e-9)
        snaps.append(dict(tf=tf, O=O.astype(np.float32).copy(),
                          C=Cn.astype(np.float32).copy(), L=Ln.astype(np.float32).copy()))
        if f % 12 == 0 or f == T-1:
            print(f"      frame {f+1:2d}/{T}  t={tf:.2f}  visible cells={int(vis.sum())}  "
                  f"anoxic={ (O[INSIDE]<P.ANOXIC_THR).mean()*100:4.0f}%")
    metrics = {k: np.array(v) for k, v in metrics.items()}
    np.savez(os.path.join(FIG, "F_metrics.npz"), **metrics)
    return (pos, tax, tapp), fils, snaps, metrics


# ---- styling helpers -----------------------------------------------
def style3d(ax, title):
    ax.set_facecolor("#0a0e14"); ax.set_xlim(-R, R); ax.set_ylim(-R, R); ax.set_zlim(0, R)
    ax.set_box_aspect((2, 2, 1.2)); ax.grid(False); ax.set_title(title, color="#e7ecf3", fontsize=11, pad=2)
    for a in (ax.xaxis, ax.yaxis, ax.zaxis):
        a.set_pane_color((0.04, 0.055, 0.08, 1.0)); a.line.set_color("#2a3344")
    ax.tick_params(colors="#5b6678", labelsize=6); ax.set_xticks([-R, 0, R]); ax.set_yticks([-R, 0, R]); ax.set_zticks([0, R])


def cube(ax):
    c = [-R, R]; pts = [(x, y, z) for x in c for y in c for z in [0, R]]
    for i, a in enumerate(pts):
        for b in pts[i+1:]:
            if sum(1 for k in range(3) if a[k] != b[k]) == 1:
                ax.plot(*zip(a, b), color="#2a3344", lw=0.7, alpha=0.6)


def draw_cloud(ax, F, cmap, alpha=0.32):
    val = np.clip(F[INSIDE][SUB], 0, 1)
    rgba = plt.get_cmap(cmap)(val); rgba[:, 3] = alpha
    ax.scatter(IX[SUB], IY[SUB], IZ[SUB], c=rgba, s=10, depthshade=False, edgecolors="none")
    cube(ax)


ORDER = ["Actinomyces", "other", "Leptotrichia", "Fusobacterium", "Capnocytophaga",
         "Neisseriaceae", "Streptococcus", "Haemophilus/Aggr."]

# species -> spatial niche, for the legend panel
ZONE = {
    "Corynebacterium":   "filament backbone",
    "Actinomyces":       "base (early coloniser)",
    "other":             "base (sparse)",
    "Fusobacterium":     "annulus bridge",
    "Leptotrichia":      "annulus filament",
    "Capnocytophaga":    "annulus (capnophilic)",
    "Porphyromonas":     "corncob, anaerobe",
    "Neisseriaceae":     "perimeter clusters",
    "Streptococcus":     "corncob cocci (perimeter)",
    "Haemophilus/Aggr.": "corncob outer layer",
}
LEGEND_ORDER = ["Corynebacterium", "Actinomyces", "other", "Fusobacterium", "Leptotrichia",
                "Capnocytophaga", "Porphyromonas", "Neisseriaceae", "Streptococcus", "Haemophilus/Aggr."]


def draw_species_legend(ax):
    ax.set_facecolor("#0a0e14"); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.text(0.02, 0.97, "species", color="#e7ecf3", fontsize=12, weight="bold", va="top")
    for i, t in enumerate(LEGEND_ORDER):
        y = 0.88 - i*0.092
        ax.scatter(0.07, y, s=90, c=P.TAXA[t]["color"], edgecolors="none")
        ax.text(0.14, y+0.012, t, color="#e7ecf3", fontsize=9.5, va="center", style="italic")
        ax.text(0.14, y-0.028, ZONE[t], color="#8b96a8", fontsize=7.5, va="center")


def plot_biomass_cells(ax, pos, tax, fils, tf):
    """Draw the visible filaments + cells (by taxon) at colonisation time tf."""
    for pts, tp in fils:
        k = int((tp <= tf).sum())
        if k > 1:
            ax.plot(pts[:k, 0], pts[:k, 1], pts[:k, 2],
                    color=P.TAXA["Corynebacterium"]["color"], lw=0.5, alpha=0.42)
    vis = tax_tapp_cache[0] <= tf
    for taxon in ORDER:
        m = vis & (tax == taxon)
        if m.any():
            sz = 12 if taxon in ("Streptococcus", "Haemophilus/Aggr.", "Neisseriaceae") else 5
            ax.scatter(pos[m, 0], pos[m, 1], pos[m, 2], c=P.TAXA[taxon]["color"], s=sz,
                       depthshade=False, edgecolors="none")
    cube(ax)


tax_tapp_cache = [None]   # set by renderers so plot_biomass_cells can mask by appearance time


def phase(tf):
    if tf < 0.12: return "early colonisers settle on the pellicle"
    if tf < 0.35: return "Corynebacterium filaments grow outward"
    if tf < 0.55: return "corncobs assemble at the oxic tips"
    if tf < 0.78: return "micro-aerophilic annulus fills"
    return "mature consortium; anoxic core established"


# ---- single GIF: 3D species succession with a names panel ----
def gif_species(pos, tax, tapp, fils, snaps):
    tax_tapp_cache[0] = tapp
    fig = plt.figure(figsize=(9.4, 6.2), facecolor="#0a0e14")
    gs = fig.add_gridspec(1, 2, width_ratios=[3, 1])
    ax = fig.add_subplot(gs[0, 0], projection="3d")
    axl = fig.add_subplot(gs[0, 1])

    def update(f):
        s = snaps[f]; tf = s["tf"]; ax.cla(); axl.cla()
        plot_biomass_cells(ax, pos, tax, fils, tf)
        style3d(ax, f"species succession\n{phase(tf)}   (t = {tf:.2f})")
        ax.view_init(22, 30 + f*(150/T))
        draw_species_legend(axl)
        return []

    anim = animation.FuncAnimation(fig, update, frames=T, interval=110, blit=False)
    out = os.path.join(FIG, "F_succession_species_3d.gif")
    anim.save(out, writer=animation.PillowWriter(fps=9)); plt.close(fig)
    return out


# ---- single GIF: one colour-coded field, dynamic ----
def gif_single_field(snaps, key, cmap, label, fname):
    fig = plt.figure(figsize=(6.6, 6.2), facecolor="#0a0e14")
    ax = fig.add_subplot(111, projection="3d")
    sm = cm.ScalarMappable(cmap=cmap); sm.set_array([0, 1])
    cb = fig.colorbar(sm, ax=ax, fraction=0.04, pad=0.02)
    cb.set_label(label, color="#8b96a8")
    cb.ax.yaxis.set_tick_params(color="#5b6678"); plt.setp(plt.getp(cb.ax.axes, "yticklabels"), color="#8b96a8")

    def update(f):
        s = snaps[f]; ax.cla()
        draw_cloud(ax, s[key], cmap, alpha=0.34)
        style3d(ax, f"{label.split('(')[0].strip()}  (dynamic)   t = {s['tf']:.2f}")
        ax.view_init(22, 30 + f*(150/T))
        return []

    anim = animation.FuncAnimation(fig, update, frames=T, interval=110, blit=False)
    out = os.path.join(FIG, fname)
    anim.save(out, writer=animation.PillowWriter(fps=9)); plt.close(fig)
    return out# ---- GIF 1: top-down succession (colonisation waves + O2 + CO2) ----
def gif_topview(pos, tax, tapp, fils, snaps):
    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(15, 5.4), facecolor="#0a0e14")
    th = np.linspace(0, 2*np.pi, 80)

    def update(f):
        s = snaps[f]; tf = s["tf"]
        for ax in (a1, a2, a3):
            ax.clear(); ax.set_facecolor("#0a0e14"); ax.set_aspect("equal")
            ax.set_xlim(-R*1.05, R*1.05); ax.set_ylim(-R*1.05, R*1.05); ax.axis("off")
            ax.plot(R*np.cos(th), R*np.sin(th), color="#3f72b0", lw=1, ls="--", alpha=0.7)
        # colonisation
        for pts, tp in fils:
            k = int((tp <= tf).sum())
            if k > 1:
                a1.plot(pts[:k, 0], pts[:k, 1], color=P.TAXA["Corynebacterium"]["color"], lw=0.5, alpha=0.3)
        vismask = tapp <= tf
        for taxon in ORDER:
            m = vismask & (tax == taxon)
            if m.any():
                sz = 15 if taxon in ("Streptococcus", "Haemophilus/Aggr.", "Neisseriaceae") else 6
                a1.scatter(pos[m, 0], pos[m, 1], c=P.TAXA[taxon]["color"], s=sz, edgecolors="none")
        a1.set_title(f"colonisation (top view)\n{phase(tf)}", color="#e7ecf3", fontsize=11)
        # O2 / CO2 column means
        Otop = np.where(INSIDE, s["O"], np.nan).mean(axis=2)
        Ctop = np.where(INSIDE, s["C"], np.nan)
        Ctop = np.nansum(np.where(np.isnan(Ctop), 0, Ctop), axis=2)/np.maximum(INSIDE.sum(axis=2), 1)
        a2.imshow(np.nanmean(np.where(INSIDE, s["O"], np.nan), axis=2).T, origin="lower",
                  extent=[-R, R, -R, R], cmap=CMAP["O2"], vmin=0, vmax=1)
        a2.set_title(f"O$_2$ (top view)   t = {tf:.2f}", color="#e7ecf3", fontsize=11)
        a3.imshow(Ctop.T, origin="lower", extent=[-R, R, -R, R], cmap=CMAP["CO2"], vmin=0, vmax=1)
        a3.set_title(f"CO$_2$ (top view)   t = {tf:.2f}", color="#e7ecf3", fontsize=11)
        return []

    anim = animation.FuncAnimation(fig, update, frames=T, interval=110, blit=False)
    out = os.path.join(FIG, "F_succession_topview.gif")
    anim.save(out, writer=animation.PillowWriter(fps=9)); plt.close(fig)
    return out


# ---- GIF 2: 3D multipanel succession (biomass + O2 + CO2 + lactate + names + dynamics) ----
def gif_multipanel(pos, tax, tapp, fils, snaps, metrics):
    tax_tapp_cache[0] = tapp
    fig = plt.figure(figsize=(16, 9), facecolor="#0a0e14")
    gs = fig.add_gridspec(2, 3)
    ab = fig.add_subplot(gs[0, 0], projection="3d")
    ao = fig.add_subplot(gs[0, 1], projection="3d")
    ac = fig.add_subplot(gs[0, 2], projection="3d")
    al = fig.add_subplot(gs[1, 0], projection="3d")
    axn = fig.add_subplot(gs[1, 1])
    axd = fig.add_subplot(gs[1, 2])

    def update(f):
        s = snaps[f]; tf = s["tf"]; az = 30 + f*(150/T)
        for ax in (ab, ao, ac, al):
            ax.cla()
        plot_biomass_cells(ab, pos, tax, fils, tf); style3d(ab, f"species succession\n{phase(tf)}")
        draw_cloud(ao, s["O"], CMAP["O2"]);      style3d(ao, "O$_2$ (dynamic)")
        draw_cloud(ac, s["C"], CMAP["CO2"]);     style3d(ac, "CO$_2$")
        draw_cloud(al, s["L"], CMAP["lactate"]); style3d(al, "lactate")
        for ax in (ab, ao, ac, al):
            ax.view_init(22, az)
        axn.cla(); draw_species_legend(axn)
        axd.cla(); axd.set_facecolor("#0a0e14")
        tfa = metrics["tf"]; k = f+1
        axd.plot(tfa, metrics["anox"]*100, color="#3a4763", lw=1.2)
        axd.plot(tfa[:k], metrics["anox"][:k]*100, color="#2ad4e6", lw=2.6)
        axd.scatter([tf], [metrics["anox"][f]*100], c="#ffffff", s=24, zorder=5)
        axd.set_xlim(0, 1); axd.set_ylim(0, max(metrics["anox"].max()*100*1.1, 1))
        axd.set_title("anoxic volume (%)", color="#e7ecf3", fontsize=11)
        axd.set_xlabel("colonisation time t", color="#8b96a8")
        axd.tick_params(colors="#5b6678", labelsize=8)
        for sp in axd.spines.values(): sp.set_color("#2a3344")
        fig.suptitle(f"Dynamic succession and coupled fields   t = {tf:.2f}",
                     color="#e7ecf3", fontsize=15, y=0.98)
        return []

    anim = animation.FuncAnimation(fig, update, frames=T, interval=110, blit=False)
    out = os.path.join(FIG, "F_succession_multipanel.gif")
    anim.save(out, writer=animation.PillowWriter(fps=9)); plt.close(fig)
    return out


def plot_dynamics(metrics):
    tf = metrics["tf"]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.4), facecolor="white")
    # anoxic fraction + colonisation
    a1b = a1.twinx()
    a1.plot(tf, metrics["anox"]*100, color="#2233e0", lw=2.4, label="anoxic volume fraction")
    a1b.plot(tf, metrics["ncells"], color="#8a93a0", lw=1.8, ls="--", label="visible cells")
    a1.set_xlabel("colonisation time  t (normalised)"); a1.set_ylabel("anoxic volume (%)", color="#2233e0")
    a1b.set_ylabel("visible cells", color="#666")
    a1.set_title("Oxygen depletes as the consortium colonises")
    a1.tick_params(axis="y", colors="#2233e0")
    # normalised mean fields
    def nrm(x): return x/max(x.max(), 1e-9)
    a2.plot(tf, metrics["meanO2"], color="#ff8a3d", lw=2.4, label="mean O$_2$")
    a2.plot(tf, nrm(metrics["meanCO2"]), color="#2ad4e6", lw=2, label="mean CO$_2$ (norm.)")
    a2.plot(tf, nrm(metrics["meanLac"]), color="#35d24b", lw=2, label="mean lactate (norm.)")
    a2.set_xlabel("colonisation time  t (normalised)"); a2.set_ylabel("mean concentration")
    a2.set_title("Nutrient and oxygen succession"); a2.legend(fontsize=8)
    fig.suptitle("Model F: dynamic succession analysis", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = os.path.join(FIG, "F_dynamics.png"); fig.savefig(out, dpi=130); plt.close(fig)
    return out


if __name__ == "__main__":
    print("[F] dynamic succession with evolving coupled fields ...")
    (pos, tax, tapp), fils, snaps, metrics = run_succession()
    for fn in (gif_topview(pos, tax, tapp, fils, snaps),
               gif_multipanel(pos, tax, tapp, fils, snaps, metrics),
               gif_species(pos, tax, tapp, fils, snaps),
               gif_single_field(snaps, "O", CMAP["O2"], "O$_2$ (fraction of saliva)", "F_succession_O2_3d.gif"),
               gif_single_field(snaps, "C", CMAP["CO2"], "CO$_2$ (normalised)", "F_succession_CO2_3d.gif"),
               gif_single_field(snaps, "L", CMAP["lactate"], "lactate (normalised)", "F_succession_lactate_3d.gif"),
               plot_dynamics(metrics)):
        print(f"    wrote {os.path.basename(fn)}")
