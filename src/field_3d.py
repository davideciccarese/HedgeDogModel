"""
field_3d.py  --  3D volumetric view of the O2 and metabolite fields
====================================================================
Model A is solved on a radial coordinate and is, by construction, radially
symmetric.  Here we lift the steady radial profiles into a 3-D cube and render
the hedgehog as a transparent dome (upper hemisphere) sitting on the tooth
plane, so the spatial distribution of each chemical species can be seen
directly.

For a field f, every voxel inside the dome is coloured with the species hue and
given an opacity proportional to the local concentration.  High-concentration
regions glow; depleted regions become transparent.  The result:

  * O2       : a bright translucent shell at the rim, a hollow (anoxic) core.
  * CO2      : a glowing band through the sub-rim annulus and interior.
  * lactate  : a thin bright peak just inside the rim, where Streptococcus
               ferment, decaying inward as it is consumed.

Outputs (all in ../figures/):
  D_field_O2_3d.gif        rotating O2 field, dome in a cube
  D_field_CO2_3d.gif       rotating CO2 field
  D_field_lactate_3d.gif   rotating lactate field
  D_fields_multipanel.gif  top-down O2 view + the three rotating 3-D fields
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.colors import to_rgb
import reaction_diffusion as rd
import params as P

FIG = os.path.join(os.path.dirname(__file__), "..", "figures")
os.makedirs(FIG, exist_ok=True)
RNG = np.random.default_rng(11)
R = rd.R

# ---- steady-state radial profiles from Model A ----------------------
_snaps = rd.simulate()
_, _O2, _LAC, _CO2 = _snaps[-1]
RG = rd.r
PROFILES = {
    "O$_2$":     dict(prof=_O2 / max(_O2.max(), 1e-9),  hue="#ff8a3d", gamma=1.9),
    "CO$_2$":    dict(prof=_CO2 / max(_CO2.max(), 1e-9), hue="#2ad4e6", gamma=1.5),
    "lactate":   dict(prof=_LAC / max(_LAC.max(), 1e-9), hue="#35d24b", gamma=1.3),
}

# ---- randomized volumetric samples in the upper hemisphere (the dome) ----
# Random (not lattice) points avoid the comb/striping artifact and read as mist.
_cand = RNG.uniform([-R, -R, 0], [R, R, R], size=(60000, 3))
_rc = np.sqrt((_cand**2).sum(1))
_keep = _rc <= R
VX, VY, VZ = _cand[_keep, 0], _cand[_keep, 1], _cand[_keep, 2]
VR = _rc[_keep]


def voxel_rgba(prof, hue, gamma=1.5, amax=0.40, thr=0.035):
    """Map a radial profile to per-voxel colour + opacity; drop near-empty voxels."""
    val = np.interp(VR, RG, prof, right=prof[-1])
    a = np.clip(val, 0, 1) ** gamma * amax
    keep = a > thr
    rgb = np.array(to_rgb(hue))
    rgba = np.empty((keep.sum(), 4))
    rgba[:, :3] = rgb
    rgba[:, 3] = a[keep]
    return VX[keep], VY[keep], VZ[keep], rgba


def dome_lines(n=70):
    """Faint Corynebacterium filaments for spatial context."""
    segs = []
    for _ in range(n):
        cz = RNG.uniform(0.05, 1.0)
        rr = np.sqrt(1 - cz * cz); th = RNG.uniform(0, 2 * np.pi)
        d = np.array([rr * np.cos(th), rr * np.sin(th), cz])
        segs.append((np.array([0, 0, 0]), d * R * 0.98))
    return segs


def draw_cube(ax):
    c = [-R, R]
    pts = [(x, y, z) for x in c for y in c for z in [0, R]]  # box: z in [0,R]
    edges = []
    for i, a in enumerate(pts):
        for b in pts[i + 1:]:
            if sum(1 for k in range(3) if a[k] != b[k]) == 1:
                edges.append((a, b))
    for a, b in edges:
        ax.plot(*zip(a, b), color="#2a3344", lw=0.8, alpha=0.7)
    # base (tooth) disc
    th = np.linspace(0, 2 * np.pi, 60)
    ax.plot(R * np.cos(th), R * np.sin(th), np.zeros_like(th), color="#3f72b0", lw=1.2, alpha=0.8)


def style3d(ax, title):
    ax.set_facecolor("#0a0e14")
    ax.set_xlim(-R, R); ax.set_ylim(-R, R); ax.set_zlim(0, R)
    ax.set_box_aspect((2, 2, 1.2))
    ax.set_title(title, color="#e7ecf3", fontsize=12, pad=2)
    ax.grid(False)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.set_pane_color((0.04, 0.055, 0.08, 1.0))
        axis.line.set_color("#2a3344")
    ax.tick_params(colors="#5b6678", labelsize=7)
    try:
        ax.set_xticks([-R, 0, R]); ax.set_yticks([-R, 0, R]); ax.set_zticks([0, R])
    except Exception:
        pass


def add_field(ax, name):
    info = PROFILES[name]
    xs, ys, zs, rgba = voxel_rgba(info["prof"], info["hue"], gamma=info["gamma"])
    for a, b in dome_lines():
        ax.plot(*zip(a, b), color=P.TAXA["Corynebacterium"]["color"], lw=0.5, alpha=0.16)
    ax.scatter(xs, ys, zs, c=rgba, s=9, marker="o", depthshade=False, edgecolors="none")
    draw_cube(ax)


# ---------------------------------------------------------------------
def single_gif(name, fname):
    fig = plt.figure(figsize=(6.4, 6.0), facecolor="#0a0e14")
    ax = fig.add_subplot(111, projection="3d")
    add_field(ax, name)
    style3d(ax, f"{name} field   (transparent: bright = high concentration)")
    fig.text(0.5, 0.04, "hedgehog dome in a cube  ·  base plane = tooth surface",
             color="#8b96a8", ha="center", fontsize=9)

    def update(f):
        ax.view_init(elev=22, azim=f * (360 / FRAMES))
        return []

    anim = animation.FuncAnimation(fig, update, frames=FRAMES, interval=80, blit=False)
    out = os.path.join(FIG, fname)
    anim.save(out, writer=animation.PillowWriter(fps=10))
    plt.close(fig)
    return out


def top_disc(prof, cmap):
    gg = np.linspace(-R, R, 240)
    X, Y = np.meshgrid(gg, gg)
    rr = np.sqrt(X**2 + Y**2)
    img = np.interp(rr, RG, prof, right=prof[-1])
    img[rr > R] = np.nan
    return img


def multipanel_gif():
    fig = plt.figure(figsize=(11.5, 10.5), facecolor="#0a0e14")
    # top-left: 2-D top-down view (kept from the original orientation)
    ax0 = fig.add_subplot(2, 2, 1)
    ax0.set_facecolor("#0a0e14")
    im = ax0.imshow(top_disc(PROFILES["O$_2$"]["prof"], "inferno"), origin="lower",
                    extent=[-R, R, -R, R], cmap="inferno", vmin=0, vmax=1)
    ax0.set_title("view from the top: O$_2$ (oxic rim, anoxic core)", color="#e7ecf3", fontsize=11)
    ax0.set_xlabel("x (\u03bcm)", color="#8b96a8"); ax0.set_ylabel("y (\u03bcm)", color="#8b96a8")
    ax0.tick_params(colors="#5b6678", labelsize=7)
    cb = fig.colorbar(im, ax=ax0, fraction=0.046, pad=0.04)
    cb.ax.yaxis.set_tick_params(color="#5b6678")
    plt.setp(plt.getp(cb.ax.axes, "yticklabels"), color="#8b96a8")

    ax_o2 = fig.add_subplot(2, 2, 2, projection="3d"); add_field(ax_o2, "O$_2$");   style3d(ax_o2, "O$_2$ field (3D)")
    ax_co = fig.add_subplot(2, 2, 3, projection="3d"); add_field(ax_co, "CO$_2$");  style3d(ax_co, "CO$_2$ field (3D)")
    ax_la = fig.add_subplot(2, 2, 4, projection="3d"); add_field(ax_la, "lactate"); style3d(ax_la, "lactate field (3D)")
    fig.suptitle("Metabolite distribution in the hedgehog consortium",
                 color="#e7ecf3", fontsize=14, y=0.97)

    def update(f):
        az = f * (360 / FRAMES)
        for ax in (ax_o2, ax_co, ax_la):
            ax.view_init(elev=22, azim=az)
        return []

    anim = animation.FuncAnimation(fig, update, frames=FRAMES, interval=80, blit=False)
    out = os.path.join(FIG, "D_fields_multipanel.gif")
    anim.save(out, writer=animation.PillowWriter(fps=10))
    plt.close(fig)
    return out


FRAMES = 72

if __name__ == "__main__":
    print("[D] rendering 3D metabolite fields ...")
    o2 = single_gif("O$_2$",   "D_field_O2_3d.gif")
    co = single_gif("CO$_2$",  "D_field_CO2_3d.gif")
    la = single_gif("lactate", "D_field_lactate_3d.gif")
    mp = multipanel_gif()
    for f in (o2, co, la, mp):
        print(f"    wrote {os.path.basename(f)}")
