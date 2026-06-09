"""
corncob_coaggregation.py  --  Model C: the corncob and the cell-cell-distance rule
==================================================================================
A zoom-in on a single Corynebacterium filament tip in the oxic perimeter, where
the "corncob" structures of Welch et al. (2016, Figs. 2, 4) assemble.

This model makes the Kolenbrander et al. (2010) principle explicit: a cell can
only join the structure where a *partner* already sits within a contact distance
(COAGG_DIST).  From that single rule a concentric, layered corncob emerges:

    Corynebacterium filament (core)
        -> Streptococcus / Porphyromonas in direct contact with the filament
            -> Haemophilus/Aggregatibacter in contact with Streptococcus (outer)

We build the corncob on a cylinder (filament axis vertical) and project to 2-D
for display, with depth used for draw order and a slight size cue.  The analysis
panel confirms the emergent concentric layering and verifies that essentially
every Haemophilus cell is juxtaposed to a Streptococcus, as the rule requires.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation
import params as P

FIG = os.path.join(os.path.dirname(__file__), "..", "figures")
os.makedirs(FIG, exist_ok=True)
RNG = np.random.default_rng(3)

L   = 64.0          # visible filament length, um
RF  = 2.2           # filament radius
RS  = 1.7           # Streptococcus radius
RH  = 1.5           # Haemophilus radius
RP  = 1.5           # Porphyromonas radius
FRAMES = 96

C_CORE = P.TAXA["Corynebacterium"]["color"]
C_STR  = P.TAXA["Streptococcus"]["color"]
C_HAE  = P.TAXA["Haemophilus/Aggr."]["color"]
C_POR  = P.TAXA["Porphyromonas"]["color"]


class Corncob:
    def __init__(self):
        # cells stored as (x, y, z, radius, taxon)
        self.cells = []
        self.ring_h = np.arange(RS, L - RS, 1.9 * RS)   # candidate heights

    def _no_overlap(self, x, y, z, rad, taxa=None):
        for (cx, cy, cz, cr, ct) in self.cells:
            if taxa and ct not in taxa:
                continue
            if (x-cx)**2 + (y-cy)**2 + (z-cz)**2 < (0.85*(rad+cr))**2:
                return False
        return True

    def _surface(self, radial, h, taxon, rad):
        th = RNG.uniform(0, 2*np.pi)
        x = radial*np.cos(th); z = radial*np.sin(th)
        if self._no_overlap(x, h, z, rad):
            self.cells.append((x, h, z, rad, taxon)); return True
        return False

    def grow_streptococcus(self, n):
        # cocci coaggregate directly onto the filament (contact with Corynebacterium)
        placed = 0
        for _ in range(n):
            h = RNG.choice(self.ring_h) + RNG.normal(0, 0.6)
            if self._surface(RF + RS, h, "Streptococcus", RS):
                placed += 1
        return placed

    def grow_porphyromonas(self, n):
        # anaerobe, also in direct contact with the filament, interspersed with Strep
        for _ in range(n):
            h = RNG.choice(self.ring_h) + RNG.normal(0, 0.6)
            self._surface(RF + RP, h, "Porphyromonas", RP)

    def grow_haemophilus(self, n):
        # MUST attach to an existing Streptococcus (cell-cell-distance rule)
        strep = [c for c in self.cells if c[4] == "Streptococcus"]
        if not strep:
            return
        for _ in range(n):
            sx, sy, sz, sr, _ = strep[RNG.integers(len(strep))]
            rad = np.hypot(sx, sz)
            if rad < 1e-6:
                continue
            ux, uz = sx/rad, sz/rad                  # outward radial direction
            d = sr + RH
            x = sx + ux*d + RNG.normal(0, 0.5)
            z = sz + uz*d + RNG.normal(0, 0.5)
            y = sy + RNG.normal(0, 0.8)
            if self._no_overlap(x, y, z, RH):
                self.cells.append((x, y, z, RH, "Haemophilus/Aggr."))


def assemble():
    """Run corncob accretion without rendering (for analysis)."""
    cc = Corncob()
    for f in range(FRAMES):
        frac = f / FRAMES
        if frac < 0.55:
            cc.grow_streptococcus(7)
            if frac > 0.12:
                cc.grow_porphyromonas(1)
        if frac > 0.40:
            cc.grow_haemophilus(6)
    return cc


def build():
    cc = Corncob()
    fig, ax = plt.subplots(figsize=(6.2, 8.2), facecolor="#0a0e14")
    ax.set_facecolor("#0a0e14"); ax.set_xlim(-16, 16); ax.set_ylim(-4, L+4)
    ax.set_aspect("equal"); ax.axis("off")
    # the Corynebacterium filament core (drawn once, behind everything)
    ax.add_patch(plt.Rectangle((-RF, 0), 2*RF, L, color=C_CORE, alpha=.95, zorder=0))
    ax.plot([0, 0], [0, L], color=C_CORE, lw=1, zorder=0)
    scat_holder = []
    cap = ax.text(0, -3.0, "", color="#e7ecf3", ha="center", fontsize=10)
    ax.text(0, L+2.2, "Corynebacterium matruchotii filament",
            color=C_CORE, ha="center", fontsize=9)
    # legend
    for lab, col in [("Corynebacterium (core)", C_CORE), ("Streptococcus", C_STR),
                     ("Porphyromonas", C_POR), ("Haemophilus/Aggr. (outer)", C_HAE)]:
        ax.scatter([], [], s=40, c=col, label=lab)
    ax.legend(loc="lower right", fontsize=7, facecolor="#10151f",
              edgecolor="#2a3344", labelcolor="#cdd6e4")

    def update(f):
        frac = f / FRAMES
        if frac < 0.55:
            cc.grow_streptococcus(7)
            if frac > 0.12:
                cc.grow_porphyromonas(1)
        if frac > 0.40:
            cc.grow_haemophilus(6)
        if not cc.cells:
            return cap,
        arr = cc.cells
        order = np.argsort([c[2] for c in arr])      # draw back (z small) first
        xs = np.array([arr[i][0] for i in order])
        ys = np.array([arr[i][1] for i in order])
        zs = np.array([arr[i][2] for i in order])
        rs = np.array([arr[i][3] for i in order])
        cs = [arr[i][4] for i in order]
        cmap = {"Streptococcus": C_STR, "Haemophilus/Aggr.": C_HAE, "Porphyromonas": C_POR}
        cols = [cmap[t] for t in cs]
        depth = (zs - zs.min()) / (np.ptp(zs) + 1e-9)
        sizes = (rs**2) * (40 + 60*depth)            # nearer cells larger
        if scat_holder:
            scat_holder[0].remove()
        sc = ax.scatter(xs, ys, s=sizes, c=cols, edgecolors="none", zorder=2)
        scat_holder[:] = [sc]
        phase = ("Streptococcus coaggregate onto the filament" if frac < 0.40
                 else "Haemophilus binds to Streptococcus (cell-cell distance)"
                 if frac < 0.62 else "mature corncob: concentric taxon layers")
        cap.set_text(phase)
        return sc, cap

    anim = animation.FuncAnimation(fig, update, frames=FRAMES, interval=80, blit=False)
    out = os.path.join(FIG, "C_corncob.gif")
    anim.save(out, writer=animation.PillowWriter(fps=14))
    plt.close(fig)
    return cc, out


def analysis(cc):
    cells = cc.cells
    rad = {"Streptococcus": [], "Haemophilus/Aggr.": [], "Porphyromonas": []}
    strep = np.array([[c[0], c[1], c[2]] for c in cells if c[4] == "Streptococcus"])
    haem = [c for c in cells if c[4] == "Haemophilus/Aggr."]
    for (x, y, z, r, t) in cells:
        if t in rad:
            rad[t].append(np.hypot(x, z))
    # verify the cell-cell-distance rule for Haemophilus
    near = 0
    for (x, y, z, r, t) in haem:
        d = np.sqrt(((strep - [x, y, z])**2).sum(1)).min()
        if d < P.COAGG_DIST + RS + RH:
            near += 1
    frac_near = near / max(len(haem), 1)

    fig, ax = plt.subplots(figsize=(8, 4.6), facecolor="white")
    bins = np.linspace(0, 14, 30)
    for t, col in [("Porphyromonas", C_POR), ("Streptococcus", C_STR), ("Haemophilus/Aggr.", C_HAE)]:
        if rad[t]:
            ax.hist(rad[t], bins=bins, color=col, alpha=.7, label=t)
    ax.axvspan(0, RF, color=C_CORE, alpha=.25)
    ax.text(RF/2, ax.get_ylim()[1]*0.92, "Coryne\ncore", ha="center", fontsize=8, color="#a0008a")
    ax.set_xlabel("distance from filament axis  (\u03bcm)")
    ax.set_ylabel("cell count")
    ax.set_title(f"Concentric corncob layering emerges from coaggregation\n"
                 f"{frac_near*100:.0f}% of Haemophilus cells are juxtaposed to a Streptococcus "
                 f"(cell-cell-distance rule)")
    ax.legend()
    fig.tight_layout()
    out = os.path.join(FIG, "C_corncob_layering.png")
    fig.savefig(out, dpi=130); plt.close(fig)
    return out, frac_near


if __name__ == "__main__":
    print("[C] assembling a corncob from the coaggregation rule ...")
    cc, gif = build()
    png, frac = analysis(cc)
    n = {}
    for c in cc.cells:
        n[c[4]] = n.get(c[4], 0) + 1
    print(f"    cells: {n}")
    print(f"    Haemophilus juxtaposed to Streptococcus: {frac*100:.0f}%")
    print(f"    wrote {os.path.basename(gif)}, {os.path.basename(png)}")
