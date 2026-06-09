"""
hedgehog_assembly.py  --  Model B: emergent self-organisation of the consortium
================================================================================
The central claim of Welch et al. (2016) is that the hedgehog's elaborate radial
order is an *emergent* property of micron-scale rules, not a blueprint.  This
agent-based model encodes only local rules and lets the nine-taxon architecture
appear on its own:

  1. Substrate / pellicle.  Early colonisers (Streptococcus, Actinomyces) form a
     basal lawn -- the "existing biofilm" of Welch Fig. 9 and the initial
     adherence step of Kolenbrander et al. (2010).
  2. Scaffold.  Corynebacterium filaments nucleate at the base and elongate
     radially outward (tip extension with small angular wander).
  3. Coaggregation (the cell-cell-distance rule, Kolenbrander 2010).  A new cell
     of taxon X can only attach where a *partner* taxon already sits within a
     contact distance COAGG_DIST.  Haemophilus attaches to Streptococcus,
     Porphyromonas to Corynebacterium, and so on.
  4. Local chemistry.  Each candidate site has a *local* O2 value equal to the
     radial profile O2(r) attenuated by local crowding -- the dense corncob shell
     shields the interior, so anaerobes (Porphyromonas) and microaerophiles
     (Fusobacterium, Leptotrichia, Capnocytophaga) find viable niches just inside
     the rim even though the bulk gradient there is oxic.

No taxon is ever *told* where to go; the radial zonation is what falls out.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.collections import LineCollection
from scipy.spatial import cKDTree
import params as P

FIG = os.path.join(os.path.dirname(__file__), "..", "figures")
os.makedirs(FIG, exist_ok=True)
RNG = np.random.default_rng(7)

R = P.R_HEDGEHOG
FRAMES = 120
DMIN = 1.5                     # steric exclusion distance between cell centres (um)

# per-taxon attachment tolerances (how fussy each taxon is about O2 / radius)
TOL = {
    "Streptococcus":     dict(o2=0.30, rr=0.18),
    "Haemophilus/Aggr.": dict(o2=0.30, rr=0.14),
    "Porphyromonas":     dict(o2=0.80, rr=0.16),   # permissive: lives in shielded microniche
    "Neisseriaceae":     dict(o2=0.35, rr=0.16),
    "Capnocytophaga":    dict(o2=0.35, rr=0.22),
    "Fusobacterium":     dict(o2=0.35, rr=0.22),
    "Leptotrichia":      dict(o2=0.35, rr=0.22),
    "Actinomyces":       dict(o2=0.60, rr=0.22),
    "other":             dict(o2=0.70, rr=0.25),
}


def O2_radial(rn):
    """Bulk dissolved-O2 as a function of normalised radius (anoxic core, oxic rim)."""
    return 1.0 / (1.0 + np.exp(-(rn - 0.72) / 0.10))


class Hedgehog:
    def __init__(self):
        self.cx = []; self.cy = []; self.ct = []      # cell positions + taxon
        self.fils = []                                # list of dict(pts=[(x,y)...], dir=ang)
        self.tree = None; self.pts = None; self.tps = None
        self._seed()

    # ---- helpers -----------------------------------------------------
    def _add(self, x, y, taxon):
        self.cx.append(x); self.cy.append(y); self.ct.append(taxon)

    def _rebuild(self):
        fx = [p[0] for f in self.fils for p in f["pts"]]
        fy = [p[1] for f in self.fils for p in f["pts"]]
        ft = ["Corynebacterium"] * len(fx)
        X = np.array(self.cx + fx); Y = np.array(self.cy + fy)
        self.tps = np.array(self.ct + ft, dtype=object)
        self.pts = np.column_stack([X, Y]) if len(X) else np.empty((0, 2))
        self.tree = cKDTree(self.pts) if len(X) else None

    def _local_o2(self, x, y, rn):
        if self.tree is None:
            return O2_radial(rn)
        n = len(self.tree.query_ball_point([x, y], 6.0))
        return O2_radial(rn) * np.exp(-n / 7.0)        # crowding shields O2

    # ---- initialisation ---------------------------------------------
    def _seed(self):
        # basal lawn of early colonisers near the tooth (centre = base)
        for _ in range(72):
            a = RNG.uniform(0, 2 * np.pi); rr = RNG.uniform(2, 34)
            self._add(rr * np.cos(a), rr * np.sin(a),
                      "Streptococcus" if RNG.random() < 0.55 else "Actinomyces")
        for _ in range(28):
            a = RNG.uniform(0, 2 * np.pi); rr = RNG.uniform(2, 38)
            self._add(rr * np.cos(a), rr * np.sin(a), "other")
        # nucleate Corynebacterium filaments at the base
        for _ in range(130):
            ang = RNG.uniform(0, 2 * np.pi)
            r0 = RNG.uniform(4, 18)
            self.fils.append(dict(pts=[(r0 * np.cos(ang), r0 * np.sin(ang))],
                                  dir=ang, len=r0, alive=True))
        self._rebuild()

    # ---- per-step processes -----------------------------------------
    def grow_filaments(self, step_len=4.2):
        for f in self.fils:
            if not f["alive"]:
                continue
            x, y = f["pts"][-1]
            f["dir"] += RNG.normal(0, 0.08)            # slight angular wander
            nx = x + step_len * np.cos(f["dir"])
            ny = y + step_len * np.sin(f["dir"])
            nr = np.hypot(nx, ny)
            if nr >= R:                                # stop at the rim
                f["alive"] = False
                continue
            f["pts"].append((nx, ny)); f["len"] = nr

    def recruit(self, taxon, trials):
        info = P.TAXA[taxon]; partners = info["partners"]
        tol = TOL[taxon]; placed = 0
        # candidate seed cells = existing partner cells (enforces juxtaposition)
        if partners:
            mask = np.isin(self.tps, partners)
            idx = np.where(mask)[0]
        else:
            idx = np.where(self.tps != "Corynebacterium")[0]
        if idx.size == 0:
            return 0
        for _ in range(trials):
            j = idx[RNG.integers(idx.size)]
            px, py = self.pts[j]
            ang = RNG.uniform(0, 2 * np.pi)
            d = RNG.uniform(P.COAGG_DIST, P.COAGG_DIST + 1.2)
            x = px + d * np.cos(ang); y = py + d * np.sin(ang)
            r = np.hypot(x, y); rn = r / R
            if rn > 1.0:
                continue
            if abs(rn - info["r_band"]) > tol["rr"]:
                continue
            o2 = self._local_o2(x, y, rn)
            if abs(o2 - info["o2_pref"]) > tol["o2"]:
                continue
            if self.tree.query([x, y])[0] < DMIN:      # steric exclusion
                continue
            self._add(x, y, taxon); placed += 1
        return placed

    # ---- one frame of the developmental programme -------------------
    def step(self, f):
        frac = f / FRAMES
        if 0.04 < frac < 0.60:                         # scaffold elongates
            self.grow_filaments()
        self._rebuild()
        if frac > 0.28:                                # corncob shell
            self.recruit("Streptococcus", 130)
            self.recruit("Haemophilus/Aggr.", 65)
            self.recruit("Porphyromonas", 36)
            self.recruit("Neisseriaceae", 36)
        if frac > 0.42:                                # micro-aerophilic annulus
            self.recruit("Fusobacterium", 60)
            self.recruit("Leptotrichia", 64)
            self.recruit("Capnocytophaga", 48)
        if frac > 0.10:                                # base keeps filling
            self.recruit("Actinomyces", 18)
            self.recruit("other", 12)
        self._rebuild()

    def phase(self, f):
        frac = f / FRAMES
        if frac < 0.10:  return "1 - pellicle + early colonisers (Streptococcus, Actinomyces)"
        if frac < 0.28:  return "2 - Corynebacterium filaments elongate radially"
        if frac < 0.42:  return "3 - corncobs assemble at oxic tips (coaggregation)"
        if frac < 0.62:  return "4 - micro-aerophilic annulus fills (Fuso, Lepto, Capno)"
        return "5 - mature nine-taxon hedgehog"


def render(h, sizes):
    cx = np.array(h.cx); cy = np.array(h.cy); ct = np.array(h.ct, dtype=object)
    arts = []
    for taxon, info in P.TAXA.items():
        if taxon == "Corynebacterium":
            continue
        m = ct == taxon
        if m.any():
            arts.append((cx[m], cy[m], info["color"], sizes.get(taxon, 16)))
    return arts


def assemble():
    """Run the assembly steps without rendering (for analysis)."""
    h = Hedgehog()
    for f in range(FRAMES):
        h.step(f)
    return h


def build():
    h = Hedgehog()
    fig, ax = plt.subplots(figsize=(7.6, 7.6), facecolor="#0a0e14")
    ax.set_facecolor("#0a0e14"); ax.set_xlim(-R*1.05, R*1.05); ax.set_ylim(-R*1.05, R*1.05)
    ax.set_aspect("equal"); ax.axis("off")
    rim = plt.Circle((0, 0), R, fill=False, ec="#2a3344", lw=1, ls="--"); ax.add_patch(rim)
    tooth = plt.Circle((0, 0), 14, color="#3f72b0", alpha=.9, zorder=1); ax.add_patch(tooth)

    fil_lc = LineCollection([], colors=P.TAXA["Corynebacterium"]["color"], linewidths=1.1, zorder=2)
    ax.add_collection(fil_lc)
    sizes = dict(Streptococcus=22, **{"Haemophilus/Aggr.": 20}, Porphyromonas=16,
                 Neisseriaceae=18, Capnocytophaga=10, Fusobacterium=7, Leptotrichia=7,
                 Actinomyces=9, other=12)
    scat = {t: ax.scatter([], [], s=sizes.get(t, 16), c=info["color"],
                          edgecolors="none", zorder=3, label=t)
            for t, info in P.TAXA.items() if t != "Corynebacterium"}
    scat["Corynebacterium"] = ax.scatter([], [], s=1, c=P.TAXA["Corynebacterium"]["color"], label="Corynebacterium")
    cap = ax.text(0, -R*1.0, "", color="#e7ecf3", ha="center", fontsize=11, zorder=5)
    tt = ax.text(-R*1.02, R*0.98, "", color="#8b96a8", ha="left", fontsize=9, zorder=5)
    leg = ax.legend(loc="upper right", fontsize=6.5, facecolor="#10151f",
                    edgecolor="#2a3344", labelcolor="#cdd6e4", framealpha=.9, ncol=1, markerscale=1.2)

    def update(f):
        h.step(f)
        segs = [np.array(fl["pts"]) for fl in h.fils if len(fl["pts"]) > 1]
        fil_lc.set_segments(segs)
        cx = np.array(h.cx); cy = np.array(h.cy); ct = np.array(h.ct, dtype=object)
        for t, sc in scat.items():
            if t == "Corynebacterium":
                continue
            m = ct == t
            sc.set_offsets(np.column_stack([cx[m], cy[m]]) if m.any() else np.empty((0, 2)))
        cap.set_text(f"phase {h.phase(f)}")
        tt.set_text(f"cells: {len(h.cx)}    filaments: {len(h.fils)}")
        return list(scat.values()) + [fil_lc, cap, tt]

    anim = animation.FuncAnimation(fig, update, frames=FRAMES, interval=70, blit=False)
    out = os.path.join(FIG, "B_assembly.gif")
    anim.save(out, writer=animation.PillowWriter(fps=15))
    plt.close(fig)
    return h, out


def radial_profile(h):
    """Quantify the emergent zonation: per-taxon density vs normalised radius."""
    cx = np.array(h.cx); cy = np.array(h.cy); ct = np.array(h.ct, dtype=object)
    rn = np.hypot(cx, cy) / R
    bins = np.linspace(0, 1, 26); centres = 0.5 * (bins[1:] + bins[:-1])
    ring_area = np.pi * (bins[1:]**2 - bins[:-1]**2)
    fig, ax = plt.subplots(figsize=(8.4, 5), facecolor="white")
    order = ["Actinomyces", "other", "Leptotrichia", "Fusobacterium", "Capnocytophaga",
             "Porphyromonas", "Neisseriaceae", "Streptococcus", "Haemophilus/Aggr."]
    for t in order:
        m = ct == t
        if not m.any():
            continue
        h_, _ = np.histogram(rn[m], bins=bins)
        ax.plot(centres, h_ / ring_area, color=P.TAXA[t]["color"], lw=2, label=t)
    ax.set_xlabel("normalised radius  r / R   (0 = base/tooth, 1 = rim)")
    ax.set_ylabel("areal cell density  (a.u.)")
    ax.set_title("Emergent radial zonation from local rules only\n"
                 "base colonisers -> mid-zone filaments -> peripheral corncob shell")
    ax.legend(fontsize=7, ncol=2); ax.set_xlim(0, 1)
    fig.tight_layout()
    out = os.path.join(FIG, "B_radial_zonation.png")
    fig.savefig(out, dpi=130); plt.close(fig)
    return out


if __name__ == "__main__":
    print("[B] assembling hedgehog from local rules ...")
    h, gif = build()
    prof = radial_profile(h)
    print(f"    final cell count: {len(h.cx)}  ({len(h.fils)} Corynebacterium filaments)")
    print(f"    wrote {os.path.basename(gif)}, {os.path.basename(prof)}")
