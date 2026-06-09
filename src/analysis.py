"""
analysis.py  --  one dashboard with the quantitative analysis from every model
===============================================================================
Pulls the key quantitative result from each model into a single 2x3 figure:

  A  radial chemistry (O2 / lactate / CO2) with Thiele modulus + penetration depth
  B  emergent radial zonation (taxon density vs radius)            [Model B]
  C  concentric corncob layering + cell-cell-distance check        [Model C]
  D  azimuthal heterogeneity: coupled 3D vs radially symmetric     [Model E]
  E  oxygen depletion over colonisation time                       [Model F]
  F  nutrient and oxygen succession (mean fields vs time)          [Model F]

Run after the models (run_all.py does this); Model F metrics are read from
figures/F_metrics.npz, or recomputed if that file is absent.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import params as PA
import reaction_diffusion as A
import hedgehog_assembly as B
import corncob_coaggregation as C
import reaction_diffusion_3d as E

FIG = os.path.join(os.path.dirname(__file__), "..", "figures")


def panel_A(ax):
    snaps = A.simulate(); _, O2, Lac, CO2 = snaps[-1]; r = A.r; R = A.R
    phi = R*np.sqrt(A.P.QMAX_O2/(A.P.D_O2*(A.P.KM_O2+A.P.O2_SAT)))
    delta = np.sqrt(2*A.P.D_O2*A.P.O2_SAT/A.P.QMAX_O2)
    anox = r[O2 < A.P.ANOXIC_THR]; core = anox.max() if anox.size else 0
    ax.plot(r, O2, color="#ff8a3d", lw=2.3, label="O$_2$")
    ax.plot(r, Lac/Lac.max()*0.9, color="#35d24b", lw=1.8, label="lactate")
    ax.plot(r, CO2/CO2.max()*0.9, color="#2ad4e6", lw=1.8, label="CO$_2$")
    ax.axvspan(0, core, color="#2b3bff", alpha=0.08)
    ax.set_title(f"A · radial chemistry\n\u03c6={phi:.1f}, \u03b4\u2248{delta:.0f}\u03bcm, anoxic r<{core:.0f}\u03bcm", fontsize=10)
    ax.set_xlabel("radius r (\u03bcm)"); ax.set_ylabel("norm. conc."); ax.legend(fontsize=7); ax.set_xlim(0, R)


def panel_B(ax):
    h = B.assemble()
    cx = np.array(h.cx); cy = np.array(h.cy); ct = np.array(h.ct, dtype=object)
    rn = np.hypot(cx, cy)/B.R
    bins = np.linspace(0, 1, 26); ctr = 0.5*(bins[1:]+bins[:-1])
    area = np.pi*(bins[1:]**2 - bins[:-1]**2)
    order = ["Actinomyces", "other", "Leptotrichia", "Fusobacterium", "Capnocytophaga",
             "Porphyromonas", "Neisseriaceae", "Streptococcus", "Haemophilus/Aggr."]
    for t in order:
        m = ct == t
        if m.any():
            hi, _ = np.histogram(rn[m], bins=bins)
            ax.plot(ctr, hi/area, color=PA.TAXA[t]["color"], lw=1.6, label=t)
    ax.set_title("B · emergent radial zonation", fontsize=10)
    ax.set_xlabel("r / R  (base \u2192 rim)"); ax.set_ylabel("areal density (a.u.)")
    ax.legend(fontsize=5.5, ncol=2); ax.set_xlim(0, 1)


def panel_C(ax):
    cc = C.assemble(); cells = cc.cells
    strep = np.array([[c[0], c[1], c[2]] for c in cells if c[4] == "Streptococcus"])
    haem = [c for c in cells if c[4] == "Haemophilus/Aggr."]
    rad = {"Porphyromonas": [], "Streptococcus": [], "Haemophilus/Aggr.": []}
    for (x, y, z, rr, t) in cells:
        if t in rad:
            rad[t].append(np.hypot(x, z))
    near = sum(1 for (x, y, z, rr, t) in haem
               if np.sqrt(((strep-[x, y, z])**2).sum(1)).min() < PA.COAGG_DIST+C.RS+C.RH)
    frac = near/max(len(haem), 1)
    bins = np.linspace(0, 14, 30)
    for t, col in [("Porphyromonas", "#2233e0"), ("Streptococcus", "#35d24b"), ("Haemophilus/Aggr.", "#ff9a1f")]:
        if rad[t]:
            ax.hist(rad[t], bins=bins, color=col, alpha=.7, label=t)
    ax.axvspan(0, C.RF, color="#ff36c8", alpha=.25)
    ax.set_title(f"C · corncob layering\n{frac*100:.0f}% of Haemophilus juxtaposed to Streptococcus", fontsize=10)
    ax.set_xlabel("distance from filament axis (\u03bcm)"); ax.set_ylabel("cells"); ax.legend(fontsize=6.5)


def panel_D(ax):
    cells, fils = E.build_biomass()
    Wr = E.deposit(cells, E.RESP_W); Wf = E.deposit(cells, E.FERM_W)
    O, Lf, Cf = E.solve(Wr, Wf)
    izm = int(0.5*E.R/E.dx); r_ring = 0.78*E.R
    th = np.linspace(0, 2*np.pi, 240)
    ix = np.clip(((r_ring*np.cos(th)+E.R)/E.dx).round().astype(int), 0, E.NX-1)
    iy = np.clip(((r_ring*np.sin(th)+E.R)/E.dx).round().astype(int), 0, E.NY-1)
    vals = O[ix, iy, izm]; cv = vals.std()/max(vals.mean(), 1e-9)
    ax.plot(np.degrees(th), vals, color="#ff8a3d", lw=2, label="coupled 3D solve")
    ax.axhline(vals.mean(), color="#888", ls="--", lw=1.4, label="symmetric model")
    ax.set_title(f"D · azimuthal heterogeneity\nO$_2$ on a ring, CV = {cv*100:.0f}%", fontsize=10)
    ax.set_xlabel("azimuth (degrees)"); ax.set_ylabel("O$_2$ fraction"); ax.legend(fontsize=7); ax.set_xlim(0, 360)


def load_metrics():
    p = os.path.join(FIG, "F_metrics.npz")
    if os.path.exists(p):
        d = np.load(p); return {k: d[k] for k in d.files}
    import succession_3d as F
    (_, _, _), _, _, m = F.run_succession(); return m


def panel_E(ax, mt):
    axb = ax.twinx()
    ax.plot(mt["tf"], mt["anox"]*100, color="#2233e0", lw=2.3)
    axb.plot(mt["tf"], mt["ncells"], color="#8a93a0", lw=1.6, ls="--")
    ax.set_title("E · oxygen depletes as colonisation proceeds", fontsize=10)
    ax.set_xlabel("time t (norm.)"); ax.set_ylabel("anoxic volume (%)", color="#2233e0")
    axb.set_ylabel("visible cells", color="#777"); ax.tick_params(axis="y", colors="#2233e0")


def panel_F(ax, mt):
    nrm = lambda x: x/max(x.max(), 1e-9)
    ax.plot(mt["tf"], mt["meanO2"], color="#ff8a3d", lw=2.3, label="mean O$_2$")
    ax.plot(mt["tf"], nrm(mt["meanCO2"]), color="#2ad4e6", lw=1.8, label="CO$_2$ (norm.)")
    ax.plot(mt["tf"], nrm(mt["meanLac"]), color="#35d24b", lw=1.8, label="lactate (norm.)")
    ax.set_title("F · nutrient and oxygen succession", fontsize=10)
    ax.set_xlabel("time t (norm.)"); ax.set_ylabel("mean concentration"); ax.legend(fontsize=7)


def main():
    print("[dashboard] assembling consolidated analysis ...")
    mt = load_metrics()
    fig, axes = plt.subplots(2, 3, figsize=(16, 9), facecolor="white")
    print("    panel A (chemistry) ..."); panel_A(axes[0, 0])
    print("    panel B (zonation) ...");  panel_B(axes[0, 1])
    print("    panel C (corncob) ...");   panel_C(axes[0, 2])
    print("    panel D (heterogeneity, 3D solve) ..."); panel_D(axes[1, 0])
    panel_E(axes[1, 1], mt)
    panel_F(axes[1, 2], mt)
    fig.suptitle("Hedgehog consortium: consolidated analysis across all models",
                 fontsize=16, y=0.98)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = os.path.join(FIG, "analysis_dashboard.png")
    fig.savefig(out, dpi=140); plt.close(fig)
    print(f"    wrote {os.path.basename(out)}")
    return out


if __name__ == "__main__":
    main()
