"""
reaction_diffusion.py  --  Model A: the chemical gradients behind the zonation
==============================================================================
The hedgehog is, to good approximation, radially symmetric, so the steady
chemistry can be solved on a 1-D radial coordinate and rotated for display.

We solve three coupled reaction-diffusion equations on r in [0, R] (spherical
Laplacian, symmetry at r=0, saliva-exposed rim at r=R):

    dO2/dt  = D_O2 * lap(O2) - w(r) * Qmax * O2 / (Km + O2)        (respiration)
    dLac/dt = D_L  * lap(Lac) + s(r) * Lp        - kL * Lac        (Strep ferment -> sink)
    dCO2/dt = D_C  * lap(CO2) + s(r) * Cp                          (Strep ferment)

w(r) and s(r) are the peripheral activity weights: the metabolically dense
corncob shell of facultative cocci both consumes O2 and excretes acid/CO2, so
both respiration and fermentation peak just inside the rim.  This single fact
reproduces the observed ecology (Welch et al. 2016): O2 is consumed before it
penetrates, leaving an *anoxic interior*, while a CO2-rich, acidified, low-O2
*annulus* forms beneath the corncob shell -- exactly the niche occupied by the
capnophilic, microaerophilic filaments (Fusobacterium, Leptotrichia,
Capnocytophaga).

The relevant dimensionless group is the Thiele modulus
    phi = R * sqrt(Qmax / (D_O2 * (Km + O2_SAT)))
and the penetration depth (zero-order limit) delta ~ sqrt(2 D O2_SAT / Qmax).
phi >> 1  =>  reaction-limited  =>  shallow penetration  =>  anoxic core.
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

N   = 160
R   = P.R_HEDGEHOG
dr  = R / N
r   = (np.arange(N) + 0.5) * dr            # cell centres, avoids r=0 singularity
dt  = 0.20 * dr * dr / max(P.D_O2, P.D_CO2, P.D_LAC)
T_END = 16.0
STEPS = int(T_END / dt)
FRAMES = 90
rec_every = max(1, STEPS // FRAMES)

# peripheral activity weight: dense, active corncob shell near the rim
w = 0.25 + 0.75 * np.exp(-((r - R) / (0.18 * R))**2)     # respiration weight
s = np.exp(-((r - 0.96 * R) / (0.10 * R))**2)            # fermentation source


def lap_sph(C, Crim):
    """Spherical radial Laplacian with Neumann at r=0 and Dirichlet C=Crim at rim."""
    Cext = np.empty(N + 2)
    Cext[1:-1] = C
    Cext[0] = C[0]                       # symmetry: zero flux at centre
    Cext[-1] = 2 * Crim - C[-1]          # ghost node enforcing C(R)=Crim
    d2 = (Cext[2:] - 2 * Cext[1:-1] + Cext[:-2]) / dr**2
    d1 = (Cext[2:] - Cext[:-2]) / (2 * dr)
    return d2 + (2.0 / r) * d1


def simulate():
    O2  = np.full(N, P.O2_SAT)
    Lac = np.zeros(N)
    CO2 = np.zeros(N)
    snaps = []
    for step in range(STEPS):
        upt = w * P.QMAX_O2 * O2 / (P.KM_O2 + O2)
        O2  += dt * (P.D_O2  * lap_sph(O2,  P.O2_SAT) - upt)
        Lac += dt * (P.D_LAC * lap_sph(Lac, 0.0) + s * P.LAC_PROD - 1.0 * Lac)
        CO2 += dt * (P.D_CO2 * lap_sph(CO2, 0.0) + s * P.CO2_PROD - 0.15 * CO2)
        np.clip(O2, 0, None, out=O2)
        if step % rec_every == 0:
            snaps.append((step * dt, O2.copy(), Lac.copy(), CO2.copy()))
    return snaps


def disc(profile):
    """Rotate a radial profile into a 2-D disc image for display."""
    g = np.linspace(-R, R, 220)
    X, Y = np.meshgrid(g, g)
    rr = np.sqrt(X**2 + Y**2)
    img = np.interp(rr, r, profile, right=profile[-1])
    img[rr > R] = np.nan
    return img


def make_gif(snaps):
    fig, (axd, axp) = plt.subplots(1, 2, figsize=(11, 5.2), facecolor="#0a0e14",
                                   gridspec_kw=dict(width_ratios=[1, 1.15]))
    for ax in (axd, axp):
        ax.set_facecolor("#0a0e14")

    im = axd.imshow(disc(snaps[-1][1]), origin="lower", extent=[-R, R, -R, R],
                    cmap="magma", vmin=0, vmax=P.O2_SAT)
    axd.set_title("Dissolved O$_2$  (anoxic core forms within seconds)",
                  color="#e7ecf3", fontsize=11)
    axd.set_xlabel("x  (\u03bcm)", color="#8b96a8"); axd.set_ylabel("y  (\u03bcm)", color="#8b96a8")
    axd.tick_params(colors="#8b96a8")
    cb = fig.colorbar(im, ax=axd, fraction=0.046, pad=0.04)
    cb.set_label("O$_2$  (fraction of saliva value)", color="#8b96a8")
    cb.ax.yaxis.set_tick_params(color="#8b96a8")
    plt.setp(plt.getp(cb.ax.axes, "yticklabels"), color="#8b96a8")

    # zones
    axp.axvspan(0,        0.40 * R, color="#2b3bff", alpha=0.10)
    axp.axvspan(0.40 * R, 0.85 * R, color="#4aa0ff", alpha=0.08)
    axp.axvspan(0.85 * R, R,        color="#ffb24a", alpha=0.10)
    axp.text(0.20 * R, 1.02, "anoxic\nbase", color="#9fb6ff", ha="center", fontsize=8)
    axp.text(0.62 * R, 1.02, "micro-aerophilic\nCO$_2$ annulus", color="#bfe0ff", ha="center", fontsize=8)
    axp.text(0.93 * R, 1.02, "oxic\nrim", color="#ffd089", ha="center", fontsize=8)
    lO, = axp.plot([], [], color="#ff8a3d", lw=2.4, label="O$_2$")
    lL, = axp.plot([], [], color="#35d24b", lw=2.0, label="lactate")
    lC, = axp.plot([], [], color="#2ad4e6", lw=2.0, label="CO$_2$")
    axp.set_xlim(0, R); axp.set_ylim(0, 1.12)
    axp.set_xlabel("radius r  (\u03bcm)", color="#8b96a8")
    axp.set_ylabel("normalised concentration", color="#8b96a8")
    axp.tick_params(colors="#8b96a8")
    for sp in axp.spines.values(): sp.set_color("#2a3344")
    leg = axp.legend(loc="center right", facecolor="#10151f", edgecolor="#2a3344", labelcolor="#e7ecf3")
    ttl = fig.suptitle("", color="#e7ecf3", fontsize=12)

    def update(i):
        t, O2, Lac, CO2 = snaps[i]
        im.set_data(disc(O2))
        lO.set_data(r, O2)
        lL.set_data(r, Lac / max(Lac.max(), 1e-9) * 0.9)
        lC.set_data(r, CO2 / max(CO2.max(), 1e-9) * 0.9)
        ttl.set_text(f"Reaction-diffusion gradients in the hedgehog   t = {t:5.1f} s")
        return im, lO, lL, lC, ttl

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    anim = animation.FuncAnimation(fig, update, frames=len(snaps), interval=70, blit=False)
    out = os.path.join(FIG, "A_o2_gradient.gif")
    anim.save(out, writer=animation.PillowWriter(fps=14))
    plt.close(fig)
    return out


def make_static(snaps):
    """Steady-state profiles + Thiele-modulus / penetration-depth analysis."""
    t, O2, Lac, CO2 = snaps[-1]
    phi = R * np.sqrt(P.QMAX_O2 / (P.D_O2 * (P.KM_O2 + P.O2_SAT)))
    delta = np.sqrt(2 * P.D_O2 * P.O2_SAT / P.QMAX_O2)
    # measured anoxic-core radius
    anoxic = r[O2 < P.ANOXIC_THR]
    core = anoxic.max() if anoxic.size else 0.0

    fig, ax = plt.subplots(figsize=(8, 5), facecolor="white")
    ax.plot(r, O2, color="#ff8a3d", lw=2.5, label="O$_2$")
    ax.plot(r, Lac / Lac.max() * 0.9, color="#35d24b", lw=2, label="lactate (Strep)")
    ax.plot(r, CO2 / CO2.max() * 0.9, color="#2ad4e6", lw=2, label="CO$_2$ (Strep)")
    ax.axhline(P.ANOXIC_THR, color="#888", ls=":", lw=1)
    ax.axvspan(0, core, color="#2b3bff", alpha=0.08)
    ax.set_xlabel("radius r (\u03bcm)"); ax.set_ylabel("normalised concentration")
    ax.set_title(f"Steady-state chemistry\nThiele modulus \u03c6 = {phi:.1f}   |   "
                 f"O$_2$ penetration \u03b4 \u2248 {delta:.0f} \u03bcm   |   anoxic core r < {core:.0f} \u03bcm")
    ax.legend(); ax.set_xlim(0, R); ax.set_ylim(0, 1.05)
    fig.tight_layout()
    out = os.path.join(FIG, "A_radial_profiles.png")
    fig.savefig(out, dpi=130); plt.close(fig)
    return out, phi, delta, core


if __name__ == "__main__":
    print("[A] simulating reaction-diffusion gradients ...")
    snaps = simulate()
    g = make_gif(snaps)
    s, phi, delta, core = make_static(snaps)
    print(f"    Thiele modulus phi      = {phi:.2f}")
    print(f"    O2 penetration depth    = {delta:.1f} um")
    print(f"    measured anoxic core r  < {core:.0f} um  (of R={R:.0f} um)")
    print(f"    wrote {os.path.basename(g)}, {os.path.basename(s)}")
