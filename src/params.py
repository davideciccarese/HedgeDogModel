"""
params.py
=========
Shared definitions for the hedgehog-consortium models.

Taxa, colours and ecological traits are taken from the radial nine-taxon model of
supragingival plaque described by Mark Welch, Rossetti, Rieken, Dewhirst & Borisy
(PNAS 2016, 113:E791-E800; PMID 26811460, "Biogeography of a human oral microbiome
at the micron scale"), and the colonisation / coaggregation framework of Kolenbrander,
Palmer, Periasamy & Jakubovics (Nat Rev Microbiol 2010, 8:471-480; PMID 20514044,
"Oral multispecies biofilm development and the key role of cell-cell distance").

Colours match the source schematic so all three models and the 3D viewer share one
visual language.
"""

# ----------------------------------------------------------------------
# Taxa.  o2_pref is the centre of each taxon's preferred dissolved-O2 band
# (1 = fully oxygenated perimeter, 0 = anoxic core).  r_band is the preferred
# radial position as a fraction of hedgehog radius R (0 = base/tooth, 1 = rim).
# morph: 'filament' (rod chains forming the scaffold/annulus) or 'coccus'/'rod'.
# partners: taxa it coaggregates with (drives the cell-cell-distance assembly rule).
# ----------------------------------------------------------------------
TAXA = {
    "Corynebacterium":       dict(color="#ff36c8", o2_pref=0.55, r_band=0.50, morph="filament",
                                  role="structural backbone; radiates from the base",
                                  partners=["Streptococcus", "Porphyromonas"]),
    "Streptococcus":         dict(color="#35d24b", o2_pref=0.85, r_band=0.95, morph="coccus",
                                  role="corncob cocci; ferments sugars -> lactate, acetate, CO2, H2O2",
                                  partners=["Corynebacterium", "Haemophilus", "Actinomyces"]),
    "Haemophilus/Aggr.":     dict(color="#ff9a1f", o2_pref=0.90, r_band=1.00, morph="coccus",
                                  role="outermost corncob layer; binds to Streptococcus",
                                  partners=["Streptococcus"]),
    "Porphyromonas":         dict(color="#2233e0", o2_pref=0.30, r_band=0.92, morph="coccus",
                                  role="anaerobe in direct contact with Corynebacterium in corncobs",
                                  partners=["Corynebacterium"]),
    "Neisseriaceae":         dict(color="#8a2bbf", o2_pref=0.80, r_band=0.90, morph="coccus",
                                  role="aerobic clusters at the periphery",
                                  partners=["Streptococcus"]),
    "Capnocytophaga":        dict(color="#e61f1f", o2_pref=0.45, r_band=0.70, morph="rod",
                                  role="CO2-requiring (capnophilic) rods in/around the annulus",
                                  partners=["Fusobacterium"]),
    "Fusobacterium":         dict(color="#f4e016", o2_pref=0.35, r_band=0.68, morph="filament",
                                  role="bridging filament; thrives in low-O2 / high-CO2 annulus",
                                  partners=["Corynebacterium", "Streptococcus", "Leptotrichia"]),
    "Leptotrichia":          dict(color="#2ad4e6", o2_pref=0.40, r_band=0.66, morph="filament",
                                  role="mid-zone filament intermingled with Fusobacterium",
                                  partners=["Corynebacterium", "Fusobacterium"]),
    "Actinomyces":           dict(color="#eef0f3", o2_pref=0.50, r_band=0.18, morph="rod",
                                  role="early coloniser; patches in the base region",
                                  partners=["Streptococcus", "Corynebacterium"]),
    "other":                 dict(color="#8a93a0", o2_pref=0.40, r_band=0.20, morph="coccus",
                                  role="sparse rods/cocci thinly populating the base",
                                  partners=[]),
}

# ----------------------------------------------------------------------
# Physical parameters (micron / second units).  Order-of-magnitude values
# for a mature supragingival biofilm; the hedgehog itself is tens to a few
# hundred microns in radius (Welch et al. 2016).
# ----------------------------------------------------------------------
R_HEDGEHOG   = 150.0     # consortium radius, um
D_O2         = 1.8e3     # O2 diffusivity in biofilm, um^2/s (~0.6x bulk water)
D_LAC        = 1.0e3     # lactate diffusivity, um^2/s
D_CO2        = 1.9e3     # CO2 diffusivity, um^2/s
O2_SAT       = 1.0       # dissolved O2 at the saliva-exposed rim (normalised)
QMAX_O2      = 2.5      # max areal O2 uptake (normalised units / s)
KM_O2        = 0.05      # half-saturation for O2 uptake (normalised)
LAC_PROD     = 2.0e-2    # peripheral lactate production rate
CO2_PROD     = 3.0e-2    # peripheral CO2 production rate
ANOXIC_THR   = 0.05      # O2 fraction below which a zone counts as "anoxic"

# Cell-cell distance for coaggregation-mediated attachment, um
# (Kolenbrander et al. 2010 - juxtaposition / direct contact is required).
COAGG_DIST   = 1.5

def hex_to_rgb(h):
    h = h.lstrip("#")
    return tuple(int(h[i:i+2], 16) / 255 for i in (0, 2, 4))
