# The hedgehog consortium: ecology and spatial organization

A reading of the radial multi-genus structure of supragingival dental plaque, with three
computational models that recover its organization from local rules.

Sources for the biology:

1. Mark Welch JL, Rossetti BJ, Rieken CW, Dewhirst FE, Borisy GG. *Biogeography of a human
   oral microbiome at the micron scale.* PNAS 2016, 113(6):E791 to E800. PMID 26811460.
2. Kolenbrander PE, Palmer RJ Jr, Periasamy S, Jakubovics NS. *Oral multispecies biofilm
   development and the key role of cell to cell distance.* Nat Rev Microbiol 2010, 8:471 to 480.
   PMID 20514044.

---

## 1. What the structure is

Supragingival plaque is one of the best characterized multispecies biofilms, and Kolenbrander
and colleagues treat it as the paradigm case for how mixed bacterial communities assemble on a
surface bathed in flowing host secretions [2]. Mark Welch and colleagues imaged that community
at micron resolution using spectral fluorescence in situ hybridization, guided by metagenomics,
and found that it is not a homogeneous smear of taxa [1]. Instead, much of the plaque is built
from a recurring, radially arranged consortium that they named the hedgehog.

The hedgehog is organized around filaments of the genus *Corynebacterium*, principally
*Corynebacterium matruchotii*. These filaments form a structural backbone that radiates outward
from a base anchored to the existing tooth biofilm. Around that backbone, eight other taxa occupy
distinct positions: *Streptococcus*, *Haemophilus* / *Aggregatibacter*, *Porphyromonas*,
Neisseriaceae, *Capnocytophaga*, *Fusobacterium*, *Leptotrichia*, and *Actinomyces*, with sparse
additional rods and cocci grouped as "other." A single hedgehog ranges from a few tens to a few
hundreds of microns in radius. The two-dimensional fan diagram that this repository starts from is
best understood as a wedge cut through this radially symmetric object.

## 2. The organizing principle: a gradient and a set of niches

The central ecological observation is that the taxa are not placed at random. Their positions
track function. Welch and colleagues report two clear patterns [1]:

- **Redox stratification.** Anaerobic taxa sit toward the interior, while facultative and obligate
  aerobes occupy the periphery. The structure is, in effect, a small oxygen reactor. Oxygen enters
  from the saliva-exposed rim and is consumed faster than it can diffuse inward, so the core is
  anoxic and the rim is oxic.
- **Metabolic proximity.** Producers and consumers of shared metabolites are placed close together.
  Streptococci at the periphery ferment dietary sugars and release lactate, acetate, carbon dioxide,
  and hydrogen peroxide, generating a local microenvironment that is acidic, carbon dioxide rich,
  and low in oxygen. Just inside the corncob shell this produces an annulus that suits the
  carbon dioxide requiring *Capnocytophaga* and the microaerophilic filaments *Fusobacterium* and
  *Leptotrichia*.

The peripheral motif is the corncob: a *Corynebacterium* filament tip sheathed in cocci. Welch and
colleagues describe streptococci and *Porphyromonas* in direct contact with the filament, and
*Haemophilus* / *Aggregatibacter* in contact with the streptococci rather than the filament [1].
Clusters of Neisseriaceae also sit at the periphery. The base of the structure is dominated by
*Corynebacterium* and only thinly populated by other cells, with *Actinomyces* appearing in patches
in and around the base.

The colonization logic behind these placements comes from Kolenbrander and colleagues [2]. Three
of their points matter here:

- **Adhesion is mandatory.** Planktonic cells are swallowed before they can grow, so persistence
  requires attachment to a surface or to other cells.
- **Coaggregation builds the community.** Specific, partner to partner recognition lets later
  arrivals dock onto cells that are already present. The set of who recognizes whom shapes the
  order in which a community is built and where each cell can sit.
- **Cell to cell distance is decisive.** Metabolic cooperation and signalling fall off sharply with
  separation, so the precise micron scale arrangement of partners is not incidental, it is the thing
  that makes the community work.

Welch and colleagues argue that the hedgehog is therefore an emergent structure: its global order is
not encoded anywhere, it falls out of local adhesion, coaggregation, metabolism, and the oxygen
gradient acting at the scale of single cells [1]. The three models in this repository are built to
test that claim by encoding only local rules and asking whether the global structure appears.

## 3. Model A: the gradient that sets the zones

`src/reaction_diffusion.py` solves the steady chemistry on a radial coordinate, treating the
hedgehog as radially symmetric. Three coupled reaction diffusion equations are integrated to steady
state: oxygen with Michaelis Menten respiration that is weighted toward the dense, active peripheral
shell; lactate produced by peripheral fermentation and consumed in the interior; and carbon dioxide
produced peripherally and slowly washed at the rim.

The relevant dimensionless number is the Thiele modulus, which compares the rate of consumption to
the rate of diffusive supply:

    phi = R * sqrt( Qmax / ( D_O2 * (Km + O2_rim) ) )

With order of magnitude biofilm parameters and a radius of 150 microns, the model gives phi of about
5.5 and an oxygen penetration depth of about 38 microns. Because phi is well above one, the system is
reaction limited: oxygen is consumed in a thin peripheral rind and never reaches the centre. The
measured anoxic core extends to about 66 microns. This single result reproduces the redox
stratification reported in the imaging: a thin oxic rim over a large anoxic interior. The lactate
profile peaks just inside the rim, where the streptococci that make it sit, and decays inward as it
is consumed, which is the spatial signature of producers and consumers placed next to each other.
Carbon dioxide accumulates across the sub-rim annulus, defining the niche that the capnophilic and
microaerophilic taxa occupy.

Outputs: `figures/A_o2_gradient.gif`, `figures/A_radial_profiles.png`.

## 4. Model B: the structure assembles itself

`src/hedgehog_assembly.py` is an agent based model that encodes only local rules and then watches the
nine taxon architecture appear. The developmental program follows the colonization sequence of
Kolenbrander and colleagues and the summary hypothesis of Welch and colleagues:

1. An early colonizer lawn of *Streptococcus* and *Actinomyces* forms at the base.
2. *Corynebacterium* filaments nucleate at the base and elongate radially outward.
3. New cells attach only where a partner taxon already sits within a contact distance, the cell to
   cell distance rule. Each candidate site is also checked against a local oxygen value equal to the
   radial profile attenuated by local crowding, so the dense corncob shell shields the cells just
   beneath it and creates the low oxygen microniche that anaerobes and microaerophiles require.

No taxon is ever told its target radius. The radial zonation in `figures/B_radial_zonation.png` is
what emerges: a base population of *Actinomyces*, "other," and early *Streptococcus*; a mid to outer
annulus of *Fusobacterium*, *Leptotrichia*, and *Capnocytophaga*; and a peripheral corncob shell of
*Streptococcus*, *Haemophilus*, *Porphyromonas*, and Neisseriaceae. The mechanism that produces the
anaerobic microniche at the oxic periphery is the same crowding based shielding that lets
*Porphyromonas*, an anaerobe, live inside a corncob in an otherwise oxygenated zone.

Outputs: `figures/B_assembly.gif`, `figures/B_radial_zonation.png`.

## 5. Model C: the corncob and the cell to cell distance rule

`src/corncob_coaggregation.py` zooms in on a single filament tip and makes the coaggregation rule
explicit. Streptococci dock directly onto the *Corynebacterium* filament. *Porphyromonas* also docks
onto the filament and is interspersed with the streptococci. *Haemophilus* can only attach where a
*Streptococcus* already sits within a contact distance, never directly to the filament. Steric
exclusion prevents cells from overlapping.

From that one rule a concentric, layered corncob appears, shown in `figures/C_corncob_layering.png`:
a *Corynebacterium* core, a *Streptococcus* and *Porphyromonas* shell in direct contact with it, and
a *Haemophilus* outer layer. The model also verifies the rule quantitatively: about 99 percent of the
*Haemophilus* cells are juxtaposed to a *Streptococcus*, exactly the contact dependent arrangement
that Welch and colleagues observed in the images and that Kolenbrander and colleagues identified as
the organizing constraint.

Outputs: `figures/C_corncob.gif`, `figures/C_corncob_layering.png`.

## 6. What the three models show together

The models are deliberately layered. Model A explains why the zones exist: a reaction limited oxygen
gradient with a high Thiele modulus produces a thin oxic rim and an anoxic core, and the
fermentation byproducts of the peripheral cocci set up the chemical microenvironment of the annulus.
Model B shows that the whole nine taxon structure can self organize from adhesion, coaggregation, and
that local chemistry, without any global blueprint. Model C isolates the single rule, contact
dependent attachment, that builds the peripheral corncob and explains the contact arrangement of its
layers. Read together they support the emergence thesis: a complex, reproducible spatial structure
arising from interactions at the scale of single cells.

## 7. Caveats

These are minimal, didactic models, not fitted simulations. The parameters are order of magnitude
values chosen to put the system in the correct regime, and the geometry is idealized as radially
symmetric. The crowding based oxygen shielding in Model B is a simple proxy for the true coupling
between local respiration and local oxygen, which Model A treats more carefully but only in one
dimension. The taxa traits in `src/params.py` are simplifications of the niches described in the
literature. The aim is to show that the reported organization is consistent with, and recoverable
from, a small set of physically reasonable local rules, not to claim a unique quantitative fit.

## 8. References

1. Mark Welch JL, Rossetti BJ, Rieken CW, Dewhirst FE, Borisy GG. Biogeography of a human oral
   microbiome at the micron scale. Proc Natl Acad Sci USA. 2016;113(6):E791 to E800.
   doi:10.1073/pnas.1522149113. PMID 26811460.
2. Kolenbrander PE, Palmer RJ Jr, Periasamy S, Jakubovics NS. Oral multispecies biofilm development
   and the key role of cell to cell distance. Nat Rev Microbiol. 2010;8:471 to 480.
   doi:10.1038/nrmicro2381. PMID 20514044.
