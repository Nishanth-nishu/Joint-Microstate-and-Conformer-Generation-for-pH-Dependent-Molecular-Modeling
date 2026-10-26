# CoarsenConf (Reidenbach & Krishnapriyan 2023)

## Citation
CoarsenConf: Equivariant Coarsening with Aggregated Attention for Molecular Conformer Generation. arXiv:2306.14852.
Code: https://github.com/ASK-Berkeley/CoarsenConf. Checked 2026-10-26.

## Fulltext status
full_text_read (arXiv v1 HTML; last ~5k chars not read).

## Method
SE(3)-equivariant hierarchical VAE learning p(X | R) where R is an RDKit approximate conformer.
"We coarsen molecules into a single bead for each connected component resulting from severing all rotatable bonds."
Decoding is autoregressive bead-by-bead with aggregated attention back to fine-grained atoms.

## Inputs/outputs — protonation / charge / H handling
Node features "such as element type, atomic charge, and hybridization state"; edges labelled by bond type plus
auxiliary 4 Å radius edges. The paper does not state whether hydrogens are explicit => Unknown
(resolvable from the repo's featurization). No pH, no state variable.

## Data / splits
GEOM-QM9 and GEOM-DRUGS with the Ganea et al. splits (106586/13323/1000 and 243473/30433/1000). For
data-constrained training it keeps the first 5 (QM9) / 2 (DRUGS) conformers, "typically those with the largest
Boltzmann weight". NOT scaffold-separated.

## Conformer populations
Boltzmann weights are used only to select training conformers; the model predicts no energies or weights.
Objective = Kabsch-aligned coordinate MSE + KL + distance auxiliary loss.

## Relevance to p(z, x | G, h)
The coarse-bead decomposition is the closest structural analogue in this family to "generate a coarse chemical
decision, then refine geometry", but the coarse variable is purely geometric (a rigid fragment pose), never
chemical, so it is not a protonation-state mechanism. Also inherits TD/GeoMol's dependence on an RDKit conformer.
