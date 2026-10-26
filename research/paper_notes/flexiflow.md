# FlexiFlow
## Citation
Tedoldi R., Engkvist O., Bryant P., Azizpour H., Janet J.P., Tibo A., "FlexiFlow: decomposable flow matching for
generation of flexible molecular ensemble", arXiv:2511.17249, 2025-11-21.
PDF: research/papers/arxiv_2511.17249.pdf. Checked 2026-10-26.
## Fulltext status
full_text_read (arXiv HTML v1 via WebFetch; Appendix B partially truncated at 100k chars).
## Problem
Existing 3D molecule generators emit ONE conformation per molecule; drug discovery needs the molecule together
with its low-energy conformational ensemble.
## Method (mechanism level)
Flow matching extended to sample a molecule AND multiple conformations jointly, preserving equivariance and
permutation invariance. Node features h concatenate two one-hot vectors, "one encoding the atom type and the other
encoding the charge type" (Sec. 4.2.1); three shared MLPs predict logits for atom and charge types; Algorithm 1
samples an initial charge state from a uniform categorical and applies the Campbell et al. (2024) discrete update
at inference, while coordinates use a flow-matching loss. Invariant features (atom + charge type, time) condition
both coordinate sets, linked through scalar products (Eq. 13). Trained with explicit hydrogens (Sec. 5.1).
Transfers to pocket-conditioned ligand generation.
## Experiments & baselines
QM9 and GEOM-Drugs; ensemble coverage reportedly matching physics-based methods at much lower inference cost.
## Relevance to p(z,x|G,h)
The closest existing architecture on the x-side of the claim: it is the only found model that generates a
DISCRETE CHARGE STATE jointly with MULTIPLE CONFORMATIONS of the same molecule in one flow. What is missing is
exactly the project's claim: the charge categorical is a per-atom formal-charge label used for valence bookkeeping,
NOT a pH-dependent protonation microstate, and the paper never mentions protonation, tautomers or pH.
## Exact supporting quotes
- Sec. 4.2.1: "one encoding the atom type and the other encoding the charge type".
- Sec. 5.1: trained "with the hydrogens".
