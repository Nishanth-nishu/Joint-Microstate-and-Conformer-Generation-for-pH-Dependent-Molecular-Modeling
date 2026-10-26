# MolDiff and MUDiff (atom-bond consistency; unified 2D+3D)
## Citation
Peng X. et al., "MolDiff: Addressing the Atom-Bond Inconsistency Problem in 3D Molecule Diffusion Generation",
ICML 2023, arXiv:2305.07508, code github.com/pengxingang/MolDiff.
Hua C. et al., "MUDiff: Unified Diffusion for Complete Molecule Generation", LoG 2023, arXiv:2304.14621. Checked 2026-10-26.
## Fulltext status
abstract_only.
## Method (as reported)
MolDiff: diffuses atoms AND bonds simultaneously with explicit modelling of their dependency, because generating
atoms first and inferring bonds post hoc can leave no valid bond assignment for the chosen geometry. It uses a
two-stage/unequal schedule in which bond types are perturbed with higher noise than atom types and positions
(bonds recovered from atom information). Its appendix notes that continuous-diffusion variants treat one-hot
encodings as continuous values with Gaussian noise (the EDM scaling trick).
MUDiff: discrete diffusion over the 2D graph + continuous diffusion over 3D coordinates, denoised by an
equivariant graph transformer.
## Relevance to p(z,x|G,h)
MolDiff is the sharpest statement of WHY sequential discrete-then-continuous pipelines fail: the discrete choice
made without geometry may be geometrically unrealisable. That argument transfers verbatim to protonation:
picking z by rules/pKa model without geometry can yield a state whose realisable conformers are wrong. Its
unequal-noise schedule is a concrete design precedent.
