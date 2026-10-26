# EQGAT-diff
## Citation
Le T., Cremer J., Noe F., Clevert D.-A., Schutt K., "Navigating the Design Space of Equivariant Diffusion-Based
Generative Models for De Novo 3D Molecule Generation", ICLR 2024, arXiv:2309.17296.
Code: github.com/tuanle618/eqgat-diff. PDF: research/papers/arxiv_2309.17296.pdf. Checked 2026-10-26.
## Fulltext status
abstract_only (abstract, ICLR slides summary, and a 2026 survey description).
## Method (as reported)
E(3)-equivariant graph attention denoiser. Continuous diffusion on atom positions; CATEGORICAL diffusion on
chemical elements, FORMAL CHARGES and bond types; x0-parameterisation (the denoiser predicts clean coordinates,
elements, formal charges and bond types directly). Truncated-SNR time-dependent loss weighting; chemically
motivated input features (hybridisation, aromaticity, ring membership). Large-scale pretraining on PubChem3D with
implicit hydrogens then fine-tuning on explicit-hydrogen datasets.
## Relevance to p(z,x|G,h)
A second independent confirmation that formal charge is routinely generated as a categorical variable jointly with
geometry, and the clearest demonstration that implicit-H pretraining transfers to explicit-H targets - relevant
because protonation state is fundamentally an explicit-hydrogen question. No pH/protonation semantics.
