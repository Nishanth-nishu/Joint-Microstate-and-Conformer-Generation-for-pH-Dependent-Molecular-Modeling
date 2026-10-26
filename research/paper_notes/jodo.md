# JODO
## Citation
Huang H., Sun L., Du B., Lv W., "Learning Joint 2D & 3D Diffusion Models for Complete Molecule Generation",
IEEE TNNLS 2024, DOI 10.1109/tnnls.2024.3416328; arXiv:2305.12347. Code: github.com/GRAPH-0/JODO.
PDF: research/papers/arxiv_2305.12347.pdf. Checked 2026-10-26.
## Fulltext status
abstract_only (abstract + secondary summaries).
## Method (as reported)
Generates atom types, FORMAL CHARGES, bonds and 3D coordinates together with a Diffusion Graph Transformer that
predicts clean data from noisy inputs; relational attention updates node/edge features jointly with scalar and
geometric features. Supports graph-only generation and inverse design for quantum properties; uses
self-conditioning during sampling (noted elsewhere as a possible error-compounding risk).
## Relevance to p(z,x|G,h)
Shows that a single denoiser can be trained so that *subsets* of the modalities can be generated or held fixed
(graph-only mode) - the mechanism needed for conditional p(z|G,h) or p(x|z,G) ablations within one joint model.
