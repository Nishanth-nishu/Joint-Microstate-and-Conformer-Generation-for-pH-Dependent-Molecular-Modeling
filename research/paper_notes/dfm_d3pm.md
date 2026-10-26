# Discrete generative machinery: D3PM and Discrete Flow Matching
## Citation
Austin J. et al., "Structured Denoising Diffusion Models in Discrete State-Spaces" (D3PM), NeurIPS 2021,
arXiv:2107.03006. Gat I. et al., "Discrete Flow Matching", NeurIPS 2024 (spotlight), arXiv:2407.15595. Checked 2026-10-26.
## Fulltext status
abstract_only.
## Method (as reported)
D3PM: discrete diffusion with structured transition matrices beyond uniform - Gaussian-kernel-like matrices,
nearest-neighbour-in-embedding matrices, and ABSORBING-STATE matrices. Gat et al.: discrete flow matching via CTMCs
with general probability paths, non-autoregressive parallel sampling.
## Relevance to p(z,x|G,h)
D3PM's structured transition matrices are the natural place to encode chemical priors on z: transitions between
protonation microstates should be restricted to single proton add/remove moves between adjacent net charges, not
uniform jumps across the whole microstate set. No paper found does this for protonation.
