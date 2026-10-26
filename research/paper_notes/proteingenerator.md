# ProteinGenerator (RFdiffusion sequence-space)
## Citation
Lisanza S.L. et al., "Multistate and functional protein design using RoseTTAFold sequence space diffusion",
bioRxiv 2023.05.08.539766 (later Nat. Biotechnol.). Code: github.com/RosettaCommons/protein_generator. Checked 2026-10-26.
## Fulltext status
abstract_only.
## Method (as reported)
Diffusion in CONTINUOUS sequence space (noised one-hot sequence) with RoseTTAFold2 predicting the structure at
every denoising step, so the structure distribution is never modelled explicitly - an alternative coupling to
joint-denoiser designs: the discrete variable is diffused, the continuous one is a deterministic readout.
Authors argue simultaneous generation "enables coordinated guidance with both sequence and structural features".
## Relevance to p(z,x|G,h)
An architecture option where z is diffused and x is produced by a frozen predictor conditioned on the current z -
cheap to implement over a reproduced Torsional Diffusion (z diffused, TD as the x readout), and the stated
rationale (coordinated vs alternating) is the exact claim this project wants to test.
