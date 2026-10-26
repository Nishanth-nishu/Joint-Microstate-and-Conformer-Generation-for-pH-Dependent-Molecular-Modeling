# GeoMFormer (Chen et al. 2024)

## Citation
GeoMFormer: A General Architecture for Geometric Molecular Representation Learning. Tianlang Chen, Shengjie Luo,
Di He, Shuxin Zheng, Tie-Yan Liu, Liwei Wang. 2024. arXiv:2406.16853. Code (stated as forthcoming):
https://github.com/c-tl/GeoMFormer. Checked 2026-10-26.

## Fulltext status
full_text_read (arXiv v1 HTML, first ~100k chars).

## Problem / Method
NOT a conformer generator. A Transformer that maintains two interacting streams, invariant and equivariant, to
learn "both invariant and equivariant features" for molecular property and structure PREDICTION.
Atom input is "a learnable embedding vector indexed by the atom" type; 3D input is the mean-centred position
encoded with a Gaussian basis kernel.

## Inputs/outputs — protonation / charge / H handling
Formal charge: not mentioned. Explicit hydrogens: not mentioned. Both Unknown — resolvable from the code once
released, or from the appendix beyond the portion read.

## Data
OC20 (IS2RE, IS2RS), PCQM4Mv2, Molecule3D (3,899,647 molecules, HOMO-LUMO gap), N-body simulation, MD17.

## Relevance to p(z, x | G, h)
Belongs in this family as the "shared encoder for invariant + equivariant targets" architecture: if a joint model
needed one backbone that both classifies a state (invariant head) and places atoms (equivariant head), this is the
closest published template. It performs no generation and no state reasoning itself.
