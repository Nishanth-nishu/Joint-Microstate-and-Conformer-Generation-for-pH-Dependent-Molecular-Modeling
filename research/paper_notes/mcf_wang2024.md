# Molecular Conformer Fields / MCF (Wang et al. 2024)

## Citation
Swallowing the Bitter Pill: Simplified Scalable Conformer Generation (Molecular Conformer Fields).
ICML 2024. arXiv:2311.17932. Checked 2026-10-26.

## Fulltext status
full_text_read (arXiv v2 HTML).

## Method
Diffusion over fields: a PerceiverIO score network maps graph-positional coordinates to 3D coordinates, so no
torsion-angle or equivariance assumptions. "The atomic features are concatenated with graph Laplacian eigenvectors
in both context and query inputs." k = 28 (QM9) / 32 (DRUGS) leading eigenvectors of the normalized Laplacian.

## Inputs/outputs — protonation / charge / H handling
Table 6 node features: atom type one-hot over the 35 dataset elements; degree 0-6; **formal charge -1..+1**;
implicit valence 0-6; hybridization; aromaticity; number of rings 0-3. Formal charge is present but only over
the narrow {-1,0,+1} range. Whether hydrogens are explicit is not stated in the text read => Unknown
(resolvable from the released preprocessing code, which was not located).

## Data / splits
GEOM-QM9 / GEOM-DRUGS / GEOM-XL. "We split them randomly based on molecules into train/validation/test
(80%/10%/10%)" giving 106586/13323/1000 and 243473/30433/1000. NOT scaffold-separated.

## Relevance to p(z, x | G, h)
Scales well and is graph-agnostic, so a changed protomer is simply a different graph with different Laplacian
eigenvectors — no architectural obstacle. But the narrow charge encoding and the unconfirmed H treatment make it a
weaker first choice than DiTMC. No z, no pH, no populations.
