# DiTMC (Frank et al. 2025)

## Citation
Sampling 3D Molecular Conformers with Diffusion Transformers. J. Thorben Frank, Winfried Ripken, Gregor Lied,
Klaus-Robert Müller, Oliver Unke, Stefan Chmiela. NeurIPS 2025. arXiv:2506.15378.
Code: https://github.com/ML4MolSim/dit_mc. Checked 2026-10-26.

## Fulltext status
full_text_read (arXiv v2 HTML, incl. Appendix Table A8).

## Problem / Method
Adapt diffusion transformers to conformer generation. "a modular architecture that separates the processing of
3D coordinates from conditioning on atomic connectivity" (Abstract). Two graph conditioning routes: atom-wise
tokens from a MeshGraphNet-style GNN (Eq. 6) and pair-wise tokens from graph geodesic distances s(i,j),
Graphormer-style (Eq. 7): "This formulation allows conditioning on all atom pairs, including those not directly
connected by a bond." (Sec. 4.1). Three attention variants: aPE, rPE, SO(3)-equivariant PE(3). Harmonic prior.

## Data
GEOM-QM9 (133,258) and GEOM-DRUGS/AICures (304,466); GEOM-XL. "The data splits are taken from Ref. [64]"
(= GeoMol / Ganea et al.) — random molecule split, NOT scaffold-separated.

## Inputs/outputs — protonation / charge / H handling
Table A8 ("Atomic and bond features included in DiT-MC. All features are one-hot encoded"):
chirality; number of hydrogens 0-4; number of radical electrons; atom type (QM9: H,C,N,O,F / DRUGS: 35 elements
starting with H); aromaticity; degree; hybridization; implicit valence; **formal charge -5,-4,...,5, other**;
ring membership 3-8; number of rings. Bond type single/double/triple/aromatic. Features computed with RDKit.
=> explicit hydrogens are atoms, and formal charge is represented over a much wider range than TD's {-1,0,1}.
No pH input, no categorical state variable, one graph per sample.

## Experiments & baselines
AMR/COV vs TD, MCF, ET-Flow; ensemble properties (E, μ, Δε, E_min) via GFN2-xTB relaxation and Boltzmann-weighted
comparison — post-hoc evaluation only, not a model output. Post-hoc chirality correction by reflection (App. G.4).

## Limitations
"DiTMC scales quadratic in the number of atoms" (Conclusion). Conditioning is "specifically tailored to molecular
graphs as inputs", so extending beyond a fixed input graph would need a redesign (stated by the authors).

## Relevance to p(z, x | G, h)
Strongest drop-in conditional geometry generator x ~ p(x | G) for charged/variant-H inputs: it already encodes
formal charge over -5..+5 and explicit H counts, and generates full Cartesian coordinates (so ring and X-H geometry
changes caused by a protonation change are expressible — unlike torsion-only TD). What it lacks for this project is
any z: no state variable, no pH, no population output.
