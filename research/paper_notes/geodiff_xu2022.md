# GeoDiff (Xu et al. 2022)

## Citation
GeoDiff: a Geometric Diffusion Model for Molecular Conformation Generation. ICLR 2022. arXiv:2203.02923.
Code: https://github.com/MinkaiXu/GeoDiff. Checked 2026-10-26.

## Fulltext status
abstract_only. The repo's utils/chem.py was fetched but contains only bond-type constants
(`BOND_TYPES = {t: i for i, t in enumerate(BT.names.values())}`); the dataset builder was not retrieved
(raw.githubusercontent timeouts). Featurization claims below are therefore NOT asserted.

## Method (abstract level)
Equivariant Markov denoising diffusion directly on atomic Cartesian coordinates, with roto-translation invariant
likelihood via an equivariant Markov kernel.

## Inputs/outputs — protonation / charge / H handling
Unknown. GeoDiff reuses the ConfGF data pipeline in the published benchmarks (TD states "GeoDiff originally used a
small subset of the DRUGS dataset, so we retrained it using the splits from Ganea et al. (2021)"), which would imply
atom-type-only node features and no formal charge — but this was NOT verified against GeoDiff's own code in this run.
What would resolve it: read `utils/datasets.py` / `models/epsnet/*.py` node-feature construction in MinkaiXu/GeoDiff.

## Relevance to p(z, x | G, h)
Cartesian diffusion means a protomer change is expressible in the output space, but without a confirmed charge
feature the conditioning side is doubtful. Superseded on quality by ET-Flow / DiTMC / MCF.
