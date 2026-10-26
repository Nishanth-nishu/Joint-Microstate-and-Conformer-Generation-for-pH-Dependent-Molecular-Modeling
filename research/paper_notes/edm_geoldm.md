# EDM (and the EDM-lineage: GeoLDM, DiffSBDD)
## Citation
Hoogeboom E., Satorras V.G., Vignac C., Welling M., "Equivariant Diffusion for Molecule Generation in 3D",
ICML 2022, arXiv:2203.17003. Code: github.com/ehoogeboom/e3_diffusion_for_molecules. Checked 2026-10-26.
## Fulltext status
abstract_only.
## Method (as reported)
E(3)-equivariant denoising diffusion over continuous coordinates jointly with node features. On QM9 the node
features are atom type (one-hot) AND an INTEGER-VALUED ATOMIC CHARGE, i.e. the categorical/ordinal chemical label
is relaxed into continuous space and diffused with Gaussian noise alongside coordinates (the "scaling" trick of
dividing one-hots by a scalar before diffusion, later reused and criticised by MolDiff).
## Relevance to p(z,x|G,h)
Earliest and simplest coupling: single Gaussian diffusion over [coords | relaxed one-hots]. It is the baseline
parameterisation against which MiDi/EQGAT-diff/FlowMol argue for true categorical noise - the same choice the
project faces for z.
## Unknowns
How the atom count is sampled (prior over N) was not confirmed from the sources read; relevant to variable-size
protonation (adding/removing H changes atom count). Resolvable from the EDM paper Sec. on the size prior.
