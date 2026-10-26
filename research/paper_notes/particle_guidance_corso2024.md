# Particle Guidance (Corso et al. 2024)

## Citation
Particle Guidance: non-I.I.D. Diverse Sampling with Diffusion Models. Gabriele Corso, Yilun Xu,
Valentin de Bortoli, Regina Barzilay, Tommi Jaakkola. ICLR 2024. arXiv:2310.13102.
Code: https://github.com/gcorso/particle-guidance. Checked 2026-10-26.

## Fulltext status
abstract_only (arXiv abstract + ICLR listing). No mechanism claims beyond the abstract are made here.

## Method (abstract level)
Instead of sampling a batch i.i.d., a joint time-evolving potential over the whole set of particles is added to the
diffusion sampler, so the set is diverse by construction.

## Results relevant here
Applied ON TOP OF torsional diffusion for conformer generation: the abstract reports reducing the state-of-the-art
median error by 13% on average, and the ICLR listing reports median precision/coverage error reductions of
19% and 8% respectively.

## Inputs/outputs — protonation / charge / H handling
Not applicable at the method level: particle guidance is a sampler modification and inherits whatever featurization
the base model uses (here, torsional diffusion's — explicit H, formal charge in {-1,0,+1}).

## Relevance to p(z, x | G, h)
Two uses. (1) A free accuracy/diversity gain on whatever conditional geometry generator is chosen for the
sequential baseline. (2) A cautionary point for population estimation: particle guidance deliberately makes the
sample set non-i.i.d., so sample frequencies from a particle-guided sampler must NOT be read as Boltzmann
populations without correction.
