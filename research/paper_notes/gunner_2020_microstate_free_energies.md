# Standard state free energies, not pKas (Gunner et al. 2020)

## Citation
Gunner M. R., Murakami T., Rustenburg A. S., Işık M., Chodera J. D.
"Standard state free energies, not pKas, are ideal for describing small molecule protonation and
tautomeric states." J. Comput. Aided Mol. Des. 34:561 (2020). DOI 10.1007/s10822-020-00280-7.
Checked 2026-10-26.

## Fulltext status
abstract_only (lab publication page summary).

## Content (as stated)
Argues that "the standard state free energies of each microscopic protonation/tautomeric state at a
single pH is sufficient" to reconstruct all pH-dependent microscopic and macroscopic populations,
and that blind challenges should report microstate free energies at a reference pH as the minimum
reporting unit.

## Relevance to p(z,x|G,h)
Gives the correct target parameterisation for the discrete part z of the model: learn
ΔG°(microstate) at a reference pH, then obtain p(z | pH) analytically. It also implies the
population benchmark must be built from free energies/ratios, never from counts of enumerated states.
