# ET-Flow (Hassan, Shenoy et al. 2024)

## Citation
ET-Flow: Equivariant Flow-Matching for Molecular Conformer Generation. NeurIPS 2024. arXiv:2410.22388.
Code: https://github.com/shenoynikhil/ETFlow. Checked 2026-10-26.

## Fulltext status
abstract_only for the paper (arXiv abstract page) + official code read (etflow/commons/featurization.py).
Methodological statements below are restricted to the abstract and to the code.

## Method (abstract + code)
Equivariant flow matching with "equivariance and harmonic prior", operating "directly on all-atom coordinates
with minimal assumptions" (Abstract). Backbone is a TorchMD-NET / equivariant transformer (etflow/networks/torchmd_net).

## Inputs/outputs — protonation / charge / H handling (from code)
- `mol = Chem.AddHs(mol)` in featurization.py => explicit hydrogens are nodes.
- `[atom.GetFormalCharge() for atom in mol.GetAtoms()]` is computed per atom => formal charge is an input feature.
- Chiral centers extracted (`get_chiral_tensors`) and carried as chiral_index / chiral_nbr_index / chiral_tag.
- No pH, no categorical state, one graph per molecule.

## Data
GEOM-QM9 / GEOM-DRUGS with the standard Ganea et al. splits (as reported by DiTMC and SO(3)-AvgFlow, which
reproduce ET-Flow numbers on those splits). Unknown whether any scaffold split was attempted — would be resolved
by reading Sec. 5 / Appendix of the paper.

## Relevance to p(z, x | G, h)
Viable conditional geometry generator: accepts charged graphs and explicit H, generates all-atom Cartesians, and is
cheap at inference (few ODE steps). Same gap as DiTMC: no z, no pH, no population estimate.
