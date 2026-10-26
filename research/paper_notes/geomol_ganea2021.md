# GeoMol (Ganea et al. 2021)

## Citation
GeoMol: Torsional Geometric Generation of Molecular 3D Conformer Ensembles. NeurIPS 2021. arXiv:2106.07802.
Code: https://github.com/PattanaikL/GeoMol. Checked 2026-10-26.

## Fulltext status
abstract_only for the paper; official code read (model/featurization.py). Method statements below come from the
code and from descriptions in papers that re-implement it (TD, DiTMC), not from GeoMol's own methods section.

## Inputs/outputs — protonation / charge / H handling (from code)
- `mol = Chem.AddHs(mol)` (featurization.py) => explicit hydrogens.
- `atom_features.extend(one_k_encoding(atom.GetFormalCharge(), [-1, 0, 1]))` => formal charge feature, narrow
  {-1,0,+1} bucket (TD's featurization is directly derived from this file).
- `qm9_types = {'H','C','N','O','F'}`; `drugs_types` = 35 elements.
- Canonical SMILES comparison uses `Chem.RemoveHs(mol)` for graph identity checks.
- No pH, no state variable.

## Data / splits
GEOM-QM9 / GEOM-DRUGS. GeoMol defines the random 243473/30433/1000 and 106586/13323/1000 splits that essentially
all later conformer papers (TD, MCF, ET-Flow, DiTMC, SO(3)-AvgFlow, CoarsenConf) reuse. NOT scaffold-separated.

## Relevance to p(z, x | G, h)
Historically important because it is the source of the field's standard split and of the featurization that TD
inherited — including the {-1,0,+1} charge bucket that would bottleneck multiply-charged protomers. One-shot
(non-iterative) local-structure + torsion prediction makes it fast but low-diversity.
