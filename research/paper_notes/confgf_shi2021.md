# ConfGF (Shi et al. 2021)

## Citation
Learning Gradient Fields for Molecular Conformation Generation. ICML 2021. arXiv:2105.03902.
Code: https://github.com/DeepGraphLearning/ConfGF. Checked 2026-10-26.

## Fulltext status
abstract_only for the paper; official code read (confgf/dataset/dataset.py).

## Inputs/outputs — protonation / charge / H handling (from code)
- `mol = Chem.AddHs(Chem.MolFromSmiles(smiles))` => explicit hydrogens.
- The stored graph is `Data(atom_type=z, pos=pos, edge_index=edge_index, edge_type=edge_type, rdmol=..., smiles=...)`
  where `z` is built from `atom.GetAtomicNum()`. **No formal charge feature is stored.** Aromaticity and
  hybridization are computed in the loop but not placed in the saved `Data` object's node features.
- A charged and a neutral protomer that differ only by an H and a charge are therefore distinguishable only through
  the extra H node and the resulting atom-count/edge pattern, not through any charge channel.
- No pH, no state variable.

## Data / splits
GEOM-QM9 / GEOM-Drugs subsets (the "Shi et al. 2021" splits later reused by EC-Conf). NOT scaffold-separated.

## Relevance to p(z, x | G, h)
Baseline of historical interest only; the missing charge channel makes it the weakest candidate of the family for
charged/variant-H inputs.
