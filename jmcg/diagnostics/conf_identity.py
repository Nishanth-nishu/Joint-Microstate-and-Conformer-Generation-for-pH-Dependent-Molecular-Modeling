"""Diagnose why the 10 3D conformers give identical Uni-pKa free energies in the N1 pilot."""
import sys
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolAlign
from unipkainfer.pka_predictor import free_energy as fe
import jmcg.experiments.n1_pilot  # noqa: F401  applies the faithful Uni-pKa conformer protocol to fe

run = sys.argv[1]
c = pd.read_csv(f'{run}/conformer_free_energies_fold0.csv')
p = c.pivot_table(index='smiles_c', columns='conformer_index', values='conformer_free_energy')
spread3d = p.iloc[:, :10].max(axis=1) - p.iloc[:, :10].min(axis=1)
print('microstates', len(p), 'with any 3D spread > 1e-6:', int((spread3d > 1e-6).sum()))
print('3D==2D copy (all 10):', int(p.iloc[:, :10].eq(p.iloc[:, 10], axis=0).all(axis=1).sum()))
# most flexible microstates: check coordinate diversity of what unipkainfer generates
smis = sorted(p.index, key=lambda s: -Chem.rdMolDescriptors.CalcNumRotatableBonds(Chem.MolFromSmiles(s)))[:5]
for s in smis:
    coords = fe._smiles_to_3d_coords(s, 10, 'mmff', num_threads=4)
    mol = Chem.AddHs(Chem.MolFromSmiles(s))
    heavy = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    d = [np.sqrt(((coords[0][heavy] - x[heavy]) ** 2).sum(1).mean()) for x in coords[1:]]
    print(s, 'rot', Chem.rdMolDescriptors.CalcNumRotatableBonds(mol), 'raw coord RMS diff vs conf0:',
          np.round(d, 3), 'G spread', float(spread3d[s]))
