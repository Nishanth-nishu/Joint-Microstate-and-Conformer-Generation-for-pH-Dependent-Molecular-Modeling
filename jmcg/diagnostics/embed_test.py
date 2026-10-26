"""Isolate why EmbedMultipleConfs/MMFF gives identical conformers (RDKit version, threads, MMFF)."""
import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

print('rdkit', rdBase.rdkitVersion)
smi = 'COc1ccc(CNC(=O)[C@@H](NC(=O)[C@H](Cc2ccc(OC)c(OC)c2)NC(=O)Cc2cccc(Oc3ccccc3)c2)C(C)C)c(O)c1'


def spread(mol, ids):
    x0 = mol.GetConformer(int(ids[0])).GetPositions()
    return [round(float(np.sqrt(((mol.GetConformer(int(i)).GetPositions() - x0) ** 2).sum(1).mean())), 3) for i in ids[1:]]


for threads in (1, 4):
    for mmff in (False, True):
        mol = Chem.AddHs(Chem.MolFromSmiles(smi))
        ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=10, randomSeed=0, numThreads=threads))
        if mmff:
            AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=threads)
        print(f'threads={threads} mmff={mmff} n={len(ids)} ids={ids} rms-vs-conf0={spread(mol, ids)}', flush=True)
p = AllChem.ETKDGv3(); p.randomSeed = 0
mol = Chem.AddHs(Chem.MolFromSmiles(smi)); ids = list(AllChem.EmbedMultipleConfs(mol, 10, p))
print('ETKDGv3 params object:', spread(mol, ids))
