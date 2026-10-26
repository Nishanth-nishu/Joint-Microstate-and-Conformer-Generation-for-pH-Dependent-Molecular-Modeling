"""Which random seeds give diverse conformers in this RDKit build?"""
import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

print('rdkit', rdBase.rdkitVersion, flush=True)
smi = 'COc1ccc(CNC(=O)[C@@H](NC(=O)[C@H](Cc2ccc(OC)c(OC)c2)NC(=O)Cc2cccc(Oc3ccccc3)c2)C(C)C)c(O)c1'
for seed in (0, 1, 42, -1):
    mol = Chem.AddHs(Chem.MolFromSmiles(smi))
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=10, randomSeed=seed))
    x0 = mol.GetConformer(ids[0]).GetPositions()
    d = [round(float(np.sqrt(((mol.GetConformer(i).GetPositions() - x0) ** 2).sum(1).mean())), 2) for i in ids[1:]]
    print(f'seed={seed} n={len(ids)} rms-vs-conf0={d}', flush=True)
