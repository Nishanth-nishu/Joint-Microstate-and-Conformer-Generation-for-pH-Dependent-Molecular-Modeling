#!/bin/bash
#SBATCH --job-name=jmcg_e2dbg
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# Hypothesis test: MMFFGetMoleculeProperties mutates aromaticity flags of the Mol it is given, so features.py validates
# against MMFF-aromaticity bond types / bounds and rejects valid conformers (2-pyridone rt_bonds, 2-hydroxypyridine 1-2).
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage
echo "job $SLURM_JOB_ID start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
cd $P; source venv/bin/activate; export PYTHONPATH=$P OMP_NUM_THREADS=1
python - <<'PY'
from rdkit import Chem
from rdkit.Chem import AllChem
for s in ['O=c1cccc[nH]1', 'Oc1ccccn1', 'c1ccccc1O', 'O=c1[nH]cnc2[nH]cnc12']:
    m = Chem.AddHs(Chem.MolFromSmiles(s))
    before = [(b.GetIdx(), str(b.GetBondType()), b.GetIsAromatic()) for b in m.GetBonds()]
    AllChem.MMFFGetMoleculeProperties(m, mmffVariant='MMFF94s')
    after = [(b.GetIdx(), str(b.GetBondType()), b.GetIsAromatic()) for b in m.GetBonds()]
    ch = [(x, y) for x, y in zip(before, after) if x != y]
    print(s, 'bonds changed by MMFFGetMoleculeProperties:', len(ch), ch[:4])
PY
python - <<'PY'
import pandas as pd
F = pd.read_csv('results/e2/features/features.csv')
f = F.failures.fillna('')
print('states', len(F), 'with phi', F.phi1.notna().sum())
for k in ['rt_bonds', 'bonds_angles', 'rt_isomeric', 'clash', 'stereo', 'connectivity', 'rt_atoms', 'MMFF94s parameters missing', 'embedding failed']:
    print(k, f.str.contains(k).sum())
print('aromatic SMILES among no-phi:', F[F.phi1.isna()].smiles.str.contains('[a-z]', regex=True).sum(), 'of', F.phi1.isna().sum())
PY
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
