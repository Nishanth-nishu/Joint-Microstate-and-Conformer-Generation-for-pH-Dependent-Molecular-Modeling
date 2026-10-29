#!/bin/bash
#SBATCH --job-name=jmcg_e2dbg
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# After the MMFF-copy fix: which bonds_angles failures remain, and the worst 1-2 / 1-3 pair per failing conformer.
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage
echo "job $SLURM_JOB_ID start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
rsync -a $STAGE/jmcg/ $P/jmcg/
cd $P; source venv/bin/activate; export PYTHONPATH=$P OMP_NUM_THREADS=1
python - <<'PY'
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdDistGeom
from jmcg.joint import features as ft
from jmcg.condgen import validate as va
for s in ['O=c1cccc[nH]1', 'Oc1ccccn1', 'C/C(O)=C/C(C)=O', 'O=c1[nH]cnc2[nH]cnc12']:
    r = ft.ensemble_features(s)
    print(s, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items() if k != 'failures'}, sorted(set(r['failures'])))
# worst pairs for 2-hydroxypyridine after MMFF (fixed pipeline)
s = 'Oc1ccccn1'
m = Chem.AddHs(Chem.MolFromSmiles(s)); bm = rdDistGeom.GetMoleculeBoundsMatrix(m); topo = Chem.GetDistanceMatrix(m)
mf = Chem.Mol(m)
p = AllChem.ETKDGv3(); p.randomSeed = ft.conformer_seeds(0, 1)[0]; AllChem.EmbedMolecule(mf, p)
props = AllChem.MMFFGetMoleculeProperties(mf, mmffVariant='MMFF94s'); ff = AllChem.MMFFGetMoleculeForceField(mf, props); ff.Minimize(maxIts=2000)
pos = mf.GetConformer().GetPositions(); d = np.linalg.norm(pos[:, None] - pos[None], axis=-1)
rows = []
for i in range(len(pos)):
    for j in range(i + 1, len(pos)):
        if topo[i, j] in (1, 2):
            dev = max(bm[j, i] - d[i, j], d[i, j] - bm[i, j])
            rows.append((round(dev, 3), int(topo[i, j]), m.GetAtomWithIdx(i).GetSymbol() + str(i), m.GetAtomWithIdx(j).GetSymbol() + str(j), round(d[i, j], 3), round(bm[j, i], 3), round(bm[i, j], 3)))
for x in sorted(rows, reverse=True)[:6]: print('dev,topo,a,b,d,lb,ub', x)
PY
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
