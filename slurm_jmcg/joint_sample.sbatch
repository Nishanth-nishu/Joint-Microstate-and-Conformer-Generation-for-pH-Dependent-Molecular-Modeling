#!/bin/bash
#SBATCH --job-name=jmcg_jsamp
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# M5 demo: joint sampler on the 6 smoke parents at pH 7.4, arms S0 / S1 / J1 (E2 final coefficients), seed 0.
# Check: S0 state populations equal the sequential state stage's p_mean at pH 7.4 (same G_U, same formula).
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage
echo "job $SLURM_JOB_ID start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
rsync -a --checksum $STAGE/jmcg/ $P/jmcg/; md5sum $P/jmcg/joint/*.py | sed "s|$P/||"
mkdir -p $P/results/joint/smoke_states; rsync -a $STAGE/research/experiments/M3_seq/smoke/states/ $P/results/joint/smoke_states/
cd $P; source venv/bin/activate; export PYTHONPATH=$P OMP_NUM_THREADS=1
for arm in S0 S1 J1; do
  python -m jmcg.joint.sample --state_dir results/joint/smoke_states --coef results/e2/fit/e2_result.json --arm $arm \
    --ph 7.4 --out results/joint/sample_$arm --n_pairs 20 --seed 0 > results/joint/sample_$arm.log 2>&1 || { tail -20 results/joint/sample_$arm.log; exit 1; }
  echo "== $arm"; cat results/joint/sample_$arm/run_info.json
done
python - <<'PY'
import pandas as pd
P = pd.read_csv('results/joint/smoke_states/populations.csv'); P = P[(P.pH - 7.4).abs() < 1e-9]
for arm in ('S0', 'S1', 'J1'):
    J = pd.read_csv(f'results/joint/sample_{arm}/joint_state_populations.csv')
    m = J.merge(P[['state_id', 'p_mean']], on='state_id', how='left').fillna({'p_mean': 0.0})
    print(arm, "max |p_joint - p_seq(state stage)| =", round(float((m.p_joint - m.p_mean).abs().max()), 6), "(state stage averages fold populations; sampler softmaxes the fold-mean G)")
    top = J.sort_values('p_joint', ascending=False).groupby('parent').head(2)
    print(top[["parent", "smiles", "charge", "p_joint", "phi_ok", "coupled_score_used"]].to_string(index=False))
    S = pd.read_csv(f'results/joint/sample_{arm}/sampled_pairs.csv'); print(arm, 'pairs', len(S), 'with conformer', int(S.has_conformer.sum()))
PY
python - <<'PY'
# S0 identity check at matched averaging: softmax of fold-mean G over each parent's enumerated states (pH 7.4)
import numpy as np, pandas as pd
from jmcg.seq import score as sc
ens = pd.read_csv('results/joint/smoke_states/ensembles.csv')
G = pd.read_csv('results/joint/smoke_states/energies.csv').groupby('smiles').G.mean()
J = pd.read_csv('results/joint/sample_S0/joint_state_populations.csv').set_index('state_id')
d = 0.0
for pid, g in ens.groupby('parent'):
    p = sc.populations(G[g.smiles].to_numpy(), g.charge.to_numpy(), 7.4)
    d = max(d, float(np.abs(p - J.p_joint[g.state_id].to_numpy()).max()))
print('S0 vs softmax(fold-mean G) max abs diff', d)
PY
D=$STAGE/research/experiments/M5_joint_sampler; mkdir -p $D
for arm in S0 S1 J1; do mkdir -p $D/$arm; cp results/joint/sample_$arm/* $D/$arm/; done
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
