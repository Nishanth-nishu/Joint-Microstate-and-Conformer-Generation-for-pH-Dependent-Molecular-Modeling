#!/bin/bash
#SBATCH --job-name=jmcg_e2feat
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=03:00:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# E2 (pre-registered, research/experiments/E2_joint/PREREG.md): conformer-ensemble features phi and 2D control psi for
# every state of the 221 measured Tautobase pairs. ETKDGv3 x30 + MMFF94s, validated. jmcg venv, CPU.
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage
echo "job $SLURM_JOB_ID on $(hostname) start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
rsync -a $STAGE/jmcg/ $P/jmcg/; rsync -a $STAGE/data_manifests/ $P/data_manifests/; rsync -a $STAGE/research/experiments/E2_joint/ $P/research/experiments/E2_joint/
cd $P; source venv/bin/activate; export PYTHONPATH=$P OMP_NUM_THREADS=1
echo "PREREG sha $(sha256sum research/experiments/E2_joint/PREREG.md | cut -c1-16)"
OUT=$P/results/e2/features; mkdir -p $OUT
python - <<'PY'
from jmcg.joint import features as ft
for s in ['C/C(O)=C/C(C)=O', 'CC(=O)CC(C)=O', 'Oc1ccccn1', 'O=c1cccc[nH]1']:
    r = ft.ensemble_features(s); r.update(ft.graph_features(s)); print(s, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in r.items()})
PY
t0=$(date +%s); python -m jmcg.joint.e2_tautomer features --data data_manifests/canonical_v1 --out $OUT --procs 8
echo "features seconds $(( $(date +%s) - t0 ))"
D=$STAGE/research/experiments/E2_joint/features; mkdir -p $D; cp $OUT/features.csv $D/
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
