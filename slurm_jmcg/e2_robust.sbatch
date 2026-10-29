#!/bin/bash
#SBATCH --job-name=jmcg_e2rob
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# E2 exploratory robustness (post hoc, NOT pre-registered) + check that the refactored fit reproduces job 45572 exactly.
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage
echo "job $SLURM_JOB_ID start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
rsync -a $STAGE/jmcg/ $P/jmcg/
cd $P; source venv/bin/activate; export PYTHONPATH=$P OMP_NUM_THREADS=1
F=$P/results/e2/features_v2/features.csv; EN=results/seq/e1/full/energies.csv
python -m jmcg.joint.e2_tautomer fit --data data_manifests/canonical_v1 --features $F --energies $EN --out $P/results/e2/fit_refactor_check > /dev/null
cmp $P/results/e2/fit/e2_result.json $P/results/e2/fit_refactor_check/e2_result.json && echo "refactor reproduces 45572: e2_result.json identical"
cmp $P/results/e2/fit/oof_predictions.csv $P/results/e2/fit_refactor_check/oof_predictions.csv && echo "oof_predictions.csv identical"
python -m jmcg.joint.e2_tautomer robust --data data_manifests/canonical_v1 --features $F --energies $EN --out $P/results/e2/robust
mkdir -p $STAGE/research/experiments/E2_joint/robust; cp $P/results/e2/robust/* $STAGE/research/experiments/E2_joint/robust/
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
