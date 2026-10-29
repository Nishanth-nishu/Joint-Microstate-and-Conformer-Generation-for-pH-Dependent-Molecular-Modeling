#!/bin/bash
#SBATCH --job-name=jmcg_m3e1
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=10:00:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# Milestone 3/4 E1 (2026-10-28): sequential-baseline state component on all experimental labels with defined states:
# 356 external macro-pKa transitions (oracle vs enumerated, two Uni-pKa templates) and 522 Tautobase pairs.
# Released Uni-pKa, 5 folds, CPU. Smoke (5 records) first; the job aborts on any error.
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage
echo "job $SLURM_JOB_ID on $(hostname) start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
rsync -a $STAGE/jmcg/ $P/jmcg/; rsync -a $STAGE/slurm_jmcg/ $P/slurm_jmcg/; rsync -a $STAGE/data_manifests/ $P/data_manifests/
cd $P; source venv/bin/activate; export PYTHONPATH=$P OMP_NUM_THREADS=8
echo "jmcg HEAD $(git rev-parse --short HEAD) + working tree; data manifest sha $(sha256sum data_manifests/canonical_v1/manifest.json | cut -c1-16)"
M=$P/models/unipka; OUT=$P/results/seq/e1; mkdir -p $OUT
python -m jmcg.seq.test_seq
python -m jmcg.seq.e1_states --data data_manifests/canonical_v1 --model_dir $M --out $OUT/smoke --limit 5 --threads 8 \
  --n1_transitions $P/results/n1_pilot/full_v2/transitions.csv > $OUT/smoke.log 2>&1 || { tail -30 $OUT/smoke.log; exit 1; }
echo "smoke ok"; t0=$(date +%s)
python -m jmcg.seq.e1_states --data data_manifests/canonical_v1 --model_dir $M --out $OUT/full --threads 8 \
  --n1_transitions $P/results/n1_pilot/full_v2/transitions.csv
echo "full seconds $(( $(date +%s) - t0 ))"
D=$STAGE/research/experiments/M3_seq/e1; mkdir -p $D
rsync -a --include '*/' --include '*.json' --include '*.csv' --include '*.log' --exclude '*' $OUT/ $D/
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
