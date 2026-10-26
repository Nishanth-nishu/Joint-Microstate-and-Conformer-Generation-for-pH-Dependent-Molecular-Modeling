#!/bin/bash
#SBATCH --job-name=jmcg_diag
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage
echo "job $SLURM_JOB_ID on $(hostname) start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
echo "=== STAGE 1: sync"
rsync -a $STAGE/jmcg/ $P/jmcg/
cd $P; source venv/bin/activate; export PYTHONPATH=$P
echo "=== STAGE 2: diagnose"
python -m jmcg.diagnostics.embed_test
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
