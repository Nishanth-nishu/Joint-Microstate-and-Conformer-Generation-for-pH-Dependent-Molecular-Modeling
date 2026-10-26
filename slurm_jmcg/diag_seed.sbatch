#!/bin/bash
#SBATCH --job-name=jmcg_diag
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage
echo "job $SLURM_JOB_ID on $(hostname) start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
rsync -a $STAGE/jmcg/ $P/jmcg/; cd $P; export PYTHONPATH=$P
echo "=== STAGE 1: RDKit 2026.09 (project venv)"
venv/bin/python -m jmcg.diagnostics.embed_seed_test
echo "=== STAGE 2: RDKit 2024.03 (throwaway venv)"
[[ -x venv_rd2024/bin/python ]] || { python3.10 -m venv venv_rd2024; venv_rd2024/bin/pip install -q --no-index --find-links /home2/nishanth.r/wheels_jmcg "rdkit==2024.3.6" numpy pillow; }
venv_rd2024/bin/python -m jmcg.diagnostics.embed_seed_test
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
