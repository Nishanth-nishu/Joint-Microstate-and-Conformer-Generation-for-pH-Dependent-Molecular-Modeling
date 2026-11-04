#!/bin/bash
#SBATCH --job-name=jmcg_t2_audit
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=01:00:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# EXP1-V2-T2 task 1: read-only audit of the saved T1 E1 outputs + the TD evaluate_confs.py code that defines the metrics.
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage; TD=/scratch/nishanth.r/tordiff; R=$TD/torsional-diffusion
echo "job $SLURM_JOB_ID start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
source $TD/venv/bin/activate
D=$STAGE/research/experiments/EXP1_V2_T2; mkdir -p $D
sha256sum $R/evaluate_confs.py
python $STAGE/experiments/exp1/t2_audit_e1.py --in_dir $P/results/e1_confs/n200 --t1_dir $P/results/exp1/T1_e1 --out $D/audit_e1.json
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
