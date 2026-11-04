#!/bin/bash
#SBATCH --job-name=jmcg_t2_geomval
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=03:00:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# EXP1-V2-T2 task 3.4-3.5 / 5: identity, stereo, clashes and local-geometry error decomposition of every E1 arm
# (TD, L(proj)fix, L(noproj)fix, constrained xTB, MMFF, ETKDGv3), seeds 0 and 1. Read only on the saved outputs.
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage; TD=/scratch/nishanth.r/tordiff; R=$TD/torsional-diffusion
echo "job $SLURM_JOB_ID start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
source $TD/venv/bin/activate
D=$STAGE/research/experiments/EXP1_V2_T2; mkdir -p $D
for s in 0 1; do
  python $STAGE/experiments/exp1/t2_geom_validate.py --in_dir $P/results/e1_confs/n200 --t1_dir $P/results/exp1/T1_e1 \
    --true_mols $R/data/QM9/test_mols.pkl --seed $s --workers 16 --out $D/geom_validate_seed$s.json
done
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
