#!/bin/bash
#SBATCH --job-name=jmcg_t1_leak
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=01:00:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# EXP1-V2-T1: E1 test molecules vs the paired training/validation caches (read only).
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage; TD=/scratch/nishanth.r/tordiff
set -euo pipefail
source $TD/venv/bin/activate
python $STAGE/experiments/exp1/t1_leak_check.py $P/results/e1_confs/n200/subset.csv $TD/data/QM9/cache_paired.train $TD/data/QM9/cache_paired.val
