#!/bin/bash
#SBATCH --job-name=jmcg_cginspect
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# Read-only inventory for the wrapper-fidelity check (TD scripts, QM9 test split, checkpoints).
TD=/scratch/nishanth.r/tordiff; R=$TD/torsional-diffusion
echo "job $SLURM_JOB_ID start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
ls $TD $R $TD/workdir $TD/workdir/pretrained; ls -la $TD/workdir/pretrained/*
ls -la $R/data 2>&1 | head -30
find $TD -maxdepth 5 \( -name 'test_mols*' -o -name '*.csv' -o -name 'split*.npy' -o -name '*.pkl' \) 2>/dev/null | grep -v venv | head -40
echo "=== generate_confs.py args"; grep -n "add_argument" $R/generate_confs.py
echo "=== evaluate_confs.py"; cat $R/evaluate_confs.py
echo "=== generate_confs.py"; cat $R/generate_confs.py
echo "=== model_parameters"; for d in $TD/workdir/pretrained/*; do echo "-- $d"; cat $d/model_parameters.yml 2>/dev/null | head -40; done
nproc; free -g | head -2; nvidia-smi --query-gpu=name,memory.used,memory.total --format=csv 2>&1 | head
