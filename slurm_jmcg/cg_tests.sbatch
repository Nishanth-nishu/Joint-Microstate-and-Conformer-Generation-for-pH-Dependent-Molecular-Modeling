#!/bin/bash
#SBATCH --job-name=jmcg_cgtest
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# Conditional generation: unit tests only (states + backend), no generation run.
P=/scratch/nishanth.r/jmcg; TD=/scratch/nishanth.r/tordiff
echo "job $SLURM_JOB_ID start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
rsync -a /home2/nishanth.r/jmcg_stage/jmcg/condgen/ $P/jmcg/condgen/
source $TD/venv/bin/activate; cd $P
export PYTHONPATH=$P TD_REPO=$TD/torsional-diffusion OMP_NUM_THREADS=4 CUDA_VISIBLE_DEVICES=""
python -m jmcg.condgen.test_states 2>/dev/null | tail -1
python -m jmcg.condgen.test_backend --model_dir $TD/workdir/pretrained/qm9_1order --states_csv jmcg/condgen/smoke_states.csv 2>&1 | grep -E '^(PASS|FAIL|test_backend)'
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
