#!/bin/bash
#SBATCH --job-name=jmcg_v2smoke_chk
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# Stage 4 only of exp1_v2_smoke_cpu.sbatch (job 46744 trained, then v2_smoke.py failed reading model_parameters.yml:
# torch.__version__ was dumped as a python object). Re-runs the smoke checks on the existing 46744 run; no retraining.
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage; TD=/scratch/nishanth.r/tordiff
echo "job $SLURM_JOB_ID on $(hostname) start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
rsync -a --checksum $STAGE/flexitors/ $P/flexitors/; rsync -a --checksum $STAGE/experiments/exp1/ $P/experiments/exp1/
source $TD/venv/bin/activate; cd $P
export FLEXITORS_TD_ROOT=$TD/torsional-diffusion PYTHONPATH=$P OMP_NUM_THREADS=10 MKL_NUM_THREADS=10 OPENBLAS_NUM_THREADS=10
CV=$TD/data/QM9/cache_paired.val; W=$P/workdir/exp1; O=$P/results/exp1/v2_smoke
D=$STAGE/research/experiments/EXP1_V2_smoke; mkdir -p $D
cp $O/train.log $O/cpu_resume_*.log $D/ 2>/dev/null || true
cp $W/v2_smoke_s0/{epochs.jsonl,train_log.jsonl,model_parameters.yml} $D/
cp $W/v2_cpu_resume/epochs.jsonl $D/cpu_resume_epochs.jsonl; cp $W/v2_cpu_resume/train_log.jsonl $D/cpu_resume_train_log.jsonl
sha256sum $W/v2_smoke_s0/best_model.pt $W/v2_smoke_s0/last_model.pt | sed "s|$P/||" > $D/checkpoints.sha256
echo "=== 4: smoke checks $(date -u '+%T UTC')"
CUDA_VISIBLE_DEVICES="" python experiments/exp1/v2_smoke.py --log_dir $W/v2_smoke_s0 --cache_val $CV --out $O/smoke_report.json > $O/smoke.log 2>&1 || { tail -40 $O/smoke.log; exit 1; }
cp $O/smoke_report.json $O/smoke.log $D/
cat $O/smoke_report.json
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
