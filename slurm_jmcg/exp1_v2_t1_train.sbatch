#!/bin/bash
#SBATCH --job-name=jmcg_v2_t1
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --cpus-per-task=14
#SBATCH --mem=20G
#SBATCH --time=2-12:00:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# EXP1-V2-T1 training, one arm (ARM=proj|noproj). Design pre-declared in research/experiments/EXP1_V2_train/DESIGN.md.
# CPU, reduced budget (10 epochs); GPU quota exhausted. TD data/checkout read only. Resumable: resubmit with RESUME=1.
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage; TD=/scratch/nishanth.r/tordiff
echo "job $SLURM_JOB_ID ARM=$ARM on $(hostname) start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
case $ARM in proj) GP=on ;; noproj) GP=off ;; *) echo "bad ARM"; exit 1 ;; esac
rsync -a --checksum $STAGE/flexitors/ $P/flexitors/; rsync -a --checksum $STAGE/experiments/exp1/ $P/experiments/exp1/
(cd $P && md5sum flexitors/training/*.py flexitors/models/*.py flexitors/sampling/*.py)
source $TD/venv/bin/activate; cd $P
export FLEXITORS_TD_ROOT=$TD/torsional-diffusion PYTHONPATH=$P CUDA_VISIBLE_DEVICES=""
CT=$TD/data/QM9/cache_paired.train; CV=$TD/data/QM9/cache_paired.val; W=$P/workdir/exp1/T1_${ARM}_s0
D=$STAGE/research/experiments/EXP1_V2_train/T1_${ARM}_s0; mkdir -p $D
OMP_NUM_THREADS=4 python -m pytest -q flexitors/tests 2>&1 | tail -1
R=""; [ "${RESUME:-0}" = 1 ] && R="--resume" || rm -rf $W
OMP_NUM_THREADS=10 MKL_NUM_THREADS=10 python -m flexitors.training.loop --cache_train $CT --cache_val $CV --log_dir $W \
  --w_tau 0 --w_geo 1 --geo_norm base --l_head_init zero --geo_projector $GP --sigma_b 0.15 --seed 0 --lr 1e-3 \
  --batch_size 32 --n_epochs 10 --limit_val_iters 100 --num_workers 4 --log_every 100 $R > $D/train_$SLURM_JOB_ID.log 2>&1 &
TP=$!
while kill -0 $TP 2>/dev/null; do  # mirror logs to staging every 10 minutes
  cp $W/{epochs.jsonl,train_log.jsonl,model_parameters.yml} $D/ 2>/dev/null || true; sleep 600
done
rc=0; wait $TP || rc=$?
cp $W/{epochs.jsonl,train_log.jsonl,model_parameters.yml} $D/ || true
sha256sum $W/{best,last}_model.pt | sed "s|$P/||" > $D/checkpoints.sha256 || true
echo "train exit $rc end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
exit $rc
