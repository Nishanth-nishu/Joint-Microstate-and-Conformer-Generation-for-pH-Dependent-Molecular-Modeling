#!/bin/bash
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --cpus-per-task=8
#SBATCH --mem=20G
#SBATCH --time=4-00:00:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# EXP1-V2-T2 part D, one run (ARM=proj|noproj, SEED=0|1|2); design: research/experiments/EXP1_V2_T2/DESIGN.md.
# 20 epochs, flags identical to T1. SEED=0 continues T1's run from its last_model.pt in a copied T2 directory (T1 is
# never written). Resubmit with RESUME=1 after a failure. CPU unless the job was submitted with a GPU.
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage; TD=/scratch/nishanth.r/tordiff
echo "job $SLURM_JOB_ID ARM=$ARM SEED=$SEED on $(hostname) cpus $SLURM_CPUS_PER_TASK start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
case $ARM in proj) GP=on ;; noproj) GP=off ;; *) echo "bad ARM"; exit 1 ;; esac
case $SEED in 0|1|2) ;; *) echo "bad SEED"; exit 1 ;; esac
rsync -a --checksum $STAGE/flexitors/ $P/flexitors/; rsync -a --checksum $STAGE/experiments/exp1/ $P/experiments/exp1/
(cd $P && md5sum flexitors/training/*.py flexitors/models/*.py)
source $TD/venv/bin/activate; cd $P
export FLEXITORS_TD_ROOT=$TD/torsional-diffusion PYTHONPATH=$P
CT=$TD/data/QM9/cache_paired.train; CV=$TD/data/QM9/cache_paired.val; W=$P/workdir/exp1/T2_${ARM}_s$SEED
D=$STAGE/research/experiments/EXP1_V2_T2/train/T2_${ARM}_s$SEED; mkdir -p $D
OMP_NUM_THREADS=4 python -m pytest -q flexitors/tests 2>&1 | tail -1
R="--resume"
if [ "${RESUME:-0}" = 1 ] && [ ! -f $W/last_model.pt ]; then echo "RESUME=1 but no $W/last_model.pt; not starting from scratch"; exit 1; fi
if [ "${RESUME:-0}" != 1 ]; then
  if [ "$SEED" = 0 ]; then
    [ -e $W ] && { echo "$W exists; use RESUME=1"; exit 1; }
    T1=$P/workdir/exp1/T1_${ARM}_s0; mkdir -p $W
    cp -p $T1/{last_model.pt,best_model.pt,epochs.jsonl,train_log.jsonl,model_parameters.yml} $W/
    sha256sum $T1/{last_model.pt,best_model.pt} $W/{last_model.pt,best_model.pt} | sed "s|$P/||"
  else
    rm -rf $W; R=""
  fi
fi
OMP_NUM_THREADS=4 MKL_NUM_THREADS=4 python -m flexitors.training.loop --cache_train $CT --cache_val $CV --log_dir $W \
  --w_tau 0 --w_geo 1 --geo_norm base --l_head_init zero --geo_projector $GP --sigma_b 0.15 --seed $SEED --lr 1e-3 \
  --batch_size 32 --n_epochs 20 --limit_val_iters 100 --num_workers 4 --log_every 100 $R > $D/train_$SLURM_JOB_ID.log 2>&1 &
TP=$!
while kill -0 $TP 2>/dev/null; do  # mirror logs to staging every 10 minutes
  cp $W/{epochs.jsonl,train_log.jsonl,model_parameters.yml} $D/ 2>/dev/null || true; sleep 600
done
rc=0; wait $TP || rc=$?
cp $W/{epochs.jsonl,train_log.jsonl,model_parameters.yml} $D/ || true
sha256sum $W/{best,last}_model.pt | sed "s|$P/||" > $D/checkpoints.sha256 || true
echo "train exit $rc end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
exit $rc
