#!/bin/bash
#SBATCH --job-name=jmcg_v2_rszero
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# EXP-1 V2 re-smoke 2 after S2 FAIL in 46744 (raw) and 46998 (geo_norm base). Declared before the run (2026-10-29):
#   PRIMARY: one change vs 46998: --l_head_init zero (geometry head weights start at 0, Δ̂ = 0 at initialisation).
#     Else identical: geo_norm base, w_tau 1, w_geo 1, sigma_b 0.15, seed 0, batch 32, lr 1e-3, 400 iterations,
#     60 val batches. Criteria S1–S5 of experiments/exp1/v2_smoke.py UNCHANGED (S2: raw val L_geo < 0.9 × raw base).
#     The smoke verdict is this run's PASS field only.
#   SECONDARY (reported, not a gate, does not change the verdict): same configuration, 2,000 iterations.
# CPU (plafnet2 GPU quota exhausted). TD data and checkout read only.
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage; TD=/scratch/nishanth.r/tordiff
echo "job $SLURM_JOB_ID on $(hostname) start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
rsync -a --checksum $STAGE/flexitors/ $P/flexitors/; rsync -a --checksum $STAGE/experiments/exp1/ $P/experiments/exp1/
(cd $P && md5sum flexitors/training/*.py flexitors/models/*.py experiments/exp1/v2_smoke.py)
source $TD/venv/bin/activate; cd $P
export FLEXITORS_TD_ROOT=$TD/torsional-diffusion PYTHONPATH=$P CUDA_VISIBLE_DEVICES=""
CT=$TD/data/QM9/cache_paired.train; CV=$TD/data/QM9/cache_paired.val; W=$P/workdir/exp1
echo "=== 1: tests $(date -u '+%T UTC')"
OMP_NUM_THREADS=4 python -m pytest -q flexitors/tests 2>&1 | tail -3
run() {  # name iterations
  local N=$1 IT=$2 O=$P/results/exp1/$1 D=$STAGE/research/experiments/EXP1_V2_smoke/$1
  mkdir -p $O $D; rm -rf $W/$N; local t0=$(date +%s)
  echo "=== train $N ($IT iterations) $(date -u '+%T UTC')"
  OMP_NUM_THREADS=12 MKL_NUM_THREADS=12 python -m flexitors.training.loop --cache_train $CT --cache_val $CV --log_dir $W/$N \
    --sigma_b 0.15 --w_tau 1 --w_geo 1 --geo_norm base --l_head_init zero --n_epochs 1 --seed 0 --num_workers 4 \
    --log_every 20 --limit_train_iters $IT --limit_val_iters 60 > $O/train.log 2>&1 || { tail -40 $O/train.log; exit 1; }
  echo "train wall seconds $(( $(date +%s) - t0 ))"; grep EPOCH $O/train.log
  echo "=== smoke checks $N $(date -u '+%T UTC')"
  OMP_NUM_THREADS=12 python experiments/exp1/v2_smoke.py --log_dir $W/$N --cache_val $CV --out $O/smoke_report.json \
    > $O/smoke.log 2>&1 || { tail -40 $O/smoke.log; exit 1; }
  cat $O/smoke_report.json; echo
  cp $O/{smoke_report.json,train.log} $D/; cp $W/$N/{epochs.jsonl,train_log.jsonl,model_parameters.yml} $D/
  sha256sum $W/$N/{best,last}_model.pt | sed "s|$P/||" > $D/checkpoints.sha256
}
run v2_resmoke_zero_s0 400
run v2_resmoke_zero_s0_long 2000
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
