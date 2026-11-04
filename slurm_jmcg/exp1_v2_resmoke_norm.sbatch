#!/bin/bash
#SBATCH --job-name=jmcg_v2_resmoke
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# EXP-1 V2 re-smoke after S2 FAIL (jobs 46744/46848). Declared before the run (2026-10-29):
#   ONE change to the objective: --geo_norm base (optimise L_geo / its Δ̂ = 0 base, detached, per batch).
#   Everything else as 46744: w_tau 1, w_geo 1, sigma_b 0.15, seed 0, batch 32, lr 1e-3, 400 iterations, 60 val batches.
#   Criteria S1–S5 of experiments/exp1/v2_smoke.py UNCHANGED; S2 still compares RAW val L_geo with 0.9 × raw base.
#   Infrastructure change (no effect on the objective): DataLoader workers use 1 torch thread (oversubscription fix).
#   Secondary, reported only: overfit diagnostic (32 mols, w_tau 1, geo_norm base) vs 46855's 0.913.
# CPU (plafnet2 GPU quota exhausted). TD data and checkout read only.
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage; TD=/scratch/nishanth.r/tordiff
echo "job $SLURM_JOB_ID on $(hostname) start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
rsync -a --checksum $STAGE/flexitors/ $P/flexitors/; rsync -a --checksum $STAGE/experiments/exp1/ $P/experiments/exp1/
(cd $P && md5sum flexitors/training/*.py flexitors/models/*.py experiments/exp1/v2_smoke.py)
source $TD/venv/bin/activate; cd $P
export FLEXITORS_TD_ROOT=$TD/torsional-diffusion PYTHONPATH=$P CUDA_VISIBLE_DEVICES=""
CT=$TD/data/QM9/cache_paired.train; CV=$TD/data/QM9/cache_paired.val
W=$P/workdir/exp1; O=$P/results/exp1/v2_resmoke_norm; D=$STAGE/research/experiments/EXP1_V2_smoke/resmoke_norm
mkdir -p $W $O $D
echo "=== 1: tests $(date -u '+%T UTC')"
OMP_NUM_THREADS=4 python -m pytest -q flexitors/tests 2>&1 | tail -6
echo "=== 2: training (400 iterations, geo_norm base) $(date -u '+%T UTC')"
rm -rf $W/v2_resmoke_norm_s0; t0=$(date +%s)
OMP_NUM_THREADS=12 MKL_NUM_THREADS=12 python -m flexitors.training.loop --cache_train $CT --cache_val $CV --log_dir $W/v2_resmoke_norm_s0 \
  --sigma_b 0.15 --w_tau 1 --w_geo 1 --geo_norm base --n_epochs 1 --seed 0 --num_workers 4 --log_every 20 \
  --limit_train_iters 400 --limit_val_iters 60 > $O/train.log 2>&1 || { tail -40 $O/train.log; exit 1; }
echo "train wall seconds $(( $(date +%s) - t0 ))"; grep EPOCH $O/train.log
echo "=== 3: smoke checks $(date -u '+%T UTC')"
OMP_NUM_THREADS=12 python experiments/exp1/v2_smoke.py --log_dir $W/v2_resmoke_norm_s0 --cache_val $CV --out $O/smoke_report.json \
  > $O/smoke.log 2>&1 || { tail -40 $O/smoke.log; exit 1; }
cat $O/smoke_report.json
cp $O/{smoke_report.json,train.log} $D/; cp $W/v2_resmoke_norm_s0/{epochs.jsonl,train_log.jsonl,model_parameters.yml} $D/
sha256sum $W/v2_resmoke_norm_s0/{best,last}_model.pt | sed "s|$P/||" > $D/checkpoints.sha256
echo "=== 4: overfit diagnostic, w_tau 1, geo_norm base (secondary) $(date -u '+%T UTC')"
rm -rf $W/v2_overfit_wt1_norm
OMP_NUM_THREADS=12 python -m flexitors.training.loop --cache_train $CT --cache_val $CV --log_dir $W/v2_overfit_wt1_norm --limit_train_mols 32 \
  --batch_size 8 --n_epochs 100 --limit_val_iters 1 --w_tau 1 --w_geo 1 --geo_norm base --sigma_b 0.15 --seed 0 --num_workers 0 \
  --log_every 1000 > $D/overfit_wt1_norm.log 2>&1
cp $W/v2_overfit_wt1_norm/epochs.jsonl $D/overfit_epochs_wt1_norm.jsonl
python - <<PY
import json, numpy as np
e = [json.loads(x) for x in open('$D/overfit_epochs_wt1_norm.jsonl')]
r = [x['train_geo'] / x['train_geo_base'] for x in e]
print('overfit w_tau 1 geo_norm base: ratio first10 %.3f last10 %.3f (46855 raw: 0.913)' % (np.mean(r[:10]), np.mean(r[-10:])),
      'td/base last10 %.3f' % np.mean([x['train_td'] / x['train_td_base'] for x in e[-10:]]))
print('  ratio by epoch', [round(x, 3) for x in r[::10]])
PY
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
