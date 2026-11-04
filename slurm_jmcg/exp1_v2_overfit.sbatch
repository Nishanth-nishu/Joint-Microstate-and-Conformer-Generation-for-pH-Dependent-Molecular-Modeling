#!/bin/bash
#SBATCH --job-name=jmcg_v2_overfit
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# EXP-1 V2 smoke S2 failed (job 46744/46848: val geo 0.01302 vs Δ̂=0 base 0.01288). Diagnostic, declared before the run
# (2026-10-29): overfit 32 training molecules (batch 8, 100 epochs = 400 iterations, CPU), arms w_tau ∈ {0, 1}, w_geo 1.
# Capacity/plumbing OK iff mean train geo / geo_base over the last 10 epochs < 0.9. Training molecules only; no
# validation/test molecule informs anything (validation capped at 1 batch, reported only).
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage; TD=/scratch/nishanth.r/tordiff
echo "job $SLURM_JOB_ID on $(hostname) start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
rsync -a --checksum $STAGE/flexitors/ $P/flexitors/; rsync -a --checksum $STAGE/experiments/exp1/ $P/experiments/exp1/
source $TD/venv/bin/activate; cd $P
export FLEXITORS_TD_ROOT=$TD/torsional-diffusion PYTHONPATH=$P OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 CUDA_VISIBLE_DEVICES=""
CT=$TD/data/QM9/cache_paired.train; CV=$TD/data/QM9/cache_paired.val; W=$P/workdir/exp1
D=$STAGE/research/experiments/EXP1_V2_smoke/overfit; mkdir -p $D
for wt in 0 1; do
  rm -rf $W/v2_overfit_wt$wt
  python -m flexitors.training.loop --cache_train $CT --cache_val $CV --log_dir $W/v2_overfit_wt$wt --limit_train_mols 32 \
    --batch_size 8 --n_epochs 100 --limit_val_iters 1 --w_tau $wt --w_geo 1 --sigma_b 0.15 --seed 0 --num_workers 0 \
    --log_every 1000 > $D/train_wt$wt.log 2>&1 &
done
wait
for wt in 0 1; do cp $W/v2_overfit_wt$wt/epochs.jsonl $D/epochs_wt$wt.jsonl; done
python - <<PY
import json, numpy as np
for wt in (0, 1):
    e = [json.loads(x) for x in open('$D/epochs_wt%d.jsonl' % wt)]
    r = [x['train_geo'] / x['train_geo_base'] for x in e]
    v = np.mean(r[-10:])
    print('w_tau', wt, 'epochs', len(e), 'ratio first10 %.3f last10 %.3f' % (np.mean(r[:10]), v), 'OVERFIT_OK' if v < 0.9 else 'OVERFIT_FAIL',
          'td/base last10 %.3f' % np.mean([x['train_td'] / x['train_td_base'] for x in e[-10:]]))
    print('  ratio by epoch', [round(x, 3) for x in r[::10]])
PY
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
