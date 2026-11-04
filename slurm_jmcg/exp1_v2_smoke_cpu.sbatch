#!/bin/bash
#SBATCH --job-name=jmcg_v2smoke_cpu
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# CPU FALLBACK (GPU job 46725 held: plafnet2 GrpTRESMins gres/gpu quota exhausted). EXP-1 V2 build check + training smoke: tests (CPU), tiny CPU resume run (w_tau = 0),
# 1-epoch GPU run (w_tau = 1, w_geo = 1, sigma_b = 0.15, seed 0) on the round-2 paired QM9 cache, then the smoke
# checks of experiments/exp1/v2_smoke.py (validation molecules only). Checkpoints in $P/workdir/exp1/ (own path).
# TD data and checkout under /scratch/nishanth.r/tordiff are read, never written.
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage; TD=/scratch/nishanth.r/tordiff
echo "job $SLURM_JOB_ID on $(hostname) start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
rsync -a --checksum $STAGE/flexitors/ $P/flexitors/; rsync -a --checksum $STAGE/experiments/exp1/ $P/experiments/exp1/
(cd $P && md5sum flexitors/geometry/*.py flexitors/models/*.py flexitors/training/*.py flexitors/sampling/l_readout.py experiments/exp1/*.py)
source $TD/venv/bin/activate; cd $P
export FLEXITORS_TD_ROOT=$TD/torsional-diffusion PYTHONPATH=$P OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
nvidia-smi --query-gpu=index,name,memory.used,memory.total --format=csv || true
CT=$TD/data/QM9/cache_paired.train; CV=$TD/data/QM9/cache_paired.val
W=$P/workdir/exp1; O=$P/results/exp1/v2_smoke; mkdir -p $W $O
echo "=== 1: tests (CPU) $(date -u '+%T UTC')"
CUDA_VISIBLE_DEVICES="" python -m pytest -q flexitors/tests 2>&1 | tail -15
echo "=== 2: CPU resume check (64 train mols, 3 iters/epoch, w_tau 0)"
rm -rf $W/v2_cpu_resume
CUDA_VISIBLE_DEVICES="" python -m flexitors.training.loop --cache_train $CT --cache_val $CV --log_dir $W/v2_cpu_resume \
  --limit_train_mols 64 --limit_train_iters 3 --limit_val_iters 2 --n_epochs 1 --w_tau 0 --w_geo 1 --num_workers 0 --log_every 1 > $O/cpu_resume_a.log 2>&1
CUDA_VISIBLE_DEVICES="" python -m flexitors.training.loop --cache_train $CT --cache_val $CV --log_dir $W/v2_cpu_resume \
  --limit_train_mols 64 --limit_train_iters 3 --limit_val_iters 2 --n_epochs 2 --w_tau 0 --w_geo 1 --num_workers 0 --log_every 1 --resume > $O/cpu_resume_b.log 2>&1
grep RESUMED $O/cpu_resume_b.log
python - <<PY
import json
e = [json.loads(x)['epoch'] for x in open('$W/v2_cpu_resume/epochs.jsonl')]
t = [json.loads(x) for x in open('$W/v2_cpu_resume/train_log.jsonl')]
print('resume epochs', e, 'RESUME_OK' if e == [0, 1] else 'RESUME_FAIL')
print('w_tau=0 torsion-head grad norms', sorted({r['grad_torsion_head'] for r in t}), 'l_head', [round(r['grad_l_head'], 4) for r in t])
PY
echo "=== 3: CPU training, 1 epoch capped at 400 iterations, validation capped at 60 batches $(date -u '+%T UTC')"
rm -rf $W/v2_smoke_s0
t0=$(date +%s)
CUDA_VISIBLE_DEVICES="" OMP_NUM_THREADS=10 python -m flexitors.training.loop --cache_train $CT --cache_val $CV --log_dir $W/v2_smoke_s0 --sigma_b 0.15 \
  --w_tau 1 --w_geo 1 --n_epochs 1 --seed 0 --num_workers 6 --log_every 20 --limit_train_iters 400 --limit_val_iters 60 > $O/train.log 2>&1 || { tail -40 $O/train.log; exit 1; }
echo "train wall seconds $(( $(date +%s) - t0 ))"; head -3 $O/train.log; grep EPOCH $O/train.log
echo "=== 4: smoke checks $(date -u '+%T UTC')"
CUDA_VISIBLE_DEVICES="" OMP_NUM_THREADS=10 python experiments/exp1/v2_smoke.py --log_dir $W/v2_smoke_s0 --cache_val $CV --out $O/smoke_report.json > $O/smoke.log 2>&1 || { tail -40 $O/smoke.log; exit 1; }
cat $O/smoke_report.json
D=$STAGE/research/experiments/EXP1_V2_smoke; mkdir -p $D
cp $O/smoke_report.json $O/train.log $O/cpu_resume_*.log $D/; cp $W/v2_smoke_s0/{epochs.jsonl,train_log.jsonl,model_parameters.yml} $D/
sha256sum $W/v2_smoke_s0/best_model.pt $W/v2_smoke_s0/last_model.pt | sed "s|$P/||" > $D/checkpoints.sha256
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
