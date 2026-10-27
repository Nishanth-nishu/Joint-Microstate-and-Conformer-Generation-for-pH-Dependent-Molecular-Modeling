#!/bin/bash
#SBATCH --job-name=jmcg_tpilot_s4
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=16:00:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# Resume of job 42624 (cancelled 2026-10-27 14:46 UTC to free CPUs): stage 3 output kept, rerun stages 4-6 unchanged.
# Tautomer pilot (CPU only): T1 zero-shot Uni-pKa on Tautobase + T2 head refit. PREREG: research/experiments/T_pilot/PREREG.md
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage
echo "job $SLURM_JOB_ID on $(hostname) start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
echo "=== STAGE 1: sync + deps"
rsync -a $STAGE/jmcg/ $P/jmcg/; rsync -a $STAGE/slurm_jmcg/ $P/slurm_jmcg/
rsync -a --exclude LoQI/.git $STAGE/ext/ $P/ext/
mkdir -p $P/models && rsync -a /home2/nishanth.r/jmcg_models/unipka/ $P/models/unipka/
cd $P; source venv/bin/activate
python -c 'import torch, unipkainfer, scipy; print("torch", torch.__version__, torch.cuda.is_available())'
export PYTHONPATH=$P OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
OUT=$P/results/t_pilot; D=ext/unipka_orig/dataset; TB=ext/tautobase/Tautobase_SMIRKS.txt; M=models/unipka
test -s $OUT/zeroshot/pairs.csv && echo "stage 3 output present: $(wc -l < $OUT/zeroshot/pairs.csv) lines"
echo "=== STAGE 4: feature cache (fold 0)"
python -m jmcg.experiments.taut_features --data $D --tautobase $TB --model_dir $M --out $OUT/feat
echo "=== STAGE 5: head refits"
python -m jmcg.experiments.taut_headfit --features $OUT/feat/features_fold0.npz --head $OUT/feat/head_released_fold0.pt \
  --data $D --tautobase $TB --out $OUT/fit
echo "=== STAGE 6: pre-registered analysis"
python -m jmcg.experiments.taut_analyze --stage0 $OUT/zeroshot/pairs.csv --fit $OUT/fit --dwar $D/dwar-iBond.tsv --out $OUT/analysis.json
mkdir -p $STAGE/research/experiments/T_pilot/run1
cp $OUT/analysis.json $OUT/zeroshot/pairs.csv $OUT/fit/*.csv $OUT/fit/fit_info.json $STAGE/research/experiments/T_pilot/run1/
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
