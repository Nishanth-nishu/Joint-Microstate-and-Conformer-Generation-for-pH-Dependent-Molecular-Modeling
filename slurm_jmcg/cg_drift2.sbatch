#!/bin/bash
#SBATCH --job-name=jmcg_cgdrift2
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# Drift diagnosis step 2: drift starts in model call 0 on identical input (job 45212). Per-module outputs of call 0,
# two processes, to locate the first operation whose output differs.
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage; TD=/scratch/nishanth.r/tordiff; R=$TD/torsional-diffusion
echo "job $SLURM_JOB_ID on $(hostname) start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
rsync -a $STAGE/jmcg/condgen/ $P/jmcg/condgen/
source $TD/venv/bin/activate
cd $P
export PYTHONPATH=$P TD_REPO=$R CUDA_VISIBLE_DEVICES="" OMP_NUM_THREADS=4
MD=$TD/workdir/pretrained/qm9_default; OUT=$P/results/cond_gen/drift2; mkdir -p $OUT
pip list 2>/dev/null | grep -iE "^(torch|torch-scatter|torch-sparse|torch-cluster|torch-geometric|e3nn|numpy) " || true
echo "=== TD score model: graph construction + aggregation calls"
grep -nE "radius|scatter|knn|def forward|def build|tensor_product|o3\.|jit|set\b|sorted|dict\(" $R/diffusion/score_model.py | head -60
T="python -m jmcg.condgen.drift_trace trace --model_dir $MD --state C/C(O)=C/C(C)=O --seed 0 --hook_modules"
for cfg in "default" "hash0 PYTHONHASHSEED=0"; do
  set -- $cfg; name=$1; shift
  for rep in x y; do env "$@" $T --out $OUT/${name}_$rep.pkl > $OUT/${name}_$rep.log 2>&1 || echo "$name $rep FAILED: $(tail -3 $OUT/${name}_$rep.log)"; done
  echo "=== $name ($*)"; python -m jmcg.condgen.drift_trace compare $OUT/${name}_x.pkl $OUT/${name}_y.pkl
done
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
