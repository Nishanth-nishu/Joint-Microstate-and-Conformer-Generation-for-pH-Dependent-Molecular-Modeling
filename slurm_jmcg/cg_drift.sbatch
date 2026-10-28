#!/bin/bash
#SBATCH --job-name=jmcg_cgdrift
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:45:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# Cross-process TD drift diagnosis (2026-10-28): acac Z-enol, seed 0, qm9_default (smoke 3: 0.49 A between processes).
# Two processes per environment setting; report the first stage where the traces differ.
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage; TD=/scratch/nishanth.r/tordiff
echo "job $SLURM_JOB_ID on $(hostname) start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
rsync -a $STAGE/jmcg/condgen/ $P/jmcg/condgen/; rsync -a $STAGE/slurm_jmcg/ $P/slurm_jmcg/
source $TD/venv/bin/activate
cd $P
export PYTHONPATH=$P TD_REPO=$TD/torsional-diffusion CUDA_VISIBLE_DEVICES=""
MD=$TD/workdir/pretrained/qm9_default; OUT=$P/results/cond_gen/drift; mkdir -p $OUT
lscpu | grep -E "Model name|Flags" | cut -c1-200
python -c 'import torch; print(torch.__config__.show())' | head -20
T="python -m jmcg.condgen.drift_trace trace --model_dir $MD --state C/C(O)=C/C(C)=O --seed 0"
echo "=== within one process (2 calls)"
OMP_NUM_THREADS=4 $T --repeat 2 --out $OUT/same >/dev/null 2>&1
python -m jmcg.condgen.drift_trace compare $OUT/same.0 $OUT/same.1
run_cfg () {  # name, then env assignments / prefix
  local name=$1; shift
  for rep in x y; do env "$@" $T --out $OUT/${name}_$rep.pkl > $OUT/${name}_$rep.log 2>&1 || echo "$name $rep FAILED: $(tail -2 $OUT/${name}_$rep.log)"; done
  echo "=== $name ($*)"; python -m jmcg.condgen.drift_trace compare $OUT/${name}_x.pkl $OUT/${name}_y.pkl
}
run_cfg default    OMP_NUM_THREADS=4
run_cfg threads1   OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
run_cfg mkl_cbwr   OMP_NUM_THREADS=4 MKL_CBWR=COMPATIBLE
run_cfg t1_cbwr    OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MKL_CBWR=COMPATIBLE
run_cfg aten_def   OMP_NUM_THREADS=4 ATEN_CPU_CAPABILITY=default
run_cfg no_aslr    OMP_NUM_THREADS=4 setarch $(uname -m) -R
echo "=== verbose diff, default"; python -m jmcg.condgen.drift_trace compare $OUT/default_x.pkl $OUT/default_y.pkl --verbose | head -30
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
