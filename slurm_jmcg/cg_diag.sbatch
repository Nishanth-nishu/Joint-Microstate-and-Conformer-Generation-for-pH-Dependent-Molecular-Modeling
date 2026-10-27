#!/bin/bash
#SBATCH --job-name=jmcg_cgdiag
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# Conditional generation, smoke 1 follow-up (2026-10-28): no new states, no new checkpoint.
#  a) backend tests incl. the new H-bond clash test; b) re-validate the smoke-1 SDFs with the H-bond-aware clash rule
#  (no regeneration); c) diagnose cross-process TD drift (smoke 1: same seed, separate processes, 1e-4..2.7e-2 A).
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage; TD=/scratch/nishanth.r/tordiff
echo "job $SLURM_JOB_ID on $(hostname) start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
rsync -a $STAGE/jmcg/condgen/ $P/jmcg/condgen/; rsync -a $STAGE/slurm_jmcg/ $P/slurm_jmcg/
source $TD/venv/bin/activate
cd $P
export PYTHONPATH=$P TD_REPO=$TD/torsional-diffusion CUDA_VISIBLE_DEVICES=""
QM9_DIR=$TD/workdir/pretrained/qm9_1order  # same checkpoint as smoke 1 (sha256 4199663e...)
S1=$P/results/cond_gen/smoke1; OUT=$P/results/cond_gen/smoke1_diag; mkdir -p $OUT/repro
echo "=== a: backend tests"
OMP_NUM_THREADS=4 python -m jmcg.condgen.test_backend --model_dir $QM9_DIR --states_csv jmcg/condgen/smoke_states.csv \
  2>&1 | grep -E '^(PASS|FAIL|test_backend)' || echo "TEST_BACKEND FAILED"
echo "=== b: re-validate smoke-1 SDFs (H-bond-aware clash rule)"
OMP_NUM_THREADS=4 python -m jmcg.condgen.revalidate --run_dir $S1/qm9 --out $OUT/revalidated_qm9
OMP_NUM_THREADS=4 python -m jmcg.condgen.revalidate --run_dir $S1/rep_seed1 --out $OUT/revalidated_rep_seed1
echo "=== c: cross-process reproducibility, 2 processes per setting"
for cfg in "t4_hashrand 4 random" "t4_hash0 4 0" "t1_hashrand 1 random" "t1_hash0 1 0"; do
  set -- $cfg
  for rep in x y; do
    OMP_NUM_THREADS=$2 PYTHONHASHSEED=$3 python -m jmcg.condgen.run --model_dir $QM9_DIR --out $OUT/repro/$1_$rep \
      --states_csv $S1/repro_states.csv --n_confs 8 --seed 0 > $OUT/repro/$1_$rep.log 2>&1 || echo "run $1_$rep failed"
  done
done
python - <<EOF
import glob, os, numpy as np
from rdkit import Chem
def xyz(p): return [m.GetConformer().GetPositions() for m in Chem.SDMolSupplier(p, removeHs=False)]
for cfg in ['t4_hashrand', 't4_hash0', 't1_hashrand', 't1_hash0']:
    for f in sorted(glob.glob("$OUT/repro/%s_x/*__td.sdf" % cfg)):
        A, B = xyz(f), xyz(f.replace('_x/', '_y/'))
        d = max(float(np.abs(a - b).max()) for a, b in zip(A, B))
        print(f"{cfg:12s} {os.path.basename(f):32s} max|dx| across processes = {d:.1e} A")
EOF
echo "=== d: copy to staging"
D=$STAGE/research/experiments/CG_smoke/run1_diag; mkdir -p $D
rsync -a --include '*/' --include '*.json' --include '*.csv' --include '*.log' --exclude '*' $OUT/ $D/
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
