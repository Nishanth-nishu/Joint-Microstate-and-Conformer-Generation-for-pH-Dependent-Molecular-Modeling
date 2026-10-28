#!/bin/bash
#SBATCH --job-name=jmcg_cgdrift3
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# Drift diagnosis step 3: verify that TD's torus score-norm table differs between processes, and that importing TD
# under a fixed numpy seed (td_backend.TORUS_IMPORT_SEED) makes same-seed generation bit-identical across processes.
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage; TD=/scratch/nishanth.r/tordiff; R=$TD/torsional-diffusion
echo "job $SLURM_JOB_ID on $(hostname) start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
rsync -a $STAGE/jmcg/condgen/ $P/jmcg/condgen/; rsync -a $STAGE/slurm_jmcg/ $P/slurm_jmcg/
source $TD/venv/bin/activate
cd $P
export PYTHONPATH=$P TD_REPO=$R CUDA_VISIBLE_DEVICES="" OMP_NUM_THREADS=4
MD=$TD/workdir/pretrained/qm9_default; OUT=$P/results/cond_gen/drift3; mkdir -p $OUT
echo "=== TD diffusion/torus.py (module level)"; sed -n '1,60p' $R/diffusion/torus.py
ls -la $P/.p.npy $P/.score.npy 2>&1 || true
H='import hashlib, numpy as np; f=lambda a: hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()[:16]'
echo "=== a: plain TD import, 2 processes (expect different score_norm_ if the table is the cause)"
for r in 1 2; do python -c "$H
import sys; sys.path.insert(0, '$R'); import diffusion.torus as t
print('plain', 'score_norm_', f(t.score_norm_), 'p_', f(t.p_), 'score_', f(t.score_))"; done
echo "=== b: import via td_backend (seeded), 2 processes"
for r in 1 2; do python -c "$H
from jmcg.condgen import td_backend as tb; t = tb.td_torus
print('seeded', 'score_norm_', f(t.score_norm_), 'p_', f(t.p_), 'score_', f(t.score_), 'np state restored:', np.random.get_state()[1][0])"; done
echo "=== c: drift trace with the fix, 2 processes + 2 calls in one process"
T="python -m jmcg.condgen.drift_trace trace --model_dir $MD --state C/C(O)=C/C(C)=O --seed 0"
$T --out $OUT/fix_x.pkl > /dev/null 2>&1; $T --out $OUT/fix_y.pkl > /dev/null 2>&1
python -m jmcg.condgen.drift_trace compare $OUT/fix_x.pkl $OUT/fix_y.pkl
OMP_NUM_THREADS=1 $T --out $OUT/fix_t1.pkl > /dev/null 2>&1
echo "--- 4 threads vs 1 thread (with the fix)"; python -m jmcg.condgen.drift_trace compare $OUT/fix_x.pkl $OUT/fix_t1.pkl
$T --repeat 2 --out $OUT/fix_same > /dev/null 2>&1
echo "--- within one process"; python -m jmcg.condgen.drift_trace compare $OUT/fix_same.0 $OUT/fix_same.1
echo "=== d: backend tests"
python -m jmcg.condgen.test_backend --model_dir $MD --states_csv jmcg/condgen/smoke_states.csv 2>&1 | grep -E '^(PASS|FAIL|test_backend)'
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
