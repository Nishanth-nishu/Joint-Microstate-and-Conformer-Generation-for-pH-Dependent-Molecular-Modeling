#!/bin/bash
#SBATCH --job-name=jmcg_exp1m0
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# EXP-1 M0 (revised, EXP-1 §10.3): oracle-velocity integration test of the V2 geometry sampler on the VALIDATION paired
# cache. CPU only, no learning; gates any FlexiTors training. TD venv (torch_geometric to unpickle the cache); the TD
# data under /scratch/nishanth.r/tordiff is read, never written.
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage; TD=/scratch/nishanth.r/tordiff
echo "job $SLURM_JOB_ID on $(hostname) start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
rsync -a --checksum $STAGE/flexitors/ $P/flexitors/; rsync -a --checksum $STAGE/experiments/exp1/ $P/experiments/exp1/
md5sum $P/flexitors/geometry/projector.py $P/flexitors/geometry/l_step.py $P/experiments/exp1/m0_oracle.py $P/flexitors/tests/test_projector.py | sed "s|$P/||"
source $TD/venv/bin/activate; cd $P; export FLEXITORS_TD_ROOT=$TD/torsional-diffusion PYTHONPATH=$P OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 CUDA_VISIBLE_DEVICES=""
echo "=== tests"
python -m pytest -q flexitors/tests/test_projector.py flexitors/tests/test_internal.py flexitors/tests/test_torsion_types.py flexitors/tests/test_no_dead_code.py 2>&1 | tail -15
CACHE=$TD/data/QM9/cache_paired.val; sha256sum $CACHE | cut -c1-16
O=$P/results/exp1/m0; mkdir -p $O
echo "=== smoke (first 100 molecules)"
python experiments/exp1/m0_oracle.py --cache $CACHE --out $O/smoke --procs 16 --limit 100 --ablation_pairs 20 > $O/smoke.log 2>&1 || { tail -30 $O/smoke.log; exit 1; }
python -c "import json; r=json.load(open('$O/smoke/m0_summary.json')); print(r['verdict'], {k: (v['n_pairs'], round(v['recovery'] or 0, 3), round(v['seconds_per_pair_median'], 4)) for k, v in r['strata'].items()})"
echo "=== full validation cache $(date -u '+%T UTC')"
python experiments/exp1/m0_oracle.py --cache $CACHE --out $O/full --procs 16 --ablation_pairs 300 > $O/full.log 2>&1 || { tail -30 $O/full.log; exit 1; }
python - <<PY
import json; r = json.load(open('$O/full/m0_summary.json'))
print('VERDICT', r['verdict']); print(json.dumps(r['strata'], indent=1)); print(json.dumps(r['ablation_restore_off'], indent=1))
PY
D=$STAGE/research/experiments/EXP1_M0; mkdir -p $D
cp $O/full/m0_summary.json $D/; cp $O/smoke/m0_summary.json $D/m0_summary_smoke100.json; gzip -c $O/full/m0_pairs.csv > $D/m0_pairs.csv.gz
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
