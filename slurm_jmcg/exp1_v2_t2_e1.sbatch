#!/bin/bash
#SBATCH --job-name=jmcg_v2_t2_e1
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --cpus-per-task=8
#SBATCH --mem=22G
#SBATCH --time=1-12:00:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# EXP1-V2-T2 part D, the pre-declared final E1 evaluation (research/experiments/EXP1_V2_T2/DESIGN.md), run ONCE:
# checkpoint completeness and hashes, E1 input manifest, the fixed read-out of each of the six best_model.pt on the
# frozen-TD samples (TD seeds 0 and 1), evaluate_confs.py, t2_geom_validate.py with every arm, t2_final_analysis.py,
# v2_diag.py on validation molecules, sacct resources of the training jobs. If a training run is not complete and a
# T2 training or resume job is still queued or running, the job resubmits itself after those and does nothing else.
# Addendum D (2026-11-04): 8 CPUs / 22 G / 6 workers instead of 32 / 90 G / 32 (node shared, 2 CPUs left free);
# outputs are deterministic, so only wall time changes; cost is re-timed for the head and constrained xTB on the same
# 20 molecules with the same worker count (step 4b).
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage; TD=/scratch/nishanth.r/tordiff; R=$TD/torsional-diffusion
echo "job $SLURM_JOB_ID on $(hostname) start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
W=$P/workdir/exp1; IN=$P/results/e1_confs/n200; T1=$P/results/exp1/T1_e1; OUT=$P/results/exp1/T2_e1
X=$STAGE/research/experiments/EXP1_V2_T2; D=$X/e1; RUNS="proj_s0 noproj_s0 proj_s1 noproj_s1 proj_s2 noproj_s2"
source $TD/venv/bin/activate

# 1. all six runs complete: epochs 0..19 logged, both checkpoints present and equal to the hashes the training job wrote
incomplete=""
for r in $RUNS; do
  python - $W/T2_$r/epochs.jsonl <<'PY' && [ -f $W/T2_$r/best_model.pt ] || incomplete="$incomplete $r"
import json, sys
e = sorted({json.loads(x)['epoch'] for x in open(sys.argv[1])})
sys.exit(0 if e == list(range(20)) else 1)
PY
done
if [ -n "$incomplete" ]; then
  q=$(squeue -u $USER -h -o "%i %j" | awk '$2 ~ /^jmcg_v2_t2_(proj|noproj)_s[012](_resume)?$/ {print $1}' | paste -sd:)
  if [ -n "$q" ]; then
    echo "incomplete:$incomplete; T2 jobs still queued/running ($q): resubmitting after them, no evaluation now"
    sbatch --dependency=afterany:$q $STAGE/slurm_jmcg/exp1_v2_t2_e1.sbatch; exit 0
  fi
  echo "incomplete:$incomplete and no T2 job left; not evaluating"; exit 1
fi
for r in $RUNS; do
  (cd $P && sha256sum -c $X/train/T2_$r/checkpoints.sha256) || { echo "checkpoint hash mismatch $r"; exit 1; }
done

# 2. E1 inputs are the audited ones (manifest written by exp1_t2_manifest.sbatch)
sha256sum -c $X/e1_manifest.sha256 || { echo "E1 manifest missing or mismatch"; exit 1; }

# 3. run once
[ -e $OUT ] && { echo "$OUT exists: the final evaluation already ran (or started); refusing to run it again"; exit 1; }
mkdir -p $OUT $D; echo "job $SLURM_JOB_ID $(date -u '+%F %T UTC')" > $OUT/STARTED

source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/Modules/init/bash
rsync -a --checksum $STAGE/flexitors/ $P/flexitors/; rsync -a --checksum $STAGE/experiments/exp1/ $P/experiments/exp1/
rsync -a $STAGE/jmcg/condgen/ $P/jmcg/condgen/
cd $P; export FLEXITORS_TD_ROOT=$R PYTHONPATH=$P CUDA_VISIBLE_DEVICES="" OMP_NUM_THREADS=1
(md5sum flexitors/models/*.py experiments/exp1/t1_e1_readout.py experiments/exp1/t2_*.py) | tee $D/code.md5
OMP_NUM_THREADS=4 python -m pytest -q flexitors/tests 2>&1 | tail -1
for r in $RUNS; do sha256sum workdir/exp1/T2_$r/best_model.pt; done | tee $D/checkpoints.sha256
NW=6  # addendum D (was 32)
# peak PID FILE: wait for PID, writing its wall time, the job cgroup's peak memory and mean CPU load (cores) to FILE
CG=$(awk -F: '$2=="" || $2 ~ /memory/ {print $3; exit}' /proc/self/cgroup)
for f in /sys/fs/cgroup$CG/memory.current /sys/fs/cgroup/memory$CG/memory.usage_in_bytes; do [ -f $f ] && CUR=$f && break; done
peak() {
  local t0=$(date +%s) mx=0 m c0=$(awk '/^cpu /{print $2+$3+$4}' /proc/stat) c1 rc=0
  while kill -0 $1 2>/dev/null; do m=$(cat ${CUR:-/dev/null} 2>/dev/null || echo 0); [ "${m:-0}" -gt $mx ] && mx=$m; sleep 2; done
  wait $1 || rc=$?
  c1=$(awk '/^cpu /{print $2+$3+$4}' /proc/stat)
  echo "exit $rc wall_s $(( $(date +%s) - t0 )) cgroup_peak_GB $(awk "BEGIN{print $mx/2^30}")" \
       "node_busy_cores $(awk "BEGIN{print ($c1-$c0)/$(getconf CLK_TCK)/($(date +%s)-$t0+1)}")" | tee $2
  return $rc
}

# 4. read-out (fixed read-out, guard and torsion restoration as T1 "fix"), then TD evaluate_confs
for s in 0 1; do
  for r in $RUNS; do
    case $r in proj_*) PJ=on ;; *) PJ=off ;; esac
    python experiments/exp1/t1_e1_readout.py --confs $IN/td$s.pkl --method l --run_dir $W/T2_$r --project $PJ \
      --workers $NW --out $OUT/td${s}_T2_$r.pkl > $OUT/td${s}_T2_$r.readout.log 2>&1 &
    peak $! $OUT/td${s}_T2_$r.time; tail -3 $OUT/td${s}_T2_$r.readout.log
  done
done
for s in 0 1; do
  for r in $RUNS; do
    k=td${s}_T2_$r
    (cd $R && python evaluate_confs.py --confs $OUT/$k.pkl --test_csv $IN/subset.csv --true_mols data/QM9/test_mols.pkl \
       --dataset qm9 --n_workers $NW --out_results $OUT/$k.eval.pkl > $OUT/$k.eval.log 2>&1) || echo "$k eval FAILED"
  done
done

# 4b. addendum D: like-for-like cost, head (T2 proj s0) vs constrained xTB, first 20 molecules of TD seed 0, NW workers
set +u; module load u22/xtb/6.7.1; set -u
export TMPDIR=/scratch/nishanth.r/tmp_xtb_$SLURM_JOB_ID; mkdir -p $TMPDIR
python experiments/exp1/t1_e1_readout.py --confs $IN/td0.pkl --method l --run_dir $W/T2_proj_s0 --project on \
  --limit_mols 20 --workers $NW --out $OUT/timing_Lproj.pkl > $OUT/timing_Lproj.log 2>&1 &
peak $! $OUT/timing_Lproj.time
python experiments/exp1/t1_e1_readout.py --confs $IN/td0.pkl --method xtb --limit_mols 20 --workers $NW \
  --out $OUT/timing_xtb.pkl > $OUT/timing_xtb.log 2>&1 &
peak $! $OUT/timing_xtb.time
rm -rf $TMPDIR

# 5. geometry validation with every arm, the analysis, validation diagnostics
for s in 0 1; do
  A=""; for r in $RUNS; do A="$A --arm ${r}=$OUT/td${s}_T2_$r"; done
  python experiments/exp1/t2_geom_validate.py --in_dir $IN --t1_dir $T1 --true_mols $R/data/QM9/test_mols.pkl --seed $s \
    --workers $NW $A --out $OUT/geom_validate_seed$s.json > $OUT/geom_validate_seed$s.log
done
python experiments/exp1/t2_final_analysis.py --n200 $IN --t1 $T1 --t2 $OUT --geom_dir $OUT --train_dir $W \
  --out $OUT/final_analysis.json | tee $OUT/final_analysis.log
OMP_NUM_THREADS=8 python experiments/exp1/v2_diag.py --cache_val $TD/data/QM9/cache_paired.val \
  --runs $(for r in $RUNS; do echo -n "$W/T2_$r "; done) --n 1000 --out $OUT/v2_diag_val.json > $OUT/v2_diag_val.log 2>&1 \
  || echo "v2_diag FAILED"

# 6. training resources
J=$(awk '{print $4","$6}' $X/train_jobs.txt | paste -sd,)
sacct -j $J -o JobID%12,JobName%28,State,Start,End,Elapsed,TotalCPU,AllocCPUS,MaxRSS,ExitCode -P > $OUT/train_sacct.txt
cat $OUT/train_sacct.txt

rsync -a --include '*.json' --include '*.log' --include '*.time' --include '*.txt' --include 'STARTED' --exclude '*' $OUT/ $D/
echo "job $SLURM_JOB_ID $(date -u '+%F %T UTC')" > $OUT/DONE; cp $OUT/DONE $D/
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
