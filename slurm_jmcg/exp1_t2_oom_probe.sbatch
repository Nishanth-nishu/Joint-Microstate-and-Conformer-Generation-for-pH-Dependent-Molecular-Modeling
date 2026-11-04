#!/bin/bash
#SBATCH --job-name=jmcg_t2_oomprobe
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=02:00:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# EXP1-V2-T2 task 1.7: where do the oom_kill events of 47249 / 49899 come from, and do they change any output?
# Re-runs the seed-0 steps of 47249 with the same resources (16 CPUs, 40G), writing to a separate directory, while a
# monitor logs the cgroup oom_kill counter, cgroup memory and the process list every 2 s with phase markers.
# Then compares every output with the saved one (bit-identical RMSD matrices / coordinates expected).
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage; TD=/scratch/nishanth.r/tordiff; R=$TD/torsional-diffusion
echo "job $SLURM_JOB_ID on $(hostname) start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
rsync -a --checksum $STAGE/flexitors/ $P/flexitors/; rsync -a --checksum $STAGE/experiments/exp1/ $P/experiments/exp1/
source $TD/venv/bin/activate; cd $P
export FLEXITORS_TD_ROOT=$R PYTHONPATH=$P CUDA_VISIBLE_DEVICES="" OMP_NUM_THREADS=1
IN=$P/results/e1_confs/n200; OLD=$P/results/exp1/T1_e1; OUT=$P/results/exp1/T2_oomprobe; W=$P/workdir/exp1; mkdir -p $OUT
D=$STAGE/research/experiments/EXP1_V2_T2/oom_probe; mkdir -p $D
CG=$(awk -F: '$2=="" || $2 ~ /memory/ {print $3; exit}' /proc/self/cgroup)
for f in /sys/fs/cgroup$CG/memory.events /sys/fs/cgroup/memory$CG/memory.oom_control; do [ -f $f ] && EV=$f && break; done
for f in /sys/fs/cgroup$CG/memory.current /sys/fs/cgroup/memory$CG/memory.usage_in_bytes; do [ -f $f ] && CUR=$f && break; done
echo "cgroup $CG events ${EV:-none} current ${CUR:-none}"
PH=$OUT/phase; echo start > $PH
( while true; do
    k=$(grep -E '^oom_kill ' ${EV} | awk '{print $2}'); m=$(( $(cat ${CUR}) / 1048576 ))
    echo "$(date -u +%T) phase=$(cat $PH) oom_kill=$k mem_MB=$m nproc=$(ps --no-headers -u $USER -o pid | wc -l)"
    ps --no-headers -u $USER -o rss=,comm= --sort=-rss | head -3 | awk '{printf "   top %s MB %s\n", int($1/1024), $2}'
    sleep 2
  done ) > $D/monitor.log 2>&1 &
MON=$!
trap 'kill $MON 2>/dev/null || true' EXIT
run() { echo "$1" > $PH; echo "=== $1 $(date -u +%T)"; shift; "$@"; }
run readout_Lproj  python experiments/exp1/t1_e1_readout.py --confs $IN/td0.pkl --method l --run_dir $W/T1_proj_s0 --project on --workers 16 --out $OUT/td0_Lproj.pkl > $OUT/r1.log 2>&1
run readout_mmff   python experiments/exp1/t1_e1_readout.py --confs $IN/td0.pkl --method mmff --workers 16 --out $OUT/td0_mmff.pkl > $OUT/r2.log 2>&1
run eval_Lproj  bash -c "cd $R && python evaluate_confs.py --confs $OUT/td0_Lproj.pkl --test_csv $IN/subset.csv --true_mols data/QM9/test_mols.pkl --dataset qm9 --n_workers 16 --out_results $OUT/td0_Lproj.eval.pkl > $OUT/td0_Lproj.eval.log 2>&1"
run eval_mmff   bash -c "cd $R && python evaluate_confs.py --confs $OUT/td0_mmff.pkl --test_csv $IN/subset.csv --true_mols data/QM9/test_mols.pkl --dataset qm9 --n_workers 16 --out_results $OUT/td0_mmff.eval.pkl > $OUT/td0_mmff.eval.log 2>&1"
run compare python -m jmcg.condgen.e1_confs compare --ref $IN/td0.eval.pkl --arm $OUT/td0_Lproj.eval.pkl $OUT/td0_mmff.eval.pkl --out $OUT/compare_seed0.json > /dev/null
echo done > $PH; sleep 3
python - <<PY
import json, pickle, numpy as np
# the current code includes the rigid-molecule fix, so the L read-out is compared with the saved *fix* arm (49899)
for k, old in (('td0_Lproj', 'td0_Lprojfix'), ('td0_mmff', 'td0_mmff')):
    a, b = pickle.load(open('$OUT/%s.pkl' % k, 'rb')), pickle.load(open('$OLD/%s.pkl' % old, 'rb'))
    dx = max(float(np.abs(x.GetConformer().GetPositions() - y.GetConformer().GetPositions()).max()) for s in a for x, y in zip(a[s], b[s]))
    ea, eb = pickle.load(open('$OUT/%s.eval.pkl' % k, 'rb'))['results'], pickle.load(open('$OLD/%s.eval.pkl' % old, 'rb'))['results']
    same_keys = set(ea) == set(eb)
    dr = max(float(np.nanmax(np.abs(ea[m]['rmsd'] - eb[m]['rmsd']))) if not np.isnan(ea[m]['rmsd']).all() else 0.0 for m in ea)
    print(k, 'max coordinate diff vs saved', dx, '| eval molecules identical', same_keys, '| max RMSD diff vs saved', dr)
new = json.load(open('$OUT/compare_seed0.json'))['paired_arm_minus_ref']['$OUT/td0_Lproj.eval.pkl']['MAT_R']
print('MAT_R Lproj - TD rerun', new)
PY
echo "oom_kill per phase:"; awk '{for(i=1;i<=NF;i++){if($i~/^phase=/)p=$i; if($i~/^oom_kill=/)k=$i}} /phase=/{last[p]=k} END{for(p in last) print p, last[p]}' $D/monitor.log | sort -t= -k3 -n
grep -h "Killed\|MemoryError\|oom" $OUT/*.log || echo "no kill/memory message in step logs"
cp $OUT/*.log $D/ 2>/dev/null || true
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
