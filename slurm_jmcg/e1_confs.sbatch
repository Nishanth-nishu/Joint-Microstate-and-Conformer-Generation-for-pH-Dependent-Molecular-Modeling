#!/bin/bash
#SBATCH --job-name=jmcg_e1confs
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# E1 conformer component (2026-10-28): TD wrapper vs ETKDGv3 (no refinement) on 200 GEOM-QM9 test molecules,
# seeds 0/1, budget 2 x n_true, scored with TD evaluate_confs.py; ETKDGv3+MMFF94s as a labelled ablation (seed 0).
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage; TD=/scratch/nishanth.r/tordiff; R=$TD/torsional-diffusion
echo "job $SLURM_JOB_ID on $(hostname) start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
rsync -a $STAGE/jmcg/condgen/ $P/jmcg/condgen/
source $TD/venv/bin/activate
export PYTHONPATH=$P TD_REPO=$R CUDA_VISIBLE_DEVICES="" OMP_NUM_THREADS=8
MD=$TD/workdir/pretrained/qm9_default; OUT=$P/results/e1_confs/n200; mkdir -p $OUT
sha256sum $MD/best_model.pt | cut -c1-16
cd $P; python -m jmcg.condgen.fidelity subset --test_csv $R/data/QM9/test_smiles.csv --n 200 --seed 0 --out $OUT/subset.csv
for s in 0 1; do
  t0=$(date +%s); python -m jmcg.condgen.fidelity wrapper --model_dir $MD --subset $OUT/subset.csv --seed $s --out $OUT/td$s.pkl > $OUT/td$s.gen.log 2>&1
  echo "td$s seconds $(( $(date +%s) - t0 )) $(tail -1 $OUT/td$s.gen.log)"
  t0=$(date +%s); python -m jmcg.condgen.e1_confs etkdg --subset $OUT/subset.csv --seed $s --out $OUT/etkdg$s.pkl
  echo "etkdg$s seconds $(( $(date +%s) - t0 ))"
done
python -m jmcg.condgen.e1_confs etkdg --subset $OUT/subset.csv --seed 0 --mmff --out $OUT/etkdgmmff0.pkl
for k in td0 td1 etkdg0 etkdg1 etkdgmmff0; do
  (cd $R && python evaluate_confs.py --confs $OUT/$k.pkl --test_csv $OUT/subset.csv --true_mols data/QM9/test_mols.pkl \
     --dataset qm9 --n_workers 8 --out_results $OUT/$k.eval.pkl > $OUT/$k.eval.log 2>&1) || echo "$k eval FAILED"
  grep -E "^SUMMARY" $OUT/$k.eval.log
done
python -m jmcg.condgen.e1_confs compare --ref $OUT/td0.eval.pkl $OUT/td1.eval.pkl \
  --arm $OUT/etkdg0.eval.pkl $OUT/etkdg1.eval.pkl $OUT/etkdgmmff0.eval.pkl --out $OUT/compare.json > /dev/null
python - <<PY
import json; r = json.load(open("$OUT/compare.json")); print('n_common', r['n_common'], 'failures', r['failures'])
for a, v in r['arms'].items(): print(a.split('/')[-1], {k: round(x, 3) for k, x in v.items()})
for a, v in r['paired_arm_minus_ref'].items(): print('minus TD', a.split('/')[-1], {k: (round(x['mean'], 3), [round(c, 3) for c in x['ci95']]) for k, x in v.items() if k in ('MAT_R', 'COV_R@0.5', 'COV_R@0.25', 'MAT_P', 'COV_P@0.5')})
PY
D=$STAGE/research/experiments/M3_seq/e1_confs; mkdir -p $D
rsync -a --include '*.json' --include '*.csv' --include '*.log' --exclude '*' $OUT/ $D/
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
