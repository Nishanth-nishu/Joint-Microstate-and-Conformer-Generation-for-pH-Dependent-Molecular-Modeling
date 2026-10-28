#!/bin/bash
#SBATCH --job-name=jmcg_cgfid
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# Wrapper-fidelity check (2026-10-28): TD generate_confs.py vs our TDGenerator on N GEOM-QM9 test molecules,
# qm9_default pinned, both scored with TD evaluate_confs.py. Rule pre-registered in jmcg/condgen/fidelity.py.
# Usage: sbatch --export=ALL,N=3 cg_fidelity.sbatch  (smoke)   then N=50.
N=${N:-3}
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage; TD=/scratch/nishanth.r/tordiff; R=$TD/torsional-diffusion
echo "job $SLURM_JOB_ID on $(hostname) N=$N start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
rsync -a $STAGE/jmcg/condgen/ $P/jmcg/condgen/; rsync -a $STAGE/slurm_jmcg/ $P/slurm_jmcg/
source $TD/venv/bin/activate
export PYTHONPATH=$P TD_REPO=$R CUDA_VISIBLE_DEVICES="" OMP_NUM_THREADS=8
MD=$TD/workdir/pretrained/qm9_default
echo "checkpoint $MD/best_model.pt sha256 $(sha256sum $MD/best_model.pt | cut -c1-16)..."
ls -la $R/data/QM9/test_smiles.csv $R/data/QM9/test_mols.pkl
echo "=== TD get_seed (provenance: how TD builds the Mol)"; grep -n -A12 "^def get_seed" $R/diffusion/sampling.py || true
OUT=$P/results/cond_gen/fidelity_n$N; mkdir -p $OUT
cd $P; python -m jmcg.condgen.fidelity subset --test_csv $R/data/QM9/test_smiles.csv --n $N --seed 0 --out $OUT/subset.csv
for s in 0 1; do
  echo "=== R$s: TD generate_confs.py seed $s"; t0=$(date +%s)
  (cd $R && python generate_confs.py --model_dir $MD --test_csv $OUT/subset.csv --out $OUT/R$s.pkl --seed $s \
     --no_energy --inference_steps 20 > $OUT/R$s.gen.log 2>&1) || echo "R$s generation FAILED"
  echo "R$s gen seconds $(( $(date +%s) - t0 ))"; tail -2 $OUT/R$s.gen.log
  echo "=== W$s: wrapper seed $s"; t0=$(date +%s)
  python -m jmcg.condgen.fidelity wrapper --model_dir $MD --subset $OUT/subset.csv --seed $s --out $OUT/W$s.pkl \
     > $OUT/W$s.gen.log 2>&1 || echo "W$s generation FAILED"
  echo "W$s gen seconds $(( $(date +%s) - t0 ))"; tail -1 $OUT/W$s.gen.log
done
for k in R0 R1 W0 W1; do
  (cd $R && python evaluate_confs.py --confs $OUT/$k.pkl --test_csv $OUT/subset.csv --true_mols data/QM9/test_mols.pkl \
     --dataset qm9 --n_workers 8 --out_results $OUT/$k.eval.pkl > $OUT/$k.eval.log 2>&1) || echo "$k eval FAILED"
  grep -E "^SUMMARY|model failure|poor ground" $OUT/$k.eval.log || true
done
echo "=== compare"
python -m jmcg.condgen.fidelity compare --R0 $OUT/R0.eval.pkl --R1 $OUT/R1.eval.pkl --W0 $OUT/W0.eval.pkl \
  --W1 $OUT/W1.eval.pkl --out $OUT/compare.json
D=$STAGE/research/experiments/CG_fidelity/n$N; mkdir -p $D
rsync -a --include '*.json' --include '*.csv' --include '*.log' --exclude '*' $OUT/ $D/
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
