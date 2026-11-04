#!/bin/bash
#SBATCH --job-name=jmcg_p1resume
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=04:00:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# Resume of P1 job 46561: TD chunk 1 was OOM-killed (16 parallel TD processes in 48G); compare failed (KeyError psi1:
# extract/td_features did not write the 2D psi columns). Re-runs TD chunk 1 alone, then td_features, compare, sampler.
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage; TD=/scratch/nishanth.r/tordiff
echo "job $SLURM_JOB_ID on $(hostname) start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
rsync -a --checksum $STAGE/jmcg/ $P/jmcg/; md5sum $P/jmcg/joint/p1_unmin.py | sed "s|$P/||"
cd $P; export PYTHONPATH=$P OMP_NUM_THREADS=1; O=$P/results/p1_unmin
echo "chunk 1 log tail before rerun:"; { grep -E "^P[0-9]" $O/td_chunk_1.log || true; } | tail -2
QM9_DIR=$TD/workdir/pretrained/qm9_default
(source $TD/venv/bin/activate; export TD_REPO=$TD/torsional-diffusion CUDA_VISIBLE_DEVICES=""
 /usr/bin/time -v python -m jmcg.condgen.run --model_dir $QM9_DIR --states_csv $O/td_states_1.csv --n_confs 30 --seed 0 \
   --out $O/td > $O/td_chunk_1_rerun.log 2> $O/td_chunk_1_rerun.time) || { echo "chunk 1 rerun FAILED"; tail -20 $O/td_chunk_1_rerun.log; grep "Maximum resident" $O/td_chunk_1_rerun.time; }
grep "Maximum resident" $O/td_chunk_1_rerun.time || true
echo "sdf files: $(ls $O/td/*__td.sdf | wc -l)"
source venv/bin/activate
python -m jmcg.joint.p1_unmin td_features --panel $O/panel.csv --td_dir $O/td --out $O/features_td.csv
python -m jmcg.joint.p1_unmin compare --panel $O/panel.csv --mmff $O/features_mmff.csv --unmin $O/features_unmin.csv \
  --td $O/features_td.csv --data data_manifests/canonical_v1 --energies results/seq/e1/full/energies.csv \
  --e2_features results/e2/features_v2/features.csv --out $O/compare > $O/compare.log 2>&1 || { tail -30 $O/compare.log; exit 1; }
for arm in S1 J1; do
  python -m jmcg.joint.sample --state_dir results/joint/smoke_states --coef $O/compare/coef_unmin.json --arm $arm \
    --ph 7.4 --out $O/sample_unmin_$arm --n_pairs 20 --seed 0 --feature_path unmin > $O/sample_unmin_$arm.log 2>&1 || { tail -20 $O/sample_unmin_$arm.log; exit 1; }
  echo "== $arm"; cat $O/sample_unmin_$arm/run_info.json; tail -25 $O/sample_unmin_$arm.log; done
D=$STAGE/research/experiments/P1_unmin; mkdir -p $D
cp $O/panel.csv $O/panel_info.json $O/features_*.csv $O/compare.log $D/; cp -r $O/compare $D/
for arm in S1 J1; do mkdir -p $D/sample_unmin_$arm; cp $O/sample_unmin_$arm/* $D/sample_unmin_$arm/; done
cat $O/td_chunk_*.log | grep -E "^P[0-9]" | sort -u > $D/td_per_state.txt || true
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
