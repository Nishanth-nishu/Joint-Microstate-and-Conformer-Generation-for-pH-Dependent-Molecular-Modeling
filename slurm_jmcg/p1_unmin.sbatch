#!/bin/bash
#SBATCH --job-name=jmcg_p1unmin
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# P1 (EXPLORATORY): M5 feature-extraction blocker. Arms mmff (control, E2 definition) / unmin (ETKDGv3, no
# minimisation) / td (TD qm9_default, E/Z lock default) on one state panel; validator thresholds unchanged.
# CPU only. jmcg venv except the TD generation step (TD venv, frozen checkout read-only).
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage; TD=/scratch/nishanth.r/tordiff
echo "job $SLURM_JOB_ID on $(hostname) start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
rsync -a --checksum $STAGE/jmcg/ $P/jmcg/; md5sum $P/jmcg/joint/*.py $P/jmcg/condgen/validate.py | sed "s|$P/||"
cd $P; export PYTHONPATH=$P OMP_NUM_THREADS=1
O=$P/results/p1_unmin; mkdir -p $O
source venv/bin/activate
echo "=== tests"; python -m jmcg.joint.test_joint
echo "=== panel"; python -m jmcg.joint.p1_unmin panel --data data_manifests/canonical_v1 --e1 results/seq/e1/full \
  --smoke $STAGE/research/experiments/M3_seq/smoke/states --out $O --per_cat 150 --seed 0 --td_chunks 16
cat $O/panel_info.json
for path in mmff unmin; do t0=$(date +%s)
  python -m jmcg.joint.p1_unmin extract --panel $O/panel.csv --path $path --out $O/features_$path.csv --procs 16
  echo "$path seconds $(( $(date +%s) - t0 ))"; done
deactivate
echo "=== TD generation (16 chunks, 1 thread each) $(date -u '+%T UTC')"
QM9_DIR=$TD/workdir/pretrained/qm9_default
sha256sum $QM9_DIR/best_model.pt | grep -q "^6a9a6ee9b0e51bda" || { echo "CHECKPOINT SHA MISMATCH"; exit 1; }
(source $TD/venv/bin/activate; export TD_REPO=$TD/torsional-diffusion CUDA_VISIBLE_DEVICES=""
 for k in $(seq 0 15); do
   python -m jmcg.condgen.run --model_dir $QM9_DIR --states_csv $O/td_states_$k.csv --n_confs 30 --seed 0 \
     --out $O/td > $O/td_chunk_$k.log 2>&1 &
 done; wait) || echo "some TD chunk returned non-zero"
grep -L "^done" $O/td_chunk_*.log || true
echo "TD done $(date -u '+%T UTC'); sdf files: $(ls $O/td/*__td.sdf | wc -l)"
source venv/bin/activate
python -m jmcg.joint.p1_unmin td_features --panel $O/panel.csv --td_dir $O/td --out $O/features_td.csv
echo "=== compare"
python -m jmcg.joint.p1_unmin compare --panel $O/panel.csv --mmff $O/features_mmff.csv --unmin $O/features_unmin.csv \
  --td $O/features_td.csv --data data_manifests/canonical_v1 --energies results/seq/e1/full/energies.csv \
  --e2_features results/e2/features_v2/features.csv --out $O/compare > $O/compare.log 2>&1 || { tail -30 $O/compare.log; exit 1; }
head -c 6000 $O/compare/p1_result.json; echo
echo "=== M5 sampler, unmin path (smoke parents, pH 7.4; coefficients from the exploratory unmin refit)"
for arm in S1 J1; do
  python -m jmcg.joint.sample --state_dir results/joint/smoke_states --coef $O/compare/coef_unmin.json --arm $arm \
    --ph 7.4 --out $O/sample_unmin_$arm --n_pairs 20 --seed 0 --feature_path unmin > $O/sample_unmin_$arm.log 2>&1 || { tail -20 $O/sample_unmin_$arm.log; exit 1; }
  echo "== $arm"; cat $O/sample_unmin_$arm/run_info.json; done
D=$STAGE/research/experiments/P1_unmin; mkdir -p $D
cp $O/panel.csv $O/panel_info.json $O/features_*.csv $O/compare.log $D/; cp -r $O/compare $D/
for arm in S1 J1; do mkdir -p $D/sample_unmin_$arm; cp $O/sample_unmin_$arm/* $D/sample_unmin_$arm/; done
cp $O/td/summary.csv $D/td_summary.csv 2>/dev/null || true
cat $O/td_chunk_*.log | grep -E "^P[0-9]" > $D/td_per_state.txt || true
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
