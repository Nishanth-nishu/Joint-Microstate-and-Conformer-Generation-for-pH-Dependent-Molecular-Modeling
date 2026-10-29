#!/bin/bash
#SBATCH --job-name=jmcg_e2fit
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# E2/E3 (PREREG.md + amendment 2): recompute features with the MMFF-copy fix, fit the pre-registered arms, E6 strata on
# the full E1 outputs, recompute the E1 summary's N1 reproduction with the (record, target) key, re-run the E1 conformer
# comparison with NaN-RMSD molecules dropped. jmcg venv for all but the conformer comparison (TD venv).
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage
echo "job $SLURM_JOB_ID on $(hostname) start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
rsync -a $STAGE/jmcg/ $P/jmcg/; rsync -a $STAGE/data_manifests/ $P/data_manifests/; rsync -a $STAGE/research/experiments/E2_joint/PREREG.md $P/research/experiments/E2_joint/
cd $P; source venv/bin/activate; export PYTHONPATH=$P OMP_NUM_THREADS=1
echo "PREREG sha $(sha256sum research/experiments/E2_joint/PREREG.md | cut -c1-16)"
F=$P/results/e2/features_v2; FIT=$P/results/e2/fit; mkdir -p $F $FIT
t0=$(date +%s); python -m jmcg.joint.e2_tautomer features --data data_manifests/canonical_v1 --out $F --procs 8
echo "features seconds $(( $(date +%s) - t0 ))"
python - <<PY
import pandas as pd
F = pd.read_csv('$F/features.csv'); f = F.failures.fillna('')
print('states', len(F), 'with phi', F.phi1.notna().sum(), {k: int(f.str.contains(k).sum()) for k in ['rt_bonds', 'bonds_angles', 'rt_isomeric', 'stereo', 'connectivity', 'clash']})
PY
python -m jmcg.joint.e2_tautomer fit --data data_manifests/canonical_v1 --features $F/features.csv \
  --energies results/seq/e1/full/energies.csv --out $FIT
mkdir -p results/seq/e6; python -m jmcg.seq.e6_strata --data data_manifests/canonical_v1 --e1 results/seq/e1/full --out results/seq/e6 > results/seq/e6/stdout.txt 2>&1 || echo 'E6 FAILED'
python - <<'PY'
import json, pandas as pd
X1 = pd.read_csv('results/seq/e1/full/e1_external.csv'); X1 = X1[X1.template == X1.template.iloc[0]]
n1 = pd.read_csv('results/n1_pilot/full_v2/transitions.csv'); n1['source_record'] = [f'{s}:{i}' for s, i in zip(n1['set'], n1['idx'])]
n1['tk'] = n1.target.round(3)
x = X1.assign(tk=X1.target.round(3), set=X1.source_record.str.split(':').str[0])
m = x.merge(n1[['source_record', 'tk', 'A3D-mean']], on=['source_record', 'tk'], how='inner')
rest = x[~x.obs_id.isin(m.obs_id)]
m2 = rest.merge(n1.drop_duplicates(['set', 'tk'], keep=False)[['set', 'tk', 'A3D-mean']], on=['set', 'tk'], how='inner')
m = pd.concat([m, m2])
r = dict(n_matched=len(m), n_matched_by_set_target=len(m2), n_external=len(X1), max_abs_diff=float((m.oracle - m['A3D-mean']).abs().max()))
print('reproduction_vs_N1 (record,target key)', r)
json.dump(r, open('results/seq/e1/full/reproduction_vs_N1_v2.json', 'w'), indent=1)
PY
deactivate
O=$P/results/e1_confs/n200
(source /scratch/nishanth.r/tordiff/venv/bin/activate && python -m jmcg.condgen.e1_confs compare --ref $O/td0.eval.pkl $O/td1.eval.pkl \
  --arm $O/etkdg0.eval.pkl $O/etkdg1.eval.pkl $O/etkdgmmff0.eval.pkl --out $O/compare_v2.json > /dev/null && echo "compare_v2 ok")
D=$STAGE/research/experiments/E2_joint; mkdir -p $D/fit $D/features_v2
cp $F/features.csv $D/features_v2/; cp $FIT/* $D/fit/
mkdir -p $STAGE/research/experiments/M3_seq/e6; cp results/seq/e6/* $STAGE/research/experiments/M3_seq/e6/
cp results/seq/e1/full/reproduction_vs_N1_v2.json $STAGE/research/experiments/M3_seq/e1/full/
cp $O/compare_v2.json $STAGE/research/experiments/M3_seq/e1_confs/
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
