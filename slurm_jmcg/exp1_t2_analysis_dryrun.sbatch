#!/bin/bash
#SBATCH --job-name=jmcg_t2_anadry
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=01:00:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# Dry run of t2_final_analysis.py after addendum C (strata, near-duplicate, variability), T1 outputs standing in for
# the six T2 read-outs (geometry stand-ins copied from geom_validate_seed{s}.json). No new test-set number.
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage; TD=/scratch/nishanth.r/tordiff
set -euo pipefail
source $TD/venv/bin/activate; cd $P; export PYTHONPATH=$P
rsync -a --checksum $STAGE/experiments/exp1/ $P/experiments/exp1/; rsync -a $STAGE/jmcg/condgen/ $P/jmcg/condgen/
IN=$P/results/e1_confs/n200; T1=$P/results/exp1/T1_e1; T=$P/workdir/tmp_t2anadry_$SLURM_JOB_ID; mkdir -p $T; trap 'rm -rf $T' EXIT
for s in 0 1; do for k in 0 1 2; do for m in proj noproj; do
  for e in pkl eval.pkl pkl.log.json; do ln -s $T1/td${s}_L${m}fix.$e $T/td${s}_T2_${m}_s$k.$e; done
done; done
python - $STAGE/research/experiments/EXP1_V2_T2/geom_validate_seed$s.json $T/geom_validate_seed$s.json <<'PY'
import json, sys
g = json.load(open(sys.argv[1]))
for k in (0, 1, 2):
    for m in ('proj', 'noproj'):
        g['per_molecule'][f'{m}_s{k}'] = g['per_molecule'][f'L{m}fix']
        g['summary']['conformer_checks'][f'{m}_s{k}'] = g['summary']['conformer_checks'][f'L{m}fix']
json.dump(g, open(sys.argv[2], 'w'))
PY
done
python experiments/exp1/t2_final_analysis.py --n200 $IN --t1 $T1 --t2 $T --geom_dir $T --train_dir $P/workdir/exp1 --out $T/fa.json | tail -12
python - $T/fa.json <<'PY'
import json, sys
r = json.load(open(sys.argv[1]))
for s in ('seed0', 'seed1'):
    e = r['e1'][s]
    print(s, 'n_attempted', e['n_attempted'], 'n_common', e['n_common'])
    for k, v in e['strata'].items():
        print('  ', k, v['n'], {x: round(y['mean'], 4) for x, y in v.items() if x != 'n'})
    w = e['without_near_duplicate']; print('  near-dup removed', w['removed'], 'n', w['n'], round(w['proj_s0-td']['mean'], 4))
    print('  variability', e['variability_MAT_R'])
print({k: (v['epochs'], v['best_epoch']) for k, v in r['training'].items()})
PY
