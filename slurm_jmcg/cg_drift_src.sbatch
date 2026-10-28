#!/bin/bash
#SBATCH --job-name=jmcg_cgdriftsrc
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:05:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# Read-only: TD score-model forward source + the layer list traced by job 45218.
TD=/scratch/nishanth.r/tordiff; R=$TD/torsional-diffusion; source $TD/venv/bin/activate
sed -n '1,60p;120,215p' $R/diffusion/score_model.py
cd /scratch/nishanth.r/jmcg && PYTHONPATH=. python - <<'PY'
import pickle
A = pickle.load(open('results/cond_gen/drift2/default_x.pkl', 'rb'))
for n, ty, o, i in A['layer_out']:
    print('LAYER', n, ty, [None if x is None else x.shape for x in o], [None if x is None else x.shape for x in i])
PY
