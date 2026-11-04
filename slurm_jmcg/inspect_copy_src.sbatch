#!/bin/bash
#SBATCH --job-name=jmcg_inspect4
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# Read-only: copy flexitors (jmcg repo copy) and the TD files EXP-1 reuses to /home2/nishanth.r/jmcg_ref/src for reading;
# describe the paired cache format. Nothing under /scratch is modified.
P=/scratch/nishanth.r/jmcg; T=/scratch/nishanth.r/tordiff; D=/home2/nishanth.r/jmcg_ref/src; mkdir -p $D/td
rsync -a --exclude __pycache__ $P/flexitors $D/
for f in diffusion/likelihood.py diffusion/score_model.py diffusion/sampling.py diffusion/torus.py utils/torsion.py utils/dataset.py utils/training.py utils/standardization.py utils/utils.py utils/parsing.py utils/featurization.py train.py; do mkdir -p $D/td/$(dirname $f); cp $T/torsional-diffusion/$f $D/td/$f; done
ls $P/tools | head -50; cp $P/tools/build_paired_pickles.py $D/ 2>/dev/null
ls -la $T/data/QM9/cache_paired.val $T/data/QM9/cache_paired.train $T/data/QM9/split.npy
source $T/venv/bin/activate; cd $T/torsional-diffusion
python - <<'PY'
import pickle, torch, numpy as np
d = pickle.load(open('/scratch/nishanth.r/tordiff/data/QM9/cache_paired.val', 'rb'))
print(type(d), len(d))
x = d[0]; print(x)
for k in x.keys: 
    v = x[k]; print(k, type(v), getattr(v, 'shape', None) if not isinstance(v, list) else (len(v), type(v[0]), getattr(v[0], 'shape', None)))
print('pair_ok' in x.keys, [k for k in x.keys])
PY
