#!/bin/bash
#SBATCH --job-name=jmcg_step2
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# Step 2: jmcg venv (py3.10, offline wheels from /home2) + build the stratified pKa/tautomer set + push.
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage
echo "job $SLURM_JOB_ID on $(hostname) start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
echo "=== STAGE 1: sync code + data"
rsync -a $STAGE/jmcg/ $P/jmcg/; rsync -a $STAGE/slurm_jmcg/ $P/slurm_jmcg/
mkdir -p $P/ext && rsync -a --exclude LoQI/.git $STAGE/ext/ $P/ext/
echo "=== STAGE 2: venv"
cd $P
if [[ ! -x venv/bin/python ]]; then
  PY=$(command -v python3.10 || command -v python3); echo "base python $PY $($PY --version)"
  $PY -m venv venv
fi
source venv/bin/activate
python --version
pip install -q --no-index --find-links /home2/nishanth.r/wheels_jmcg rdkit numpy pandas scipy scikit-learn tqdm pyyaml dimorphite_dl
python -c "import rdkit, dimorphite_dl, numpy, pandas; print('rdkit', rdkit.__version__, 'numpy', numpy.__version__)"
pip freeze > requirements.jmcg.lock.txt
echo "=== STAGE 3: build stratified set"
PYTHONPATH=$P python -m jmcg.data.build_sensitive_set --ext $P/ext --out $P/data_manifests/sensitive_set_v2
mkdir -p $STAGE/data_manifests && rsync -a $P/data_manifests/ $STAGE/data_manifests/
cp requirements.jmcg.lock.txt $STAGE/
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
