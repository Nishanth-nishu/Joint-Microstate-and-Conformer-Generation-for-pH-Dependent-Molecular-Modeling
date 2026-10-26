#!/bin/bash
#SBATCH --job-name=jmcg_n1pilot
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# N1 pilot (CPU only): Uni-pKa 2D vs 3D-conformer ablation on Uni-pKa's external test sets. PREREG: research/experiments/N1_pilot/PREREG.md
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage
echo "job $SLURM_JOB_ID on $(hostname) start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
echo "=== STAGE 1: sync + deps"
rsync -a $STAGE/jmcg/ $P/jmcg/; rsync -a $STAGE/slurm_jmcg/ $P/slurm_jmcg/
rsync -a --exclude LoQI/.git $STAGE/ext/ $P/ext/
mkdir -p $P/models && rsync -a /home2/nishanth.r/jmcg_models/unipka/ $P/models/unipka/
cd $P; source venv/bin/activate
pip install -q --no-index --find-links /home2/nishanth.r/wheels_jmcg "torch==2.5.1+cpu" unipkainfer lmdb pooch platformdirs
python -c 'import torch, unipkainfer; print("torch", torch.__version__, torch.version.cuda, torch.cuda.is_available())'
pip freeze > requirements.jmcg.lock.txt; cp requirements.jmcg.lock.txt $STAGE/
export PYTHONPATH=$P OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
OUT=$P/results/n1_pilot
echo "=== STAGE 2: smoke (2 transitions per set)"
python -m jmcg.experiments.n1_pilot --data ext/unipka_orig/dataset --model_dir models/unipka --out $OUT/smoke --threads 16 --limit 2
echo "=== STAGE 3: full pilot"
python -m jmcg.experiments.n1_pilot --data ext/unipka_orig/dataset --model_dir models/unipka --out $OUT/full --threads 16
echo "=== STAGE 4: pre-registered analysis"
python -m jmcg.experiments.n1_analyze --run $OUT/full
mkdir -p $STAGE/research/experiments/N1_pilot/results
cp $OUT/full/analysis.json $OUT/full/transitions.csv $STAGE/research/experiments/N1_pilot/results/
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
