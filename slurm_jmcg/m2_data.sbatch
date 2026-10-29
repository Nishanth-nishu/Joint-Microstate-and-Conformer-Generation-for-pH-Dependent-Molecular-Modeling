#!/bin/bash
#SBATCH --job-name=jmcg_m2data
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# Milestone 2 (2026-10-28): canonical data layer v1 (identity tests, build, audit with leakage checks). CPU, jmcg venv.
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage
echo "job $SLURM_JOB_ID on $(hostname) start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
rsync -a $STAGE/jmcg/ $P/jmcg/; rsync -a $STAGE/slurm_jmcg/ $P/slurm_jmcg/; rsync -a $STAGE/data_manifests/ $P/data_manifests/
cd $P; source venv/bin/activate; export PYTHONPATH=$P
echo "jmcg HEAD $(git rev-parse --short HEAD); python $(python -V 2>&1); rdkit $(python -c 'import rdkit; print(rdkit.__version__)')"
echo "=== tests"; python -m jmcg.data.test_canonical
OUT=$P/data_manifests/canonical_v1
echo "=== build"; t0=$(date +%s); python -m jmcg.data.build_canonical --ext $P/ext --pka_curated data_manifests/sensitive_set_v2/pka_curated.csv --out $OUT
echo "build seconds $(( $(date +%s) - t0 ))"
echo "=== audit"; python -m jmcg.data.audit --data $OUT --out $OUT/audit
ls -la $OUT; rsync -a $OUT/ $STAGE/data_manifests/canonical_v1/
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
