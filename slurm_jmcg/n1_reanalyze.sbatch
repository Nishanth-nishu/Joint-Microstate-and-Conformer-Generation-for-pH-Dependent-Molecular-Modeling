#!/bin/bash
#SBATCH --job-name=jmcg_n1reana
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# N1 pilot run 2 re-analysis with the sanity gate referenced to Uni-pKa Table 1 columns (PREREG amendment 3). No model calls.
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage
echo "job $SLURM_JOB_ID on $(hostname) start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
echo "=== STAGE 1: sync code"
rsync -a $STAGE/jmcg/ $P/jmcg/; rsync -a $STAGE/slurm_jmcg/ $P/slurm_jmcg/
cd $P; source venv/bin/activate; export PYTHONPATH=$P
echo "=== STAGE 2: analysis (gate v2: Table 1 references)"
python -m jmcg.experiments.n1_analyze --run results/n1_pilot/full_v2 --out_name analysis_gate_table1.json
cp results/n1_pilot/full_v2/analysis_gate_table1.json $STAGE/research/experiments/N1_pilot/run2/
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
