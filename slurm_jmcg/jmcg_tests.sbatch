#!/bin/bash
#SBATCH --job-name=jmcg_tests
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# jmcg venv test suites (data layer, sequential baseline) + import/compile check of every jmcg module.
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage
echo "job $SLURM_JOB_ID start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -uo pipefail
rsync -a --checksum $STAGE/jmcg/ $P/jmcg/
cd $P; source venv/bin/activate; export PYTHONPATH=$P OMP_NUM_THREADS=1
python -m compileall -q jmcg && echo "compileall ok"
python -m jmcg.data.test_canonical 2>&1 | tail -3
python -m jmcg.seq.test_seq 2>&1 | tail -3
python -c "import jmcg.joint.e2_tautomer, jmcg.joint.sample, jmcg.seq.e6_strata, jmcg.seq.e1_states; print('imports ok')"
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
