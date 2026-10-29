#!/bin/bash
#SBATCH --job-name=jmcg_m3smoke
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=01:30:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# Milestone 3 smoke (2026-10-28): sequential baseline. Stage 1 (jmcg venv): tests, E1 on 10 external transitions +
# 10 tautomer pairs, state stage on 6 smoke parents at pH 7.4. Stage 2 (TD venv): condgen on the selected states.
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage; TD=/scratch/nishanth.r/tordiff
echo "job $SLURM_JOB_ID on $(hostname) start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
rsync -a $STAGE/jmcg/ $P/jmcg/; rsync -a $STAGE/slurm_jmcg/ $P/slurm_jmcg/; rsync -a $STAGE/data_manifests/ $P/data_manifests/
cd $P; echo "jmcg HEAD $(git rev-parse --short HEAD) (+ uncommitted working tree)"
OUT=$P/results/seq/smoke; rm -rf $OUT; mkdir -p $OUT; M=$P/models/unipka
source venv/bin/activate; export PYTHONPATH=$P OMP_NUM_THREADS=8
echo "=== tests"; python -m jmcg.seq.test_seq
echo "=== E1 smoke"; python -m jmcg.seq.e1_states --data data_manifests/canonical_v1 --model_dir $M --out $OUT/e1 --limit 10 \
  --n1_transitions $P/results/n1_pilot/full_v2/transitions.csv --threads 8
printf 'id,smiles\nglycine,NCC(=O)O\nhistamine,NCCc1c[nH]cn1\nacac,CC(=O)CC(C)=O\nphenol,Oc1ccccc1\nibuprofen,CC(C)Cc1ccc(cc1)[C@H](C)C(=O)O\npropylamine,CCCN\n' > $OUT/parents.csv
echo "=== state stage"; python -m jmcg.seq.state_stage --parents $OUT/parents.csv --ph 7.4 --ph_grid 2,4,10,12 --model_dir $M --out $OUT/states --threads 8
cat $OUT/states/selected_states.csv; cat $OUT/states/macro_pka.csv
deactivate
echo "=== conformer stage (TD venv)"; source $TD/venv/bin/activate
export TD_REPO=$TD/torsional-diffusion CUDA_VISIBLE_DEVICES="" OMP_NUM_THREADS=4
python -m jmcg.condgen.run --model_dir $TD/workdir/pretrained/qm9_default --out $OUT/confs --states_csv $OUT/states/selected_states.csv \
  --n_confs 8 --seed 0 --stereo_policy enumerate
D=$STAGE/research/experiments/M3_seq/smoke; mkdir -p $D
rsync -a --include '*/' --include '*.json' --include '*.csv' --include '*.log' --exclude '*' $OUT/ $D/
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
