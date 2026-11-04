#!/bin/bash
#SBATCH --job-name=jmcg_inspect2
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# Read-only: what FlexiTors code exists (models/training), round5 implementation record, Exp-1 spec, GPU state.
P=/scratch/nishanth.r/jmcg; T=/scratch/nishanth.r/tordiff
echo "job $SLURM_JOB_ID $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
cat $T/FLEXITORS_COMMIT; echo; diff -rq $P/flexitors $T/flexitors | head
cd $P/flexitors; find . -name '*.py' | xargs wc -l | sort -n | tail -40
echo "== models/"; ls -la models training geometry sampling; for f in models/*.py training/*.py; do echo "---- $f"; sed -n 1,60p $f; done
echo "== IMPLEMENTATION tail"; ls $P/round5; tail -80 $P/round5/IMPLEMENTATION.md
echo "== EXP-1 sec10"; sed -n '/## 10/,/## 11/p' $P/round5/specs/EXP-1.md | head -150
echo "== STATUS files"; find $T/results/round5 -maxdepth 2 -name 'STATUS.md' | head; ls $T/results/round5 | head -30
echo "== CTRL train timing"; ls $T/logs | grep -i ctrl | head; grep -h -i "epoch 99\|epoch: 99\|Run time\|elapsed" $T/logs/*CTRL_base* 2>/dev/null | tail -5
grep -l "CTRL_base_100ep" $T/logs/* 2>/dev/null | head -3 | while read f; do echo "-- $f"; head -3 $f; tail -3 $f; done
cat $T/workdir/qm9_CTRL_base_100ep_e100_s0/model_parameters.yml
