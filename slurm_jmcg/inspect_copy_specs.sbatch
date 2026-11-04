#!/bin/bash
#SBATCH --job-name=jmcg_inspect3
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# Read-only: copy the FlexiTors / round-5 design documents and the round-2 S4 code diff to /home2 for reading (not staged).
P=/scratch/nishanth.r/jmcg; T=/scratch/nishanth.r/tordiff; D=/home2/nishanth.r/jmcg_ref; mkdir -p $D
cp $P/round5/specs/EXP-1.md $P/round5/specs/EVAL.md $P/round5/PLAN.md $P/round5/DECISION.md $P/round5/review_P1_m2.md $P/round5/review_P2_m2.md $P/EXPERIMENTS.md $D/ 2>&1
cp $T/results/round5/STATUS.md $D/round5_STATUS.md; cp $T/results/round5/gpu_ledger.tsv $D/
ls $T/results/round5/gate0; cat $T/results/round5/gate0/*summary* 2>/dev/null | head -40
ls $P/round2 | head; ls $T/torsional-diffusion; grep -rn "lambda_embed\|lam_cond\|lamcond" $T/torsional-diffusion --include=*.py | head -30
cat $T/workdir/qm9_S4_lamcond_e100_s0/model_parameters.yml
ls $T/data/QM9; ls $T/data/QM9/qm9 | head -3; ls $T/data/QM9/qm9 | wc -l
ls $T/logs | grep -i "trainabl\|r2_train\|train" | head -20
f=$(ls $T/logs/*trainabl_1036_0.log 2>/dev/null | head -1); grep -n "Epoch\|epoch" $f | head -5; grep -n "Epoch\|epoch" $f | tail -3
nvidia-smi 2>&1 | head -5
