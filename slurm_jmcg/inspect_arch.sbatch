#!/bin/bash
#SBATCH --job-name=jmcg_inspect
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:10:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# Read-only inspection: which conformer-model code exists beyond frozen TD; TD git state; data/caches/checkpoints.
P=/scratch/nishanth.r/jmcg; T=/scratch/nishanth.r/tordiff
echo "job $SLURM_JOB_ID $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
cd $P; git log --oneline | head -5; ls; ls prior_work | head -50
echo "== nn.Module / torch in jmcg repo (excl venv,results)"
grep -rln --include=*.py -E "nn\.Module|import torch" --exclude-dir=venv --exclude-dir=results --exclude-dir=ext . | head -40
echo "== TD tree"; cd $T; ls; git log --oneline | head -15; git status --short | head -30; git diff --stat HEAD | tail -5
echo "== TD workdir dirs"; ls -d */ ; ls workdir 2>/dev/null | head -50
for d in workdir/*/; do echo "$d: $(ls $d | head -8 | tr '\n' ' ')"; done 2>/dev/null | head -60
echo "== data"; ls data 2>/dev/null | head; du -sh data/* 2>/dev/null | head -20
echo "== handoff"; ls $P/prior_work; sed -n 1,200p $P/prior_work/TORDIFF_NAX_HANDOFF.md 2>/dev/null
nvidia-smi --query-gpu=index,memory.used,memory.total,utilization.gpu --format=csv
