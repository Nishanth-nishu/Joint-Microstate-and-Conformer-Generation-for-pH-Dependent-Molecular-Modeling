#!/bin/bash
#SBATCH --job-name=jmcg_push
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# Usage: sbatch --export=ALL,MSG="commit message" sync_push.sbatch
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage
echo "job $SLURM_JOB_ID on $(hostname) start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
echo "=== STAGE 1: sync staging -> scratch"
rsync -a --exclude '*.pdf' --exclude 'papers/*.txt' $STAGE/research/ $P/research/
[[ -d $STAGE/jmcg ]] && rsync -a $STAGE/jmcg/ $P/jmcg/
[[ -d $STAGE/slurm_jmcg ]] && rsync -a $STAGE/slurm_jmcg/ $P/slurm_jmcg/
cp $STAGE/*.sbatch $P/slurm_jmcg/ 2>/dev/null || { mkdir -p $P/slurm_jmcg && cp $STAGE/*.sbatch $P/slurm_jmcg/; }
mkdir -p $P/logs && cp $STAGE/logs/*.log $P/logs/ 2>/dev/null || true
echo "=== STAGE 2: commit + push"
cd $P
grep -qx 'research/papers/' .gitignore || printf 'research/papers/\nresearch/scholar_runs/*/text/\n' >> .gitignore
git add -A
git diff --cached --quiet && echo "nothing to commit" || git commit -q -m "${MSG:-research update}"
git log --format='%h %an %s' -3
git push origin main 2>&1 | tail -2
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
