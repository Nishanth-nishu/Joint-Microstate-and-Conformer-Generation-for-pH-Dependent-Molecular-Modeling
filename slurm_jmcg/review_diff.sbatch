#!/bin/bash
#SBATCH --job-name=jmcg_review_diff
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# Read only for the repo: copy it (with .git) to a temp dir, apply sync_push.sbatch stage 1 there, `git add -A` in the
# copy, and write exactly what a commit would contain to $STAGE/logs/review_<job>.{stat,patch}.
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage; H=/scratch/nishanth.r/tmp_review_$SLURM_JOB_ID/jmcg
set -euo pipefail
mkdir -p $H; trap 'rm -rf /scratch/nishanth.r/tmp_review_$SLURM_JOB_ID' EXIT
rsync -a --exclude '*.pt' --exclude '*.pkl' --exclude '*.npz' --exclude '*.npy' --exclude venv/ $P/ $H/
rsync -a --exclude '*.pdf' --exclude 'papers/*.txt' $STAGE/research/ $H/research/
rsync -a $STAGE/jmcg/ $H/jmcg/; rsync -a $STAGE/slurm_jmcg/ $H/slurm_jmcg/; cp $STAGE/README.md $H/README.md
[[ -d $STAGE/data_manifests ]] && rsync -a $STAGE/data_manifests/ $H/data_manifests/
cp $STAGE/*.sbatch $H/slurm_jmcg/
rsync -a --checksum $STAGE/flexitors/ $H/flexitors/; rsync -a --checksum $STAGE/experiments/exp1/ $H/experiments/exp1/
cd $H; git add -A
O=$STAGE/logs/review_$SLURM_JOB_ID
git diff --cached --stat=200 HEAD > $O.stat
git diff --cached HEAD -- . ':(exclude)*.log' ':(exclude)*.log.json' ':(exclude)*.jsonl' ':(exclude)*.csv.gz' > $O.patch
git diff --cached --name-status HEAD > $O.names
du -ch $(git diff --cached --name-only --diff-filter=A HEAD) 2>/dev/null | tail -1
echo "files $(wc -l < $O.names)"; cut -f1 $O.names | sort | uniq -c
