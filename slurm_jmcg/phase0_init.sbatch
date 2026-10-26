#!/bin/bash
#SBATCH --job-name=jmcg_phase0
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# Phase 0: fresh-history copy of TorDiff-NAX into /scratch/nishanth.r/jmcg, push to the new repo; inventory the baseline.
P=/scratch/nishanth.r/jmcg
STAGE=/home2/nishanth.r/jmcg_stage
NEW=https://github.com/Nishanth-nishu/Joint-Microstate-and-Conformer-Generation-for-pH-Dependent-Molecular-Modeling.git
echo "job $SLURM_JOB_ID on $(hostname) start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail

echo "=== STAGE 1: clone + strip history"
if [[ ! -d $P/.git ]]; then
  rm -rf $P/_clone && mkdir -p $P
  git clone -q https://github.com/Nishanth-nishu/TorDiff-NAX.git $P/_clone
  rm -rf $P/_clone/.git
  shopt -s dotglob; mv $P/_clone/* $P/; rmdir $P/_clone; shopt -u dotglob
  git mv -h >/dev/null 2>&1 || true
  mkdir -p $P/prior_work && mv $P/README.md $P/prior_work/TORDIFF_NAX_HANDOFF.md
  cp $STAGE/README.md $P/README.md
  mkdir -p $P/logs && echo "logs/" >> $P/.gitignore && echo "venv/" >> $P/.gitignore
  cd $P && git init -q -b main
fi
cd $P
git config user.name Nishanth-nishu
git config user.email nishanth0962333@gmail.com
git config credential.https://github.com.helper '!f() { echo username=Nishanth-nishu; echo password=$(cat /scratch/nishanth.r/.github_token); }; f'
git remote get-url origin >/dev/null 2>&1 || git remote add origin $NEW

echo "=== STAGE 2: commit + push"
git add -A
git diff --cached --quiet || git commit -q -m "Phase 0: fresh start from TorDiff-NAX (prior Torsional Diffusion QM9 work)"
git log --format='%h %an <%ae> %s'
git push -u origin main 2>&1 | tail -3

echo "=== STAGE 3: baseline inventory (read-only)"
B=/scratch/nishanth.r/tordiff
ls $B; du -sh $B/data/* 2>/dev/null | head; ls $B/workdir 2>/dev/null | head -40
ls -la $B/workdir/qm9_default 2>/dev/null || true
cd $B/torsional-diffusion 2>/dev/null && git log --oneline -3 && ls
$B/venv/bin/python --version 2>&1 || true
$B/venv/bin/pip freeze 2>/dev/null | grep -iE 'torch|rdkit|e3nn|numpy|geometric|scatter|cluster|spyrmsd' || true
ls /scratch/nishanth.r/
df -h /scratch | tail -1
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
