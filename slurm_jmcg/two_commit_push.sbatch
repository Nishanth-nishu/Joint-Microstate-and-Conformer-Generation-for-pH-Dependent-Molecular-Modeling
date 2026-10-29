#!/bin/bash
#SBATCH --job-name=jmcg_push2
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# Two reviewable commits: (1) E/Z lock default (condgen + smoke 5), (2) M2-M5 (everything else). DRY=1 shows what each
# commit would contain and undoes the staging; DRY=0 commits and pushes. The decision log goes in twice: commit 1 gets
# it truncated before the M2-M5 entry.
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage; DRY=${DRY:-1}
echo "job $SLURM_JOB_ID DRY=$DRY start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
rsync -a --checksum --exclude '*.pdf' --exclude 'papers/*.txt' $STAGE/research/ $P/research/
rsync -a --checksum $STAGE/jmcg/ $P/jmcg/; rsync -a --checksum $STAGE/slurm_jmcg/ $P/slurm_jmcg/
rsync -a --checksum $STAGE/data_manifests/ $P/data_manifests/; cp $STAGE/README.md $P/README.md
cp $STAGE/*.sbatch $P/slurm_jmcg/; mkdir -p $P/logs; cp $STAGE/logs/*.log $P/logs/ 2>/dev/null || true
cd $P
echo "HEAD $(git rev-parse --short HEAD), branch $(git branch --show-current)"
grep -qx 'research/papers/' .gitignore || printf 'research/papers/\nresearch/scholar_runs/*/text/\n' >> .gitignore
DL=research/10_decision_log.md
n=$(grep -n '^## 2026-10-28 — M2–M5' $DL | cut -d: -f1); head -n $((n - 2)) $STAGE/$DL > $DL
L1="jmcg/condgen/run.py jmcg/condgen/states.py jmcg/condgen/td_backend.py jmcg/condgen/test_backend.py
    slurm_jmcg/cg_smoke5.sbatch research/experiments/CG_smoke/RESULTS.md research/experiments/CG_smoke/run5
    $DL"
git add -- $L1
echo "=== commit 1 (E/Z lock) staged:"; git diff --cached --stat | tail -25
echo "condgen diff vs run5/diff_vs_HEAD.patch (empty = identical code change):"
git diff --cached HEAD -- jmcg/condgen | diff - research/experiments/CG_smoke/run5/diff_vs_HEAD.patch | head -20 || true
if [[ $DRY == 0 ]]; then git commit -q -m "E/Z lock default for input-specified double bonds: run.py --ez_policy lock|free; smoke 5 regression tests; outputs identical to smoke 4"; fi
cp $STAGE/$DL $DL
echo "results/ files tracked at HEAD: $(git ls-files results | wc -l) (repo convention: results/ is committed)"
git add -A
echo "=== commit 2 (M2-M5) staged:"; git diff --cached --stat | tail -60
echo "largest staged files (KB):"; git diff --cached --name-only -z | xargs -0 -r du -k 2>/dev/null | sort -rn | head -8
big=$(git diff --cached --name-only -z | xargs -0 -r du -k 2>/dev/null | awk '$1 > 50000' | wc -l)
if [[ $DRY == 0 && $big -eq 0 ]]; then
  git commit -q -m "M2-M5: canonical data v1; sequential baseline and E1/E6; E2 joint prototype (pre-registered fit and robustness); joint sampler demo"
  git log --format='%h %an <%ae> %s' -3; git log -2 --format='%B' | grep -ci 'co-authored' || echo "no co-author trailers"
  git push origin main 2>&1 | tail -2
else
  echo "not committing commit 2 (DRY=$DRY, files > 50 MB: $big)"; git reset -q
  [[ $DRY == 1 ]] && git checkout -q -- $DL 2>/dev/null || true
fi
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
