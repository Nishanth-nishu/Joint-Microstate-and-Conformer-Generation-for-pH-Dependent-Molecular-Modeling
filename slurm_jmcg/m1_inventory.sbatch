#!/bin/bash
#SBATCH --job-name=jmcg_m1inv
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# Milestone 1 read-only inventory (2026-10-28): git state, external data, enumerators, model weights, envs, compute.
P=/scratch/nishanth.r/jmcg
echo "job $SLURM_JOB_ID start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
cd $P; echo "=== git"; git branch --show-current; git log --oneline | head -12; git status --short | head -30
echo "=== ext"; du -sh ext/* 2>/dev/null; ls ext/unipka_orig ext/unipka_orig/* 2>/dev/null | head -60
echo "=== enumerator candidates"; find ext -maxdepth 4 \( -iname '*enumer*' -o -iname '*smarts*' -o -iname '*template*' -o -iname '*dimorphite*' \) 2>/dev/null | head -30
echo "=== datasets"; for f in ext/unipka_orig/dataset/*; do echo "$f $(wc -l < $f) lines"; head -2 $f | cut -c1-300; done 2>/dev/null
ls ext/sampl6_pka ext/sampl8 ext/tautobase ext/iupac_pka 2>/dev/null | head -40
echo "=== results"; ls results results/* 2>/dev/null | head -60
echo "=== model weights"; find / -maxdepth 6 -path /proc -prune -o \( -name '*.pt' -path '*unipka*' \) -print 2>/dev/null | head; find $P -name '*.pt' -size +1M 2>/dev/null | head -20; ls ~/.cache 2>/dev/null
echo "=== venvs"; ls -d $P/venv* /scratch/nishanth.r/tordiff/venv*; $P/venv/bin/python -c "import rdkit, torch, sys; print(sys.version.split()[0], 'rdkit', rdkit.__version__, 'torch', torch.__version__)"; $P/venv/bin/pip list 2>/dev/null | grep -iE "dimorphite|openbabel|molvs|unipka|unimol|unicore|lightgbm|sklearn|scikit|xtb|crest|openmm|pandas|numpy" 
which obabel xtb crest 2>/dev/null
echo "=== GEOM"; ls /scratch/nishanth.r/tordiff/data /scratch/nishanth.r/tordiff/data/QM9 2>/dev/null | head; ls /scratch/nishanth.r/tordiff/data/DRUGS 2>/dev/null | head
echo "=== compute"; nproc; free -g | head -2; nvidia-smi -L 2>&1 | head; scontrol show node gnode118 | grep -E "CfgTRES|AllocTRES|Gres|State"
squeue -w gnode118 -o "%i %u %j %T %C %b %M" | head -20
df -h /scratch | tail -1
