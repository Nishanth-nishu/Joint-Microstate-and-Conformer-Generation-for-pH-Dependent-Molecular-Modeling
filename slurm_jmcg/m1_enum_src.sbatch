#!/bin/bash
#SBATCH --job-name=jmcg_m1enum
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# Read-only: Uni-pKa enumerator source and template file, copied to staging for review (Apache-2.0, ext/unipka_orig).
P=/scratch/nishanth.r/jmcg; D=/home2/nishanth.r/jmcg_stage/ext/unipka_orig_enumerator; mkdir -p $D
cp $P/ext/unipka_orig/enumerator/main.py $P/ext/unipka_orig/enumerator/*.tsv $P/ext/unipka_orig/LICENSE $D/
cp $P/ext/unipka_orig/scripts/infer_mean_ensemble.py $D/ 2>/dev/null
wc -l $D/*; sed -n 1,400p $D/main.py
