#!/bin/bash
#SBATCH --job-name=jmcg_t2_manifest
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# EXP1-V2-T2 part D, before the final evaluation: (1) the E1 input manifest (sha256 of every saved input the final
# evaluation reads), with any file modified after the part-A audit (2026-10-31 02:20 UTC) flagged; (2) a pipeline dry
# run of t2_geom_validate.py --arm and t2_final_analysis.py with the T1 "fix" outputs standing in for the six T2
# read-outs and the partial T2 training logs (no new test-set number: every E1 value it touches is already reported).
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage; TD=/scratch/nishanth.r/tordiff; R=$TD/torsional-diffusion
echo "job $SLURM_JOB_ID start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
IN=$P/results/e1_confs/n200; T1=$P/results/exp1/T1_e1; X=$STAGE/research/experiments/EXP1_V2_T2
F="$IN/subset.csv $R/data/QM9/test_mols.pkl"
for s in 0 1; do
  F="$F $IN/td$s.pkl $IN/td$s.eval.pkl $IN/etkdg$s.eval.pkl"
  for k in xtb mmff Lprojfix Lnoprojfix; do F="$F $T1/td${s}_$k.pkl $T1/td${s}_$k.eval.pkl"; done
done
[ -e $X/e1_manifest.sha256 ] && { echo "manifest exists, not overwriting"; exit 1; }
sha256sum $F > $X/e1_manifest.sha256; cat $X/e1_manifest.sha256
echo "modified after the audit (should be none):"; find $F -newermt '2026-10-31 02:20:00 UTC' -print
ls -l --time-style=full-iso $F

source $TD/venv/bin/activate; cd $P
rsync -a --checksum $STAGE/experiments/exp1/ $P/experiments/exp1/; rsync -a $STAGE/jmcg/condgen/ $P/jmcg/condgen/
export PYTHONPATH=$P
T=$P/workdir/tmp_t2dry_$SLURM_JOB_ID; mkdir -p $T; trap 'rm -rf $T' EXIT
for s in 0 1; do
  for k in 0 1 2; do
    for m in proj noproj; do
      for e in pkl eval.pkl; do ln -s $T1/td${s}_L${m}fix.$e $T/td${s}_T2_${m}_s$k.$e; done
      ln -s $T1/td${s}_L${m}fix.pkl.log.json $T/td${s}_T2_${m}_s$k.pkl.log.json
    done
  done
  A=""; for r in proj_s0 noproj_s0 proj_s1 noproj_s1 proj_s2 noproj_s2; do A="$A --arm ${r}=$T/td${s}_T2_$r"; done
  python experiments/exp1/t2_geom_validate.py --in_dir $IN --t1_dir $T1 --true_mols $R/data/QM9/test_mols.pkl --seed $s \
    --workers 8 $A --out $T/geom_validate_seed$s.json | tail -5
done
python experiments/exp1/t2_final_analysis.py --n200 $IN --t1 $T1 --t2 $T --geom_dir $T --train_dir $P/workdir/exp1 \
  --out $T/final_analysis.json
echo "expected: proj_s* = Lprojfix (MAT_R 0.111 / 0.106), proj_s0-Lprojfix = 0, between-seed SD = 0, geometry as part C"
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
