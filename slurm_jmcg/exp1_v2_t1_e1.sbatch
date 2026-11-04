#!/bin/bash
#SBATCH --job-name=jmcg_v2_t1_e1
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --cpus-per-task=16
#SBATCH --mem=40G
#SBATCH --time=12:00:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# EXP1-V2-T1 E1 evaluation, run ONCE after both training arms finish (design: research/experiments/EXP1_V2_train/DESIGN.md).
# Inputs: frozen-TD E1 samples td{0,1}.pkl (job 45308), best_model.pt of T1_proj_s0 / T1_noproj_s0 (validation-selected).
# Arms per seed: TD, TD+L(proj), TD+L(noproj), TD+MMFF94s(torsions restored), ETKDGv3. TD's evaluate_confs.py, paired bootstrap.
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage; TD=/scratch/nishanth.r/tordiff; R=$TD/torsional-diffusion
echo "job $SLURM_JOB_ID on $(hostname) start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
rsync -a --checksum $STAGE/flexitors/ $P/flexitors/; rsync -a --checksum $STAGE/experiments/exp1/ $P/experiments/exp1/
rsync -a $STAGE/jmcg/condgen/ $P/jmcg/condgen/
source $TD/venv/bin/activate; cd $P
export FLEXITORS_TD_ROOT=$R PYTHONPATH=$P CUDA_VISIBLE_DEVICES="" OMP_NUM_THREADS=1
IN=$P/results/e1_confs/n200; OUT=$P/results/exp1/T1_e1; W=$P/workdir/exp1; mkdir -p $OUT
D=$STAGE/research/experiments/EXP1_V2_train/e1; mkdir -p $D
sha256sum $W/T1_proj_s0/best_model.pt $W/T1_noproj_s0/best_model.pt | sed "s|$P/||"
for s in 0 1; do
  python experiments/exp1/t1_e1_readout.py --confs $IN/td$s.pkl --method l --run_dir $W/T1_proj_s0 --project on --workers 16 --out $OUT/td${s}_Lproj.pkl | tail -25
  python experiments/exp1/t1_e1_readout.py --confs $IN/td$s.pkl --method l --run_dir $W/T1_noproj_s0 --project off --workers 16 --out $OUT/td${s}_Lnoproj.pkl | tail -25
  python experiments/exp1/t1_e1_readout.py --confs $IN/td$s.pkl --method mmff --workers 16 --out $OUT/td${s}_mmff.pkl | tail -15
done
for k in td0_Lproj td0_Lnoproj td0_mmff td1_Lproj td1_Lnoproj td1_mmff; do
  (cd $R && python evaluate_confs.py --confs $OUT/$k.pkl --test_csv $IN/subset.csv --true_mols data/QM9/test_mols.pkl \
     --dataset qm9 --n_workers 16 --out_results $OUT/$k.eval.pkl > $OUT/$k.eval.log 2>&1) || echo "$k eval FAILED"
  grep -E "^SUMMARY" $OUT/$k.eval.log || true
done
for s in 0 1; do
  A="$OUT/td${s}_Lproj.eval.pkl $OUT/td${s}_Lnoproj.eval.pkl $OUT/td${s}_mmff.eval.pkl $IN/etkdg$s.eval.pkl"
  python -m jmcg.condgen.e1_confs compare --ref $IN/td$s.eval.pkl --arm $A --out $OUT/compare_seed$s.json > /dev/null
  python -m jmcg.condgen.e1_confs compare --ref $OUT/td${s}_Lnoproj.eval.pkl --arm $OUT/td${s}_Lproj.eval.pkl $OUT/td${s}_mmff.eval.pkl \
    $IN/td$s.eval.pkl $IN/etkdg$s.eval.pkl --out $OUT/ablation_seed$s.json > /dev/null
done
python - <<PY
import json
for s in (0, 1):
    for f, lab in (('compare', 'minus TD'), ('ablation', 'minus noproj')):
        r = json.load(open('$OUT/%s_seed%d.json' % (f, s)))
        print('seed', s, f, 'n_common', r['n_common'], 'failures', {k.split('/')[-1]: v for k, v in r['failures'].items()})
        if f == 'compare':
            for a, v in r['arms'].items(): print('  ', a.split('/')[-1], {k: round(x, 3) for k, x in v.items()})
        for a, v in r['paired_arm_minus_ref'].items():
            print('  ', lab, a.split('/')[-1], {k: (round(x['mean'], 4), [round(c, 4) for c in x['ci95']]) for k, x in v.items()
                  if k in ('MAT_R', 'COV_R@0.25', 'COV_R@0.5', 'MAT_P', 'COV_P@0.25', 'COV_P@0.5')})
PY
rsync -a --include '*.json' --include '*.log' --exclude '*' $OUT/ $D/
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
