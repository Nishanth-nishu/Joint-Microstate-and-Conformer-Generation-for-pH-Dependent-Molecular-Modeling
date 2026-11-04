#!/bin/bash
#SBATCH --job-name=jmcg_v2_t1_xtb
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --cpus-per-task=32
#SBATCH --mem=60G
#SBATCH --time=16:00:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# EXP1-V2-T1 addendum A (research/experiments/EXP1_V2_train/DESIGN.md): constrained GFN2-xTB control and the fixed
# L read-out (rigid molecules) on the frozen-TD E1 samples, then TD evaluate_confs and paired bootstraps.
# CPU only: plafnet2 GPU-minute quota exhausted, and xtb has no GPU path. Checkpoints unchanged (T1 best_model.pt).
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage; TD=/scratch/nishanth.r/tordiff; R=$TD/torsional-diffusion
echo "job $SLURM_JOB_ID on $(hostname) start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/Modules/init/bash; module load u22/xtb/6.7.1
set -euo pipefail
rsync -a --checksum $STAGE/flexitors/ $P/flexitors/; rsync -a --checksum $STAGE/experiments/exp1/ $P/experiments/exp1/
rsync -a $STAGE/jmcg/condgen/ $P/jmcg/condgen/
source $TD/venv/bin/activate; cd $P
export FLEXITORS_TD_ROOT=$R PYTHONPATH=$P CUDA_VISIBLE_DEVICES="" OMP_NUM_THREADS=1 TMPDIR=/scratch/nishanth.r/tmp_xtb_$SLURM_JOB_ID
mkdir -p $TMPDIR; trap 'rm -rf $TMPDIR' EXIT
xtb --version 2>&1 | grep -i "xtb version"
OMP_NUM_THREADS=4 python -m pytest -q flexitors/tests 2>&1 | tail -1
IN=$P/results/e1_confs/n200; OUT=$P/results/exp1/T1_e1; W=$P/workdir/exp1; NW=32
D=$STAGE/research/experiments/EXP1_V2_train/e1; mkdir -p $D
sha256sum $W/T1_proj_s0/best_model.pt $W/T1_noproj_s0/best_model.pt | sed "s|$P/||"
# gate: 3-molecule xtb dry run must run without errors before the full run
python experiments/exp1/t1_e1_readout.py --confs $IN/td0.pkl --method xtb --limit_mols 3 --workers $NW --out $OUT/dry_xtb.pkl | tail -20
python -c "import json,sys; s=json.load(open('$OUT/dry_xtb.pkl.log.json'))['summary']; sys.exit(1 if s['errors'] else 0)" \
  || { echo "xtb dry run had errors, stopping"; exit 1; }
for s in 0 1; do
  python experiments/exp1/t1_e1_readout.py --confs $IN/td$s.pkl --method xtb --workers $NW --out $OUT/td${s}_xtb.pkl | tail -20
  python experiments/exp1/t1_e1_readout.py --confs $IN/td$s.pkl --method l --run_dir $W/T1_proj_s0 --project on --workers $NW --out $OUT/td${s}_Lprojfix.pkl | tail -25
  python experiments/exp1/t1_e1_readout.py --confs $IN/td$s.pkl --method l --run_dir $W/T1_noproj_s0 --project off --workers $NW --out $OUT/td${s}_Lnoprojfix.pkl | tail -25
done
for k in td0_xtb td0_Lprojfix td0_Lnoprojfix td1_xtb td1_Lprojfix td1_Lnoprojfix; do
  (cd $R && python evaluate_confs.py --confs $OUT/$k.pkl --test_csv $IN/subset.csv --true_mols data/QM9/test_mols.pkl \
     --dataset qm9 --n_workers 16 --out_results $OUT/$k.eval.pkl > $OUT/$k.eval.log 2>&1) || echo "$k eval FAILED"
done
for s in 0 1; do
  python -m jmcg.condgen.e1_confs compare --ref $IN/td$s.eval.pkl --arm $OUT/td${s}_Lprojfix.eval.pkl \
    $OUT/td${s}_Lnoprojfix.eval.pkl $OUT/td${s}_xtb.eval.pkl $OUT/td${s}_mmff.eval.pkl $IN/etkdg$s.eval.pkl \
    --out $OUT/xtb_vs_td_seed$s.json > /dev/null
  python -m jmcg.condgen.e1_confs compare --ref $OUT/td${s}_xtb.eval.pkl --arm $OUT/td${s}_Lprojfix.eval.pkl \
    $OUT/td${s}_Lnoprojfix.eval.pkl $OUT/td${s}_Lproj.eval.pkl $IN/td$s.eval.pkl --out $OUT/vs_xtb_seed$s.json > /dev/null
  python -m jmcg.condgen.e1_confs compare --ref $OUT/td${s}_Lnoprojfix.eval.pkl --arm $OUT/td${s}_Lprojfix.eval.pkl \
    --out $OUT/ablation_fix_seed$s.json > /dev/null
done
python - <<PY
import json
for s in (0, 1):
    for f in ('xtb_vs_td', 'vs_xtb', 'ablation_fix'):
        r = json.load(open('$OUT/%s_seed%d.json' % (f, s)))
        print('seed', s, f, 'n_common', r['n_common'], 'failures', {k.split('/')[-1]: v for k, v in r['failures'].items()})
        if f == 'xtb_vs_td':
            for a, v in r['arms'].items(): print('  ', a.split('/')[-1], {k: round(x, 3) for k, x in v.items()})
        for a, v in r['paired_arm_minus_ref'].items():
            print('   minus ref', a.split('/')[-1], {k: (round(x['mean'], 4), [round(c, 4) for c in x['ci95']]) for k, x in v.items()
                  if k in ('MAT_R', 'COV_R@0.25', 'COV_R@0.5', 'MAT_P', 'COV_P@0.25', 'COV_P@0.5')})
PY
rsync -a --include '*.json' --include '*.log' --exclude '*' $OUT/ $D/
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
