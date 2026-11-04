#!/bin/bash
#SBATCH --job-name=jmcg_t1_signcheck
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --cpus-per-task=32
#SBATCH --mem=40G
#SBATCH --time=01:00:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# Review follow-up (addendum A): W-dihedral drift of the raw xTB / MMFF output BEFORE the exact restore, E1 seed 0.
# A sign-convention mismatch between xtb's $constrain dihedral and RDKit would show as large raw drift (≈ 2|τ|).
# Diagnostic only: writes to a separate directory, nothing in the reported comparison is replaced.
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage; TD=/scratch/nishanth.r/tordiff; R=$TD/torsional-diffusion
echo "job $SLURM_JOB_ID on $(hostname) start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
source /etc/profile.d/modules.sh 2>/dev/null || source /usr/share/Modules/init/bash; module load u22/xtb/6.7.1
set -euo pipefail
rsync -a --checksum $STAGE/flexitors/ $P/flexitors/; rsync -a --checksum $STAGE/experiments/exp1/ $P/experiments/exp1/
source $TD/venv/bin/activate; cd $P
export FLEXITORS_TD_ROOT=$R PYTHONPATH=$P CUDA_VISIBLE_DEVICES="" OMP_NUM_THREADS=1 TMPDIR=/scratch/nishanth.r/tmp_xtb_$SLURM_JOB_ID
mkdir -p $TMPDIR; trap 'rm -rf $TMPDIR' EXIT
OMP_NUM_THREADS=4 python -m pytest -q flexitors/tests 2>&1 | tail -1
IN=$P/results/e1_confs/n200; OUT=$P/results/exp1/T1_signcheck; mkdir -p $OUT
D=$STAGE/research/experiments/EXP1_V2_train/e1_signcheck; mkdir -p $D
for m in xtb mmff; do
  python experiments/exp1/t1_e1_readout.py --confs $IN/td0.pkl --method $m --workers 32 --out $OUT/td0_$m.pkl | tail -25
done
python - <<PY
import json, pickle, numpy as np
new = pickle.load(open('$OUT/td0_xtb.pkl', 'rb')); old = pickle.load(open('$P/results/exp1/T1_e1/td0_xtb.pkl', 'rb'))
d = [float(np.abs(a.GetConformer().GetPositions() - b.GetConformer().GetPositions()).max())
     for s in new for a, b in zip(new[s], old[s])]
print('rerun vs reported td0_xtb: max |coordinate difference| over all conformers', max(d))
PY
cp $OUT/*.log.json $D/
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
