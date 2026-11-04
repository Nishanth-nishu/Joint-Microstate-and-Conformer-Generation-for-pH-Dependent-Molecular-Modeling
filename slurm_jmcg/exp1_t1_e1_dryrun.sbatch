#!/bin/bash
#SBATCH --job-name=jmcg_t1_dry
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=02:00:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# Pipeline dry run of t1_e1_readout.py on 3 E1 molecules with the smoke-3 checkpoint (47102). Checks only that the
# script runs, torsions are preserved and outputs have the evaluate_confs.py format. No metric is computed.
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage; TD=/scratch/nishanth.r/tordiff
set -euo pipefail
rsync -a --checksum $STAGE/flexitors/ $P/flexitors/; rsync -a --checksum $STAGE/experiments/exp1/ $P/experiments/exp1/
source $TD/venv/bin/activate; cd $P
export FLEXITORS_TD_ROOT=$TD/torsional-diffusion PYTHONPATH=$P CUDA_VISIBLE_DEVICES="" OMP_NUM_THREADS=1
C=$P/results/e1_confs/n200/td0.pkl; O=$P/results/exp1/t1_dry; mkdir -p $O; R=$P/workdir/exp1/v2_resmoke_zero_s0
python experiments/exp1/t1_e1_readout.py --confs $C --method l --run_dir $R --project on --limit_mols 3 --workers 4 --out $O/l_on.pkl
python experiments/exp1/t1_e1_readout.py --confs $C --method l --run_dir $R --project off --limit_mols 3 --workers 4 --out $O/l_off.pkl
python experiments/exp1/t1_e1_readout.py --confs $C --method mmff --limit_mols 3 --workers 4 --out $O/mmff.pkl
python - <<PY
import pickle
for k in ('l_on', 'l_off', 'mmff'):
    d = pickle.load(open('$O/%s.pkl' % k, 'rb')); s = pickle.load(open('$C', 'rb'))
    print(k, {m: (len(v), len(s[m]), v[0].GetNumAtoms()) for m, v in d.items()})
PY
grep -h '"error"' $O/*.log.json | head -5 || true
