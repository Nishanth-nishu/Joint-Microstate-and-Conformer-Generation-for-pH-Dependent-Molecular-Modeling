#!/bin/bash
#SBATCH --job-name=jmcg_cgsmoke3
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/home2/nishanth.r/jmcg_stage/logs/%x_%j.log
# Conditional (explicit-state) conformer generation, smoke test 3 = smoke 2 unchanged except the QM9 checkpoint is PINNED to
# pretrained/qm9_default (sha256 6a9a6ee9...; smoke 1/2 auto-picked qm9_1order). Wrapper fidelity vs TD on this checkpoint: CG_fidelity. CPU only. Frozen TD checkout + venv used read-only.
# Separate from N1 / T pilot: own package jmcg/condgen, own outputs results/cond_gen/, research/experiments/CG_smoke/.
P=/scratch/nishanth.r/jmcg; STAGE=/home2/nishanth.r/jmcg_stage; TD=/scratch/nishanth.r/tordiff
echo "job $SLURM_JOB_ID on $(hostname) start $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
set -euo pipefail
rsync -a $STAGE/jmcg/condgen/ $P/jmcg/condgen/; rsync -a $STAGE/slurm_jmcg/ $P/slurm_jmcg/
source $TD/venv/bin/activate
cd $P
export PYTHONPATH=$P TD_REPO=$TD/torsional-diffusion OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK CUDA_VISIBLE_DEVICES=""
OUT=$P/results/cond_gen/smoke3; mkdir -p $OUT
echo "=== 0: env"
python -c 'import torch, rdkit, torch_geometric, e3nn; print("torch", torch.__version__, "rdkit", rdkit.__version__, "pyg", torch_geometric.__version__, "e3nn", e3nn.__version__)
from rdkit.Chem import rdDetermineBonds; print("rdDetermineBonds ok")'
(cd $TD/torsional-diffusion && git rev-parse HEAD && git status --short | head)
echo "=== 1: unit tests (states)"
python -m jmcg.condgen.test_states || echo "TEST_STATES FAILED (exit $?)"
echo "=== 2: checkpoint inventory + QM9 charge prevalence"
python -m jmcg.condgen.inventory --workdir $TD/workdir --qm9_raw $TD/data/QM9/qm9/ --qm9_split $TD/data/QM9/split.npy \
  --out $OUT/inventory.json
QM9_DIR=$TD/workdir/pretrained/qm9_default  # pinned (not auto-picked)
echo "$(sha256sum $QM9_DIR/best_model.pt)"; sha256sum $QM9_DIR/best_model.pt | grep -q "^6a9a6ee9b0e51bda" || { echo "CHECKPOINT SHA MISMATCH"; exit 1; }
DRUGS_DIR=$(python - <<EOF
import json; c=[x for x in json.load(open("$OUT/inventory.json"))["checkpoints"] if x["dataset"]=="drugs" and x["best_model_sha256"]]
print(c[0]["dir"] if c else "")
EOF
)
echo "qm9 checkpoint: ${QM9_DIR:-NONE}   drugs checkpoint: ${DRUGS_DIR:-NONE}"
test -n "$QM9_DIR"
echo "=== 2b: backend / SDF / geometry-validation tests"
python -m jmcg.condgen.test_backend --model_dir $QM9_DIR --states_csv jmcg/condgen/smoke_states.csv || echo "TEST_BACKEND FAILED (exit $?)"
echo "=== 3: smoke generation (released QM9 checkpoint), 8 confs, seed 0"
python -m jmcg.condgen.run --model_dir $QM9_DIR --out $OUT/qm9 --states_csv jmcg/condgen/smoke_states.csv \
  --n_confs 8 --seed 0 --lock_arm
if [[ -n "$DRUGS_DIR" ]]; then
  echo "=== 3b: same states with the drugs checkpoint"
  python -m jmcg.condgen.run --model_dir $DRUGS_DIR --out $OUT/drugs --states_csv jmcg/condgen/smoke_states.csv \
    --n_confs 8 --seed 0 --lock_arm || echo "drugs run failed (see above)"
fi
echo "=== 4: reproducibility (same seed twice) + seed sensitivity"
python - <<'EOF' > $OUT/repro_states.csv
import pandas as pd
d = pd.read_csv('jmcg/condgen/smoke_states.csv', comment='#')
print(d[d.id.isin(['glycine_zwitterion', 'acac_Z_enol', 'histamine_dication'])].to_csv(index=False), end='')
EOF
for tag in rep_a rep_b; do python -m jmcg.condgen.run --model_dir $QM9_DIR --out $OUT/$tag --states_csv $OUT/repro_states.csv --n_confs 8 --seed 0 > $OUT/$tag.log; done
python -m jmcg.condgen.run --model_dir $QM9_DIR --out $OUT/rep_seed1 --states_csv $OUT/repro_states.csv --n_confs 8 --seed 1 > $OUT/rep_seed1.log
python - <<EOF
import glob, os, numpy as np
from rdkit import Chem
def xyz(p): return [m.GetConformer().GetPositions() for m in Chem.SDMolSupplier(p, removeHs=False)]
for f in sorted(glob.glob("$OUT/rep_a/*.sdf")):
    b = f.replace("rep_a", "rep_b"); s1 = f.replace("rep_a", "rep_seed1")
    A, B, C = xyz(f), xyz(b), xyz(s1)
    dab = max(float(np.abs(x - y).max()) for x, y in zip(A, B))
    dac = max(float(np.abs(x - y).max()) for x, y in zip(A, C))
    print(f"{os.path.basename(f):40s} seed0 two processes max|dx|={dab:.1e} A  seed0 vs seed1 max|dx|={dac:.2f} A")
EOF
echo "=== 5: copy small outputs to staging"
D=$STAGE/research/experiments/CG_smoke/run3; mkdir -p $D
rsync -a --include '*/' --include '*.json' --include '*.csv' --include '*.sdf' --include '*.log' --exclude '*' $OUT/ $D/
du -sh $D
echo "end $(date -u '+%F %T UTC') / $(TZ=Asia/Kolkata date '+%F %T IST')"
