# Shared runtime of every Stage A flexitors job on gnode118 (plan_P3 §3, round5/DECISION.md graft 2): sourced after
# slurm/common.sh. Provides: fx_env (venv, paths, commit), fx_preflight (fail fast, exit 2), with_oom_retry (CUDA OOM
# retry inside the job, same command and seeds, exit code preserved otherwise), fx_strict_check (expected outputs and a
# SUMMARY line, exit 3 = FX_FAIL), fx_ledger (one line per GPU task into the GPU-h ledger read by status.py).
# Memory rule: <= 5000 MB per CPU on plafnet2 (round2/FIXES.md X2). Never `set -e` here: callers enable it.

fx_env() {
    activate_env
    export CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-}
    export PYTHONPATH=$PROJECT
    export FLEXITORS_RES=$RES
    export FLEXITORS_COMMIT
    FLEXITORS_COMMIT=$(cat "$PROJECT/FLEXITORS_COMMIT")
    export PYTORCH_CUDA_ALLOC_CONF=max_split_size_mb:128   # plan_P3 §3.4 (fragmentation under co-located jobs)
    export FX_LEDGER=$RES/round5/gpu_ledger.tsv
    export MODEL_QM9=$WORK/pretrained/qm9_default
    export EZ_TABLE=$PROJECT/round5/data/ez_targets.csv
    export VAL_DIR=$RES/round5/valset
    cd "$REPO"
}

fx_preflight() {  # fx_preflight NEED_GPU [input paths...]
    local need_gpu=$1; shift
    python - <<'EOF' || { echo "FX_PREFLIGHT versions"; exit 2; }
import rdkit, torch, e3nn
assert rdkit.__version__ == '2022.09.5', rdkit.__version__   # EXP-3 §9.4: generation under the cluster RDKit
assert torch.__version__.startswith('1.13.1'), torch.__version__
assert e3nn.__version__ == '0.5.1', e3nn.__version__
from flexitors.td import changed_baseline_files
assert changed_baseline_files() == [], changed_baseline_files()
EOF
    for t in .p.npy .score.npy; do [[ -s $REPO/$t ]] || { echo "FX_PREFLIGHT torus table $t"; exit 2; }; done
    for p in "$@"; do [[ -e $p ]] || { echo "FX_PREFLIGHT missing input $p"; exit 2; }; done
    local free_gb; free_gb=$(df -BG --output=avail "$RES" | tail -1 | tr -dc 0-9)
    (( free_gb > 50 )) || { echo "FX_PREFLIGHT scratch free ${free_gb} GB"; exit 2; }
    if [[ $need_gpu == 1 ]]; then cuda_check || { echo "FX_PREFLIGHT cuda"; exit 2; }; fi
    echo "FX_PREFLIGHT ok"
}

with_oom_retry() {  # with_oom_retry LOG MAX -- cmd ...   (plan_P3 §3.4)
    local log=$1 max=$2; shift 3
    local n=0 rc
    while true; do
        rc=0; "$@" > "$log.try" 2>&1 || rc=$?   # never let set -e kill the retry loop
        cat "$log.try" >> "$log"
        (( rc == 0 )) && return 0
        if grep -qE 'CUDA out of memory|CUBLAS_STATUS_ALLOC_FAILED|CUDA error: out of memory' "$log.try" \
                && (( n < max )); then
            n=$((n + 1))
            echo "FX_OOM_RETRY $n"
            command -v nvidia-smi >/dev/null && nvidia-smi --query-compute-apps=pid,used_memory --format=csv
            sleep $(( 300 * 2 ** (n - 1) ))
            continue   # same command: generation resumes from confs.partial.pkl with the same streams
        fi
        return $rc
    done
}

fx_strict_check() {  # fx_strict_check OUT_DIR LOG file...
    local out=$1 log=$2; shift 2
    for f in "$@"; do [[ -s $out/$f ]] || { echo "FX_FAIL $out missing $f"; exit 3; }; done
    grep -q '^SUMMARY' "$log" || { echo "FX_FAIL $out no SUMMARY line"; exit 3; }
}

fx_ledger() {  # fx_ledger TAG START_EPOCH N_GPU
    local tag=$1 t0=$2 ngpu=$3 t1; t1=$(date +%s)
    mkdir -p "$(dirname "$FX_LEDGER")"
    printf '%s\t%s\t%s\t%s\t%.4f\n' "${SLURM_JOB_ID:-local}" "$tag" "$t0" "$t1" \
        "$(echo "($t1 - $t0) * $ngpu / 3600" | bc -l)" >> "$FX_LEDGER"
}
