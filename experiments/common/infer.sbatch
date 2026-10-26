#!/bin/bash
#SBATCH --job-name=fx_inf
#SBATCH --partition=plafnet2
#SBATCH --account=plafnet2
#SBATCH --nodelist=gnode118
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4000M
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/nishanth.r/tordiff/logs/%x_%j.log
#
# TSV-driven Stage A inference (plan_P3 §9 "TSV-driven inference"): TASKS = path of a tasks TSV (header row), ROWS =
# row numbers joined by "+" (1-based, after the header) run in order inside this one GPU allocation (packing).
# Columns: TAG DATASET ARM SEED N KEYS BASE POST LIMIT GPUH_EST WAVE (written by make_tasks.py)
#   TAG      output dir under $RES/round5/            DATASET qm9test | qm9val | drugstest
#   ARM      O K1 K1_wrong K1_canonical K2_inf F_filter PG:<setting> | F_flip (CPU, on BASE) | O_prefix (cut 2K of BASE)
#   N        2K | 4K | pool          KEYS  key file under $RES/round5/keys/ or "-"      BASE  O / pool tag or "-"
#   POST     comma list of eval, stereo, features, validity, identity   LIMIT  molecules (0 = all; smoke uses 20)
# Every row: preflight -> generate (OOM retry, resumable) -> [splice] -> post-processing -> strict check -> GPU ledger.
PROJECT=/scratch/nishanth.r/tordiff
# shellcheck disable=SC1091
source "$PROJECT/slurm/common.sh"
# shellcheck disable=SC1091
source "$PROJECT/experiments/common/fx_common.sh"
diagnostics
set -euo pipefail
fx_env
NW=${SLURM_CPUS_PER_TASK:-8}
TASKS=${TASKS:?set TASKS}
ROWS=${ROWS:?set ROWS}
NGPU=$([[ -n "${CUDA_VISIBLE_DEVICES:-}" ]] && echo 1 || echo 0)
fx_preflight "$NGPU" "$TASKS" "$MODEL_QM9/best_model.pt" "$EZ_TABLE"

dataset_files() {  # sets CSV MOLS MODEL for a dataset name
    case $1 in
        qm9test) CSV=$QM9_TEST_CSV; MOLS=$QM9_TEST_MOLS; MODEL=$MODEL_QM9 ;;
        qm9val)  CSV=$VAL_DIR/val_smiles.csv; MOLS=$VAL_DIR/val_mols.pkl; MODEL=$MODEL_QM9 ;;
        drugstest) CSV=$DATA/DRUGS/test_smiles_corrected.csv; MOLS=$DATA/DRUGS/test_mols.pkl
                   MODEL=$WORK/pretrained/drugs_default ;;
        *) echo "FX_FAIL unknown dataset $1"; exit 3 ;;
    esac
}

run_row() {
    local TAG DATASET ARM SEED N KEYS BASE POST LIMIT GPUH WAVE
    IFS=$'\t' read -r TAG DATASET ARM SEED N KEYS BASE POST LIMIT GPUH WAVE < <(sed -n "$(( $1 + 1 ))p" "$TASKS")
    local out=$RES/round5/$TAG log; log=$out/run.log
    mkdir -p "$out"
    dataset_files "$DATASET"
    local t0; t0=$(date +%s)
    echo "ROW $1 $TAG $DATASET $ARM seed=$SEED n=$N keys=$KEYS base=$BASE post=$POST limit=$LIMIT est=$GPUH wave=$WAVE"
    local keyarg=() confs=$out/confs.pkl
    [[ $KEYS != - ]] && keyarg=(--keys "$RES/round5/keys/$KEYS")
    case $ARM in
        PG:CHOSEN) echo "FX_FAIL PG:CHOSEN must be substituted by submit_stageA.sh"; exit 3 ;;
        F_flip)
            python -m flexitors.sampling.stereo_fix --confs "$RES/round5/$BASE/confs.pkl" --ez_table "$EZ_TABLE" \
                --out_dir "$out" >> "$log" 2>&1 ;;
        O_prefix)
            python -m flexitors.sampling.splice prefix --pool "$RES/round5/$BASE/confs.pkl" --csv "$CSV" \
                --out "$confs" >> "$log" 2>&1
            cp "$RES/round5/$BASE/seeds.pkl" "$out/seeds.pkl"
            echo '{"nfe": 0, "source": "prefix of '"$BASE"'"}' > "$out/budget.json" ;;
        *)
            local armarg=(--arm "$ARM")
            [[ $ARM == PG:* ]] && armarg=(--arm PG --pg_setting "${ARM#PG:}")
            with_oom_retry "$log" 3 -- python -m flexitors.sampling.generate --model_dir "$MODEL" --csv "$CSV" \
                --out_dir "$out" "${armarg[@]}" --seed "$SEED" --n "$N" --ez_table "$EZ_TABLE" --batch_size 128 \
                --limit "$LIMIT" "${keyarg[@]}"
            fx_strict_check "$out" "$log" confs.pkl seeds.pkl budget.json provenance.json
            if [[ $BASE != - && $KEYS != - && $POST != *identity* ]]; then  # EXP-3 §9.1: spliced full set (the record), arm scored on its subset
                python -m flexitors.sampling.splice splice --o_confs "$RES/round5/$BASE/confs.pkl" --arm_confs "$confs" \
                    --keys "$RES/round5/keys/$KEYS" --out "$out/confs_spliced.pkl" >> "$log" 2>&1
            fi ;;
    esac
    local lim=(); (( LIMIT > 0 )) && [[ $KEYS == - ]] && lim=(--limit_mols "$LIMIT")  # same first rows as generate
    for p in ${POST//,/ }; do
        case $p in
            eval) python evaluate_confs.py --confs "$confs" --test_csv "$CSV" --true_mols "$MOLS" --dataset qm9 \
                      --n_workers "$NW" --out_results "$out/eval.pkl" "${lim[@]}" > "$out/eval.log" 2>&1 ;;
            stereo) python -m flexitors.eval.stereo_eval --confs "$confs" --csv "$CSV" --true_mols "$MOLS" \
                        --ez_table "$EZ_TABLE" --out "$out/stereo_eval.csv" >> "$log" 2>&1 ;;
            features) python -m flexitors.sampling.select --confs "$confs" --out "$out/features.pkl" >> "$log" 2>&1 ;;
            validity) python -m flexitors.eval.validity screen --confs "$confs" --seeds "$out/seeds.pkl" \
                          --ez_table "$EZ_TABLE" --out "$out/validity.pkl" >> "$log" 2>&1 ;;
            identity) python -m flexitors.sampling.splice splice --o_confs "$RES/round5/$BASE/confs.pkl" \
                          --arm_confs "$confs" --keys "$RES/round5/keys/$KEYS" --out "$out/confs_spliced.pkl" \
                          --identity_check "$confs" | tee -a "$log" ;;
            *) echo "FX_FAIL unknown post step $p"; exit 3 ;;
        esac
    done
    grep -h '^SUMMARY' "$log" | tail -3
    fx_ledger "$TAG" "$t0" "$NGPU"
    echo "ROW_DONE $1 $TAG"
}

for r in ${ROWS//[,+]/ }; do run_row "$r"; done
echo "done ROWS=$ROWS"
