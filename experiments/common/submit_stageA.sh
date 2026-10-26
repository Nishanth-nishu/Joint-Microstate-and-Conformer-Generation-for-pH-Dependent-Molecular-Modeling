#!/bin/bash
# Stage A submission (plan_P3 §9 "Submission", §3.6 ledger stop). Runs ON gnode118 (it reads /scratch), e.g.
#   ssh ada 'srun -n1 -c1 --mem=2G -t 00:10:00 -A plafnet2 -p plafnet2 -w gnode118 \
#            env DRY_RUN=1 bash /scratch/nishanth.r/tordiff/experiments/common/submit_stageA.sh plan'
# Modes:
#   prep            CPU: Exp 3 affected-subset key files (flexitors.data.strata --keys_dir), QM9 test (+ DRUGS if staged)
#   smoke           every new arm on 20 molecules (smoke_tasks.tsv), one GPU job (fx_smoke)
#   w1 | w2 | w3    a wave of stageA_tasks.tsv: w1 validation + TD test pools; w2 PG test pools (setting from the
#                   validation decision) + Exp 3 QM9 arms; w3 GEOM-DRUGS (needs the staged download)
#   plan            DRY_RUN of all waves with GPU-h estimates, plus the Exp 1 lines not built in milestone 2
# Rows are packed into one GPU job per (wave, seed), executed in TSV order (bases before arms). Before any sbatch the
# ledger check runs: used + planned > 25 GPU-h -> FX_BUDGET_STOP, nothing is submitted (DECISION item 10).
set -euo pipefail
PROJECT=/scratch/nishanth.r/tordiff
E=$PROJECT/experiments/common
RES=${RES:-$PROJECT/results}
DRY_RUN=${DRY_RUN:-0}
MODE=${1:?mode}
EXP1_GPU_H=${EXP1_GPU_H:-11.2}   # EXP-1 smoke 0.18 + M1 3 x 30 epochs 10.5 + Gate-0 eval 0.5 (plan_P2 §7); milestone 3

sub() {  # sub NAME GPU EST ROWS TASKS
    local name=$1 gpu=$2 est=$3 rows=$4 tasks=$5 gres=()
    (( gpu )) && gres=(--gres=gpu:1)
    local cmd=(sbatch --parsable --job-name="$name" "${gres[@]}" --export=ALL,TASKS="$tasks",ROWS="$rows" "$E/infer.sbatch")
    if [[ $DRY_RUN == 1 ]]; then echo "DRY_RUN est=${est}GPU-h ${cmd[*]}"; else "${cmd[@]}"; fi
}

wave_jobs() {  # wave_jobs TSV WAVE JOBNAME -> prints "seed<TAB>rows<TAB>est" per job
    awk -F'\t' -v w="$2" 'NR > 1 && $11 == w { s = $4; rows[s] = rows[s] (rows[s] ? "+" : "") (NR - 1); est[s] += $10 }
        END { for (s in rows) printf "%s\t%s\t%.3f\n", s, rows[s], est[s] }' "$1" | sort
}

substitute_pg() {  # copy of the task table with PG:CHOSEN replaced by the validation decision
    local dec=$RES/round5/exp5/analysis/exp5_validation_decision.json out=$RES/round5/tasks_resolved.tsv choice
    if [[ -s $dec ]]; then
        choice=$(python -c "import json,sys; print(json.load(open(sys.argv[1]))['pg']['choice'])" "$dec")
    elif [[ $DRY_RUN == 1 ]]; then choice=noninv_recall; echo "DRY_RUN PG:CHOSEN -> $choice (placeholder: no decision yet)" >&2
    else echo "FX_FAIL no validation decision $dec: run experiments/exp5/analysis.py validation first" >&2; exit 3; fi
    sed "s/PG:CHOSEN/PG:$choice/" "$E/stageA_tasks.tsv" > "$out"
    echo "$out"
}

run_wave() {  # run_wave TSV WAVE JOBNAME
    local total; total=$(awk -F'\t' -v w="$2" 'NR > 1 && $11 == w { t += $10 } END { printf "%.3f", t }' "$1")
    if [[ $DRY_RUN != 1 ]]; then
        python "$E/status.py" --res "$RES" --logs "$PROJECT/logs" --planned "$total" \
            || { echo "FX_BUDGET_STOP wave $2 ($total GPU-h planned): nothing submitted"; exit 4; }
    fi
    echo "wave $2: planned $total GPU-h"
    while IFS=$'\t' read -r seed rows est; do
        sub "$3" 1 "$est" "$rows" "$1"
    done < <(wave_jobs "$1" "$2")
}

case $MODE in
    prep)
        source "$PROJECT/slurm/common.sh"; activate_env; export PYTHONPATH=$PROJECT
        python -m flexitors.data.strata --csv "$QM9_TEST_CSV" --out "$RES/round5/keys/qm9test_strata.csv" \
            --keys_dir "$RES/round5/keys"
        if [[ -s $DATA/DRUGS/test_smiles_corrected.csv ]]; then
            python -m flexitors.data.strata --csv "$DATA/DRUGS/test_smiles_corrected.csv" \
                --out "$RES/round5/keys/drugstest_strata.csv" --keys_dir "$RES/round5/keys" --keys_prefix drugs_
        fi ;;
    smoke) run_wave "$E/smoke_tasks.tsv" smoke fx_smoke ;;
    w1) run_wave "$E/stageA_tasks.tsv" w1 fx_inf ;;
    w2|w3) run_wave "$(substitute_pg)" "$MODE" fx_inf ;;
    plan)
        DRY_RUN=1
        tot=0
        for w in w1 w2 w3; do
            t=$(awk -F'\t' -v w="$w" 'NR > 1 && $11 == w { t += $10 } END { printf "%.3f", t }' "$E/stageA_tasks.tsv")
            tot=$(echo "$tot + $t" | bc -l)
            run_wave "$(substitute_pg)" "$w" fx_inf
        done
        echo "Exp 1 (milestone 3, not built yet): smoke + M1 + Gate-0 eval ~ ${EXP1_GPU_H} GPU-h"
        printf 'STAGE A ESTIMATE: Exp 3 + Exp 5 %.2f GPU-h + Exp 1 %.2f = %.2f GPU-h (cap 25)\n' \
            "$tot" "$EXP1_GPU_H" "$(echo "$tot + $EXP1_GPU_H" | bc -l)" ;;
    *) echo "unknown mode $MODE"; exit 2 ;;
esac
