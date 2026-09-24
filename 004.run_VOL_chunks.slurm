#!/bin/bash

set -uo pipefail

WORKFLOW_DIR="${WORKFLOW_DIR:-${SLURM_SUBMIT_DIR:-$(pwd)}}"
source "${WORKFLOW_DIR}/000.config.sh"

TOTAL_TASKS=$(awk 'END {print NR - 1}' "${TASK_LIST}")

if [[ "${TOTAL_TASKS}" -lt 1 ]]; then
    echo "ERROR: No VOL tasks found."
    exit 1
fi

echo
echo "SLURM job: ${SLURM_JOB_ID:-unknown}"
echo "Node: $(hostname)"
echo "VOL calculations: ${TOTAL_TASKS}"
echo "Maximum concurrent calculations: ${MAX_CONCURRENT_TASKS}"
echo

PIDS=()
TASK_IDS=()
FAILURES=0

for ((TASK_ID=1; TASK_ID<=TOTAL_TASKS; TASK_ID++)); do
    while [[ "$(jobs -rp | wc -l)" -ge "${MAX_CONCURRENT_TASKS}" ]]; do
        sleep 2
    done

    echo "Launching VOL task ${TASK_ID}"

    bash "${WORKFLOW_DIR}/003.run_one_task.slurm" "${TASK_ID}" &

    PIDS+=("$!")
    TASK_IDS+=("${TASK_ID}")
done

for I in "${!PIDS[@]}"; do
    if wait "${PIDS[$I]}"; then
        echo "Task ${TASK_IDS[$I]} complete."
    else
        echo "Task ${TASK_IDS[$I]} failed."
        FAILURES=$((FAILURES + 1))
    fi
done

echo
echo "Tasks launched: ${#PIDS[@]}"
echo "Failures: ${FAILURES}"

if [[ "${FAILURES}" -gt 0 ]]; then
    exit 1
fi

exit 0
