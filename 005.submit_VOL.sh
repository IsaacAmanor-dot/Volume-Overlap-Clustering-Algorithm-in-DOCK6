#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/000.config.sh"

if [[ ! -s "${TASK_LIST}" ]]; then
    echo "ERROR: Missing task list."
    echo "Run bash 002.make_task_list.sh first."
    exit 1
fi

N_TASKS=$(awk 'END {print NR - 1}' "${TASK_LIST}")

if [[ "${N_TASKS}" -lt 1 ]]; then
    echo "ERROR: No VOL tasks found."
    exit 1
fi

CPUS="${MAX_CONCURRENT_TASKS}"

if [[ "${CPUS}" -gt "${N_TASKS}" ]]; then
    CPUS="${N_TASKS}"
fi

echo
echo "VOL calculations: ${N_TASKS}"
echo "Node count: 1"
echo "Concurrent calculations: ${CPUS}"
echo "Partition: ${SLURM_PARTITION}"
echo

sbatch \
    --partition="${SLURM_PARTITION}" \
    --time="${SLURM_TIME}" \
    --nodes=1 \
    --ntasks=1 \
    --cpus-per-task="${CPUS}" \
    --job-name="VOL_Island" \
    --output="${WORK_ROOT}/VOL_slurm_%j.out" \
    --export="ALL,WORKFLOW_DIR=${WORK_ROOT}" \
    "${WORK_ROOT}/004.run_VOL_chunks.slurm"
