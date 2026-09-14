#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/000.config.sh"

mkdir -p "${STATE_DIR}"

printf "task_id\tcutoff\tlabel\tcase_name\tcalculation_directory\tinput_file\n" \
    > "${TASK_LIST}"

TASK_ID=0

for CUTOFF in "${VOL_CUTOFFS[@]}"; do

    LABEL=$(cutoff_label "${CUTOFF}")

    CASE_NAME="VOL_Island_${LABEL}"

    CASE_DIR="${RUN_DIR}/${CASE_NAME}"

    INPUT_FILE="${CASE_DIR}/${CASE_NAME}.in"

    if [[ ! -d "${CASE_DIR}" ]]; then
        echo "ERROR: Calculation directory not found:"
        echo "${CASE_DIR}"
        echo
        echo "Run 001.setup_VOL_runs.sh first."
        exit 1
    fi

    if [[ ! -s "${INPUT_FILE}" ]]; then
        echo "ERROR: VOL input not found:"
        echo "${INPUT_FILE}"
        exit 1
    fi

    TASK_ID=$((TASK_ID + 1))

    printf "%d\t%s\t%s\t%s\t%s\t%s\n" \
        "${TASK_ID}" \
        "${CUTOFF}" \
        "${LABEL}" \
        "${CASE_NAME}" \
        "${CASE_DIR}" \
        "${INPUT_FILE}" \
        >> "${TASK_LIST}"

done

echo
echo "VOL task list created:"
echo "${TASK_LIST}"
echo
echo "Total VOL calculations: ${TASK_ID}"
