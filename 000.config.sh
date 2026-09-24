#!/bin/bash

WORK_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

VOL_CUTOFFS=(
    0.0
    0.1
    0.2
    0.3
    0.4
    0.5
    0.6
    0.7
    0.8
    0.9
    1.0
)

DOCK_ROOT="/gpfs/projects/rizzo/iamanor/DOCK6_Development/Similarity_Island/dock6_dev"

DOCK_BIN="${DOCK_ROOT}/bin/dock6"
DOCK_PARAMS="${DOCK_ROOT}/parameters"

VDW_DEFN_FILE="${DOCK_PARAMS}/vdw_AMBER_parm99.defn"
FLEX_DEFN_FILE="${DOCK_PARAMS}/flex.defn"

LIGAND_ATOM_FILE="/gpfs/projects/rizzo/iamanor/DOCK6_Development/Similarity_Island/TesTING_SMI_2/000_VS_For_DOCKING_SCORE/Dock_Scored_Molecules_for_SMI_Clustering_scored.mol2"

GRID_PREFIX="/gpfs/projects/AMS536/2026/group3/tutorial_isaac/Individual_Project/003_gridbox/grid"

TASK_LIST="${WORK_ROOT}/VOL_tasks.tsv"
STATUS_FILE="${WORK_ROOT}/VOL_status.tsv"

SLURM_PARTITION="rn-long-40core"
SLURM_TIME="2-00:00:00"

MAX_CONCURRENT_TASKS=11

cutoff_label()
{
    local cutoff="$1"
    echo "${cutoff//./p}"
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    echo
    echo "Volume Overlap Similarity Island workflow"
    echo
    echo "Working directory:"
    echo "${WORK_ROOT}"
    echo
    echo "DOCK binary:"
    echo "${DOCK_BIN}"
    echo
    echo "Ligand file:"
    echo "${LIGAND_ATOM_FILE}"
    echo
    echo "Grid prefix:"
    echo "${GRID_PREFIX}"
    echo
    echo "SLURM partition:"
    echo "${SLURM_PARTITION}"
    echo
    echo "Maximum concurrent calculations:"
    echo "${MAX_CONCURRENT_TASKS}"
    echo
    echo "Volume Overlap cutoffs:"

    for CUTOFF in "${VOL_CUTOFFS[@]}"; do
        echo "  ${CUTOFF}"
    done

    echo
    echo "Calculations: ${#VOL_CUTOFFS[@]}"
fi
