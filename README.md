# DOCK6 Volume Overlap Similarity Island Workflow

This repository provides a reproducible workflow for benchmarking the Volume Overlap implementation of Similarity Island clustering in DOCK6.

A ranked or scored multi-MOL2 molecular library is clustered repeatedly using different Volume Overlap cutoff values. Each cutoff is treated as an independent DOCK6 calculation.

The workflow generates the DOCK6 input files, runs cutoff calculations concurrently using SLURM, checks calculation completion, and performs quantitative comparison of clustering behavior across the tested cutoffs.

## Volume Overlap Similarity Island Method

Volume Overlap clustering is enabled in DOCK6 using:

    cluster_by_similarity_island yes
    similarity_island_scoring_type volume_overlap

The clustering threshold is controlled by:

    similarity_island_volume_overlap_cutoff CUTOFF

Volume Overlap scores use the DOCK6 Volume Overlap scoring range from 0.0 to 1.0.

A value closer to 1.0 indicates greater volumetric overlap between the molecule and the current island head.

Each unassigned molecule selected as an island head is compared against the remaining unassigned molecules. Molecules satisfying the selected cutoff are assigned to that island. The next unassigned molecule becomes the next island head and clustering continues until all molecules have been assigned.

## Workflow Files

The workflow contains:

    000.config.sh
    001.setup_VOL_runs.sh
    002.make_task_list.sh
    003.run_one_task.slurm
    004.run_VOL_chunks.slurm
    005.submit_VOL.sh
    006.check_VOL_results.sh
    007.analyze_VOL_results.py
    README.md

## 000.config.sh

Contains all shared workflow parameters.

These include:

- Volume Overlap cutoffs
- DOCK6 Similarity Island executable
- input multi-MOL2 library
- grid prefix
- DOCK6 parameter files
- SLURM configuration
- workflow directories

The default cutoff set covers the complete Volume Overlap range:

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

Users may freely modify this list.

Display the current configuration with:

    bash 000.config.sh

## 001.setup_VOL_runs.sh

Creates one independent calculation directory and DOCK6 input file for each configured Volume Overlap cutoff.

Run:

    bash 001.setup_VOL_runs.sh

For example, cutoff 0.8 produces:

    runs/VOL/VOL_Island_0p8/

containing:

    VOL_Island_0p8.in

The corresponding output files use the same cutoff-specific identifier.

## 002.make_task_list.sh

Creates the task table used by the execution workflow.

Run:

    bash 002.make_task_list.sh

The task table is written to:

    state/VOL_tasks.tsv

Each row represents one independent cutoff calculation.

## 003.run_one_task.slurm

Runs one serial Volume Overlap Similarity Island calculation.

This script is normally called automatically by the node-level worker.

## 004.run_VOL_chunks.slurm

Runs multiple independent cutoff calculations concurrently on one compute node.

Each DOCK6 calculation remains serial.

No MPI is used.

## 005.submit_VOL.sh

Submits the complete Volume Overlap experiment to SLURM.

Run:

    bash 005.submit_VOL.sh

The number of simultaneous nodes is controlled by `MAX_NODES` in `000.config.sh`.

A temporary maximum may also be supplied:

    bash 005.submit_VOL.sh 2

## 006.check_VOL_results.sh

Checks all expected cutoff calculations.

Run:

    bash 006.check_VOL_results.sh

Calculations are classified as:

    SUCCESS
    FAILED
    INCOMPLETE
    MISSING

The complete status table is written to:

    state/VOL_status.tsv

## 007.analyze_VOL_results.py

Performs quantitative analysis across all completed Volume Overlap cutoff calculations.

Run:

    python 007.analyze_VOL_results.py

The analysis generates:

    analysis/VOL_analysis_summary.csv
    analysis/VOL_island_statistics.csv
    analysis/VOL_member_statistics.csv

It also produces independent plots comparing:

- calculation runtime
- total number of islands
- island size distributions
- within-island Volume Overlap scores
- total similarity calculations

The plots are:

    VOL_runtime_vs_cutoff.png
    VOL_island_count_vs_cutoff.png
    VOL_cluster_size_vs_cutoff.png
    VOL_similarity_vs_cutoff.png
    VOL_calculations_vs_cutoff.png

Island heads are excluded from member similarity statistics because the head has a self-overlap score of 1.0. The member statistics therefore represent actual head-to-member Volume Overlap relationships.

## Running the Workflow

First inspect or modify the experiment configuration:

    bash 000.config.sh

Generate all cutoff-specific DOCK6 inputs:

    bash 001.setup_VOL_runs.sh

Create the task table:

    bash 002.make_task_list.sh

Submit the calculations:

    bash 005.submit_VOL.sh

Check calculation completion:

    bash 006.check_VOL_results.sh

After the required calculations finish, run the analysis:

    python 007.analyze_VOL_results.py

## Output Organization

Generated calculations are stored under:

    runs/VOL/

For example:

    runs/VOL/VOL_Island_0p0/
    runs/VOL/VOL_Island_0p1/
    runs/VOL/VOL_Island_0p2/
    runs/VOL/VOL_Island_0p3/
    runs/VOL/VOL_Island_0p4/
    runs/VOL/VOL_Island_0p5/
    runs/VOL/VOL_Island_0p6/
    runs/VOL/VOL_Island_0p7/
    runs/VOL/VOL_Island_0p8/
    runs/VOL/VOL_Island_0p9/
    runs/VOL/VOL_Island_1p0/

This structure keeps every input, DOCK6 output, Similarity Island summary, and molecular output associated with its exact cutoff.

## Analysis Goal

The benchmark is designed to determine how the Volume Overlap cutoff affects:

- number of Similarity Islands
- island population
- singleton frequency
- largest island size
- mean and median island size
- within-island similarity
- computational cost
- runtime

Together, these measurements help identify cutoff values that provide useful separation of molecular populations while maintaining computational efficiency.
# Volume-Overlap-Clustering-Algorithm-in-DOCK6
# Volume-Overlap-Clustering-Algorithm-in-DOCK6
