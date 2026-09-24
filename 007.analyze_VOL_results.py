#!/usr/bin/env python3

import csv
import math
import re
import statistics
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


SCRIPT_DIR = Path(__file__).resolve().parent

ANALYSIS_DIR = SCRIPT_DIR / "zzz.analysis_results"
ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)

SUMMARY_CSV = ANALYSIS_DIR / "VOL_analysis_summary.csv"
ISLAND_CSV = ANALYSIS_DIR / "VOL_island_statistics.csv"
MEMBER_CSV = ANALYSIS_DIR / "VOL_member_statistics.csv"

ISLAND_COUNT_PLOT = ANALYSIS_DIR / "VOL_island_count_vs_cutoff.png"
ISLAND_SIZE_PLOT = ANALYSIS_DIR / "VOL_island_size_vs_cutoff.png"
SINGLETON_PLOT = ANALYSIS_DIR / "VOL_singleton_fraction_vs_cutoff.png"
SIZE_DISTRIBUTION_PLOT = ANALYSIS_DIR / "VOL_island_size_distribution.png"
SIMILARITY_PLOT = ANALYSIS_DIR / "VOL_head_member_similarity_vs_cutoff.png"
CALCULATIONS_PLOT = ANALYSIS_DIR / "VOL_calculations_vs_cutoff.png"
REDUCTION_PLOT = ANALYSIS_DIR / "VOL_calculation_reduction_vs_cutoff.png"
RUNTIME_PLOT = ANALYSIS_DIR / "VOL_runtime_vs_cutoff.png"


def mean_or_nan(values):
    return statistics.mean(values) if values else math.nan


def median_or_nan(values):
    return statistics.median(values) if values else math.nan


def stdev_or_nan(values):
    if not values:
        return math.nan
    if len(values) == 1:
        return 0.0
    return statistics.pstdev(values)


def min_or_nan(values):
    return min(values) if values else math.nan


def max_or_nan(values):
    return max(values) if values else math.nan


def cutoff_label(cutoff):
    if float(cutoff).is_integer():
        return str(int(cutoff))
    return f"{cutoff:g}"


def find_output_files():
    files = []

    for path in sorted(SCRIPT_DIR.glob("VOL_Island_*.out")):
        name = path.name.lower()

        if "summary" in name:
            continue

        if path.name.startswith("VOL_slurm_"):
            continue

        files.append(path)

    return files


def parse_output(output_file):
    run = {
        "source_file": str(output_file),
        "cutoff": math.nan,
        "molecules_to_cluster": 0,
        "molecules_clustered": 0,
        "total_islands": 0,
        "total_volume_overlap_calculations": 0,
        "avg_volume_overlap_calculations_per_island": math.nan,
        "elapsed_seconds": math.nan,
        "molecules_per_second": math.nan,
        "virtual_memory_kb": 0,
        "physical_memory_kb": 0,
        "clustering_complete": False,
        "islands": [],
    }

    current_island = None
    in_calculation_table = False
    island_calculations = {}

    cutoff_re = re.compile(
        r"^Similarity Cutoff:\s*([-+]?\d+(?:\.\d+)?)"
    )

    molecules_to_cluster_re = re.compile(
        r"^Molecules to Cluster:\s*(\d+)"
    )

    island_re = re.compile(
        r"^(?:Volume_Overlap|Volume Overlap|VOL)\s+ISLAND\s+(\d+)",
        re.IGNORECASE,
    )

    island_head_re = re.compile(
        r"^Island Head:\s+(\S+)"
    )

    island_size_re = re.compile(
        r"^Island Size:\s+(\d+)"
    )

    member_re = re.compile(
        r"^\s*(\d+)\s+(\S+)\s+"
        r"([-+]?\d+(?:\.\d+)?)\s+"
        r"(HEAD|MEMBER)\s*$"
    )

    molecules_clustered_re = re.compile(
        r"^Molecules Clustered:\s*(\d+)"
    )

    total_islands_re = re.compile(
        r"^Total Islands:\s*(\d+)"
    )

    total_calculations_re = re.compile(
        r"^Total Volume[_ ]Overlap Calculations:\s*(\d+)",
        re.IGNORECASE,
    )

    avg_calculations_re = re.compile(
        r"^Average Volume[_ ]Overlap Calculations per Island:\s*"
        r"([-+]?\d+(?:\.\d+)?)",
        re.IGNORECASE,
    )

    calculation_table_re = re.compile(
        r"^\s*(\d+)\s+(\d+)\s+(\d+)\s*$"
    )

    elapsed_re = re.compile(
        r"^Total elapsed time:\s*"
        r"([-+]?\d+(?:\.\d+)?)\s+seconds"
    )

    rate_re = re.compile(
        r"^Number of molecules per second:\s*"
        r"([-+]?\d+(?:\.\d+)?)"
    )

    virtual_memory_re = re.compile(
        r"^Virtual memory used for this process:\s*"
        r"(\d+)\s+kilobytes"
    )

    physical_memory_re = re.compile(
        r"^Physical memory used for this process:\s*"
        r"(\d+)\s+kilobytes"
    )

    with output_file.open("r", errors="replace") as handle:
        for raw_line in handle:
            line = raw_line.strip()

            match = cutoff_re.match(line)
            if match:
                run["cutoff"] = float(match.group(1))
                continue

            match = molecules_to_cluster_re.match(line)
            if match:
                run["molecules_to_cluster"] = int(match.group(1))
                continue

            match = island_re.match(line)
            if match:
                current_island = {
                    "island_id": int(match.group(1)),
                    "head": "",
                    "size": 0,
                    "volume_overlap_calculations": 0,
                    "members": [],
                }

                run["islands"].append(current_island)
                in_calculation_table = False
                continue

            if current_island is not None:
                match = island_head_re.match(line)

                if match:
                    current_island["head"] = match.group(1)
                    continue

                match = island_size_re.match(line)

                if match:
                    current_island["size"] = int(match.group(1))
                    continue

                match = member_re.match(line)

                if match:
                    rank = int(match.group(1))
                    molecule = match.group(2)
                    volume_overlap = float(match.group(3))
                    role = match.group(4)

                    if role == "MEMBER":
                        current_island["members"].append(
                            {
                                "rank": rank,
                                "molecule": molecule,
                                "volume_overlap": volume_overlap,
                            }
                        )

                    continue

            if line.startswith("Similarity Island clustering complete"):
                run["clustering_complete"] = True
                current_island = None
                continue

            match = molecules_clustered_re.match(line)

            if match:
                run["molecules_clustered"] = int(match.group(1))
                continue

            match = total_islands_re.match(line)

            if match:
                run["total_islands"] = int(match.group(1))
                continue

            match = total_calculations_re.match(line)

            if match:
                run["total_volume_overlap_calculations"] = int(
                    match.group(1)
                )
                continue

            match = avg_calculations_re.match(line)

            if match:
                run[
                    "avg_volume_overlap_calculations_per_island"
                ] = float(match.group(1))
                continue

            if (
                line
                == "Similarity Island Calculations by Volume_Overlap:"
                or line
                == "Similarity Island Calculations by Volume Overlap:"
            ):
                in_calculation_table = True
                continue

            if in_calculation_table:
                if line.startswith("Island"):
                    continue

                match = calculation_table_re.match(line)

                if match:
                    island_id = int(match.group(1))
                    calculations = int(match.group(2))
                    molecules = int(match.group(3))

                    island_calculations[island_id] = {
                        "calculations": calculations,
                        "molecules": molecules,
                    }

                    continue

                if line.startswith("Rescoring"):
                    in_calculation_table = False

            match = elapsed_re.match(line)

            if match:
                run["elapsed_seconds"] = float(match.group(1))
                continue

            match = rate_re.match(line)

            if match:
                run["molecules_per_second"] = float(match.group(1))
                continue

            match = virtual_memory_re.match(line)

            if match:
                run["virtual_memory_kb"] = int(match.group(1))
                continue

            match = physical_memory_re.match(line)

            if match:
                run["physical_memory_kb"] = int(match.group(1))
                continue

    for island in run["islands"]:
        info = island_calculations.get(island["island_id"])

        if info:
            island["volume_overlap_calculations"] = info["calculations"]

            if island["size"] == 0:
                island["size"] = info["molecules"]

    return run


def validate_run(run):
    warnings = []

    parsed_islands = len(run["islands"])

    sum_sizes = sum(
        island["size"]
        for island in run["islands"]
    )

    sum_calculations = sum(
        island["volume_overlap_calculations"]
        for island in run["islands"]
    )

    if not run["clustering_complete"]:
        warnings.append(
            "clustering completion marker was not found"
        )

    if (
        run["total_islands"] > 0
        and parsed_islands > 0
        and run["total_islands"] != parsed_islands
    ):
        warnings.append(
            f"reported islands={run['total_islands']}, "
            f"parsed islands={parsed_islands}"
        )

    if (
        run["molecules_clustered"] > 0
        and sum_sizes > 0
        and run["molecules_clustered"] != sum_sizes
    ):
        warnings.append(
            f"molecules clustered={run['molecules_clustered']}, "
            f"sum island sizes={sum_sizes}"
        )

    if (
        run["total_volume_overlap_calculations"] > 0
        and sum_calculations > 0
        and run["total_volume_overlap_calculations"]
        != sum_calculations
    ):
        warnings.append(
            "reported Volume Overlap calculations="
            f"{run['total_volume_overlap_calculations']}, "
            f"sum island calculations={sum_calculations}"
        )

    return warnings


def calculate_statistics(run):
    island_sizes = [
        island["size"]
        for island in run["islands"]
        if island["size"] > 0
    ]

    member_scores = []

    for island in run["islands"]:
        member_scores.extend(
            member["volume_overlap"]
            for member in island["members"]
        )

    run["largest_island"] = (
        max(island_sizes)
        if island_sizes
        else 0
    )

    run["smallest_island"] = (
        min(island_sizes)
        if island_sizes
        else 0
    )

    run["mean_island_size"] = mean_or_nan(island_sizes)
    run["median_island_size"] = median_or_nan(island_sizes)
    run["std_island_size"] = stdev_or_nan(island_sizes)

    run["singleton_islands"] = sum(
        size == 1
        for size in island_sizes
    )

    if run["total_islands"] > 0:
        run["singleton_fraction"] = (
            run["singleton_islands"]
            / run["total_islands"]
        )
    else:
        run["singleton_fraction"] = math.nan

    run["member_volume_overlap_count"] = len(member_scores)
    run["mean_member_volume_overlap"] = mean_or_nan(member_scores)
    run["median_member_volume_overlap"] = median_or_nan(member_scores)
    run["std_member_volume_overlap"] = stdev_or_nan(member_scores)
    run["min_member_volume_overlap"] = min_or_nan(member_scores)
    run["max_member_volume_overlap"] = max_or_nan(member_scores)

    for island in run["islands"]:
        scores = [
            member["volume_overlap"]
            for member in island["members"]
        ]

        island["member_count"] = len(scores)
        island["mean_volume_overlap"] = mean_or_nan(scores)
        island["median_volume_overlap"] = median_or_nan(scores)
        island["std_volume_overlap"] = stdev_or_nan(scores)
        island["min_volume_overlap"] = min_or_nan(scores)
        island["max_volume_overlap"] = max_or_nan(scores)

    n = run["molecules_clustered"]

    if n <= 0:
        n = run["molecules_to_cluster"]

    if n > 0:
        run["unique_pair_baseline"] = (
            n * (n - 1) // 2
        )

        run["implementation_max_calculations"] = (
            n * (n + 1) // 2
        )
    else:
        run["unique_pair_baseline"] = 0
        run["implementation_max_calculations"] = 0

    observed = run["total_volume_overlap_calculations"]
    implementation_max = run["implementation_max_calculations"]

    if implementation_max > 0:
        run["calculation_fraction_of_max"] = (
            observed / implementation_max
        )

        run["calculation_reduction_percent"] = (
            100.0
            * (
                1.0
                - run["calculation_fraction_of_max"]
            )
        )
    else:
        run["calculation_fraction_of_max"] = math.nan
        run["calculation_reduction_percent"] = math.nan


def write_summary_csv(runs):
    fields = [
        "cutoff",
        "molecules_to_cluster",
        "molecules_clustered",
        "total_islands",
        "largest_island",
        "smallest_island",
        "mean_island_size",
        "median_island_size",
        "std_island_size",
        "singleton_islands",
        "singleton_fraction",
        "total_volume_overlap_calculations",
        "avg_volume_overlap_calculations_per_island",
        "unique_pair_baseline",
        "implementation_max_calculations",
        "calculation_fraction_of_max",
        "calculation_reduction_percent",
        "member_volume_overlap_count",
        "mean_member_volume_overlap",
        "median_member_volume_overlap",
        "std_member_volume_overlap",
        "min_member_volume_overlap",
        "max_member_volume_overlap",
        "elapsed_seconds",
        "molecules_per_second",
        "virtual_memory_kb",
        "physical_memory_kb",
        "source_file",
    ]

    with SUMMARY_CSV.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
        )

        writer.writeheader()

        for run in runs:
            writer.writerow(
                {
                    field: run.get(field, "")
                    for field in fields
                }
            )


def write_island_csv(runs):
    fields = [
        "cutoff",
        "island_id",
        "island_head",
        "island_size",
        "volume_overlap_calculations",
        "member_count",
        "mean_volume_overlap",
        "median_volume_overlap",
        "std_volume_overlap",
        "min_volume_overlap",
        "max_volume_overlap",
        "source_file",
    ]

    with ISLAND_CSV.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
        )

        writer.writeheader()

        for run in runs:
            for island in run["islands"]:
                writer.writerow(
                    {
                        "cutoff": run["cutoff"],
                        "island_id": island["island_id"],
                        "island_head": island["head"],
                        "island_size": island["size"],
                        "volume_overlap_calculations":
                            island["volume_overlap_calculations"],
                        "member_count":
                            island["member_count"],
                        "mean_volume_overlap":
                            island["mean_volume_overlap"],
                        "median_volume_overlap":
                            island["median_volume_overlap"],
                        "std_volume_overlap":
                            island["std_volume_overlap"],
                        "min_volume_overlap":
                            island["min_volume_overlap"],
                        "max_volume_overlap":
                            island["max_volume_overlap"],
                        "source_file":
                            run["source_file"],
                    }
                )


def write_member_csv(runs):
    fields = [
        "cutoff",
        "island_id",
        "island_head",
        "island_size",
        "rank",
        "molecule",
        "volume_overlap",
    ]

    with MEMBER_CSV.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
        )

        writer.writeheader()

        for run in runs:
            for island in run["islands"]:
                for member in island["members"]:
                    writer.writerow(
                        {
                            "cutoff": run["cutoff"],
                            "island_id": island["island_id"],
                            "island_head": island["head"],
                            "island_size": island["size"],
                            "rank": member["rank"],
                            "molecule": member["molecule"],
                            "volume_overlap":
                                member["volume_overlap"],
                        }
                    )


def plot_island_count(runs):
    x = [run["cutoff"] for run in runs]
    y = [run["total_islands"] for run in runs]

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.plot(x, y, marker="o")

    ax.set_xlabel("Volume Overlap cutoff")
    ax.set_ylabel("Number of islands")

    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(ISLAND_COUNT_PLOT, dpi=300)
    plt.close(fig)


def plot_island_sizes(runs):
    x = [run["cutoff"] for run in runs]

    fig, ax = plt.subplots(figsize=(9, 6))

    ax.plot(
        x,
        [run["mean_island_size"] for run in runs],
        marker="o",
        label="Mean",
    )

    ax.plot(
        x,
        [run["largest_island"] for run in runs],
        marker="^",
        label="Largest",
    )

    ax.set_xlabel("Volume Overlap cutoff")
    ax.set_ylabel("Molecules per island")

    ax.set_yscale("log")
    ax.grid(alpha=0.3)
    ax.legend()

    fig.tight_layout()
    fig.savefig(ISLAND_SIZE_PLOT, dpi=300)
    plt.close(fig)


def plot_singleton_fraction(runs):
    x = [run["cutoff"] for run in runs]

    y = [
        100.0 * run["singleton_fraction"]
        for run in runs
    ]

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.plot(x, y, marker="o")

    ax.set_xlabel("Volume Overlap cutoff")
    ax.set_ylabel("Singleton islands (%)")

    ax.set_ylim(0, 105)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(SINGLETON_PLOT, dpi=300)
    plt.close(fig)


def plot_size_distribution(runs):
    fig, ax = plt.subplots(figsize=(10, 6))

    data = []
    positions = []

    for run in runs:
        sizes = [
            island["size"]
            for island in run["islands"]
            if island["size"] > 0
        ]

        if sizes:
            data.append(sizes)
            positions.append(run["cutoff"])

    if not data:
        plt.close(fig)
        return

    width = 0.035

    ax.boxplot(
        data,
        positions=positions,
        widths=width,
        showfliers=False,
        manage_ticks=False,
    )

    ax.set_xlabel("Volume Overlap cutoff")
    ax.set_ylabel("Molecules per island")

    ax.set_yscale("log")
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(SIZE_DISTRIBUTION_PLOT, dpi=300)
    plt.close(fig)


def plot_similarity(runs):
    valid_runs = [
        run
        for run in runs
        if run["member_volume_overlap_count"] > 0
    ]

    if not valid_runs:
        print(
            "Skipping head-member similarity plot: "
            "no member Volume Overlap values were parsed."
        )
        return

    fig, ax = plt.subplots(figsize=(9, 6))

    for run in valid_runs:
        scores = []

        for island in run["islands"]:
            scores.extend(
                member["volume_overlap"]
                for member in island["members"]
            )

        if scores:
            ax.scatter(
                [run["cutoff"]] * len(scores),
                scores,
                alpha=0.15,
                s=10,
            )

    x = [
        run["cutoff"]
        for run in valid_runs
    ]

    ax.plot(
        x,
        [
            run["mean_member_volume_overlap"]
            for run in valid_runs
        ],
        marker="o",
        label="Mean head-member overlap",
    )

    ax.plot(
        x,
        [
            run["median_member_volume_overlap"]
            for run in valid_runs
        ],
        marker="s",
        label="Median head-member overlap",
    )

    ax.set_xlabel("Volume Overlap cutoff")
    ax.set_ylabel("Head-member Volume Overlap")

    ax.set_ylim(0, 1.05)
    ax.grid(alpha=0.3)
    ax.legend()

    fig.tight_layout()
    fig.savefig(SIMILARITY_PLOT, dpi=300)
    plt.close(fig)


def plot_calculations(runs):
    x = [run["cutoff"] for run in runs]

    observed = [
        run["total_volume_overlap_calculations"]
        for run in runs
    ]

    unique_pairs = [
        run["unique_pair_baseline"]
        for run in runs
    ]

    implementation_max = [
        run["implementation_max_calculations"]
        for run in runs
    ]

    fig, ax = plt.subplots(figsize=(9, 6))

    ax.plot(
        x,
        observed,
        marker="o",
        label="Similarity Island calculations",
    )

    ax.plot(
        x,
        unique_pairs,
        linestyle="--",
        label="Unique-pair baseline: N(N-1)/2",
    )

    ax.plot(
        x,
        implementation_max,
        linestyle=":",
        label="Implementation maximum: N(N+1)/2",
    )

    ax.set_xlabel("Volume Overlap cutoff")
    ax.set_ylabel("Similarity calculations")

    ax.grid(alpha=0.3)
    ax.legend()

    fig.tight_layout()
    fig.savefig(CALCULATIONS_PLOT, dpi=300)
    plt.close(fig)

def plot_calculation_reduction(runs):
    x = [run["cutoff"] for run in runs]

    y = [
        run["calculation_reduction_percent"]
        for run in runs
    ]

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.axhline(
        0.0,
        linewidth=1,
        linestyle="--",
    )

    ax.plot(x, y, marker="o")

    ax.set_xlabel("Volume Overlap cutoff")
    ax.set_ylabel(
        "Reduction relative to N(N-1)/2 (%)"
    )

    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(REDUCTION_PLOT, dpi=300)
    plt.close(fig)


def plot_runtime(runs):
    valid = [
        (
            run["cutoff"],
            run["elapsed_seconds"],
        )
        for run in runs
        if not math.isnan(run["elapsed_seconds"])
    ]

    if not valid:
        return

    x, y = zip(*valid)

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.plot(x, y, marker="o")

    ax.set_xlabel("Volume Overlap cutoff")
    ax.set_ylabel("Elapsed time (seconds)")

    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(RUNTIME_PLOT, dpi=300)
    plt.close(fig)


def print_summary(runs):
    print()
    print("Volume Overlap Similarity Island Analysis")
    print()

    print(
        f"{'Cutoff':>7} "
        f"{'Molecules':>10} "
        f"{'Islands':>8} "
        f"{'Largest':>8} "
        f"{'Mean':>9} "
        f"{'Singleton%':>11} "
        f"{'VOLCalcs':>12} "
        f"{'MaxCalcs':>12} "
        f"{'Fraction':>10} "
        f"{'Reduction%':>11} "
        f"{'Time(s)':>10}"
    )

    for run in runs:
        print(
            f"{cutoff_label(run['cutoff']):>7} "
            f"{run['molecules_clustered']:>10d} "
            f"{run['total_islands']:>8d} "
            f"{run['largest_island']:>8d} "
            f"{run['mean_island_size']:>9.2f} "
            f"{100.0 * run['singleton_fraction']:>11.2f} "
            f"{run['total_volume_overlap_calculations']:>12d} "
            f"{run['implementation_max_calculations']:>12d} "
            f"{run['calculation_fraction_of_max']:>10.4f} "
            f"{run['calculation_reduction_percent']:>11.2f} "
            f"{run['elapsed_seconds']:>10.3f}"
        )

    print()

    print("Analysis files:")
    print(SUMMARY_CSV)
    print(ISLAND_CSV)
    print(MEMBER_CSV)

    print()

    print("Plots:")
    print(ISLAND_COUNT_PLOT)
    print(ISLAND_SIZE_PLOT)
    print(SINGLETON_PLOT)
    print(SIZE_DISTRIBUTION_PLOT)

    if SIMILARITY_PLOT.exists():
        print(SIMILARITY_PLOT)

    print(CALCULATIONS_PLOT)
    print(REDUCTION_PLOT)
    print(RUNTIME_PLOT)

    print()


def main():
    output_files = find_output_files()

    if not output_files:
        print("ERROR: No VOL output files found.")
        sys.exit(1)

    print()
    print(
        f"Found {len(output_files)} VOL output file(s)."
    )
    print()

    runs = []

    for output_file in output_files:
        print(f"Parsing: {output_file.name}")

        run = parse_output(output_file)

        if not run["clustering_complete"]:
            print(
                "WARNING: clustering completion marker "
                "not found. Skipping."
            )
            continue

        if math.isnan(run["cutoff"]):
            print(
                "WARNING: cutoff not found. Skipping."
            )
            continue

        if not run["islands"]:
            print(
                "WARNING: no islands parsed. Skipping."
            )
            continue

        warnings = validate_run(run)

        for warning in warnings:
            print(
                f"WARNING [{output_file.name}]: "
                f"{warning}"
            )

        calculate_statistics(run)
        runs.append(run)

    if not runs:
        print(
            "ERROR: No completed VOL runs could be analyzed."
        )
        sys.exit(1)

    runs.sort(
        key=lambda run: run["cutoff"]
    )

    write_summary_csv(runs)
    write_island_csv(runs)
    write_member_csv(runs)

    plot_island_count(runs)
    plot_island_sizes(runs)
    plot_singleton_fraction(runs)
    plot_size_distribution(runs)
    plot_similarity(runs)
    plot_calculations(runs)
    plot_calculation_reduction(runs)
    plot_runtime(runs)

    print_summary(runs)


if __name__ == "__main__":
    main()
