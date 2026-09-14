#!/usr/bin/env python3

import csv
import math
import re
import statistics
import sys
from pathlib import Path

try:
    import matplotlib.pyplot as plt
except ImportError:
    print("ERROR: matplotlib is required.")
    sys.exit(1)


SCRIPT_DIR = Path(__file__).resolve().parent

RUN_DIR = SCRIPT_DIR
ANALYSIS_DIR = SCRIPT_DIR

SUMMARY_CSV = SCRIPT_DIR / "VOL_analysis_summary.csv"
ISLAND_CSV = SCRIPT_DIR / "VOL_island_statistics.csv"
MEMBER_CSV = SCRIPT_DIR / "VOL_member_statistics.csv"

RUNTIME_PLOT = SCRIPT_DIR / "VOL_runtime_vs_cutoff.png"
ISLAND_COUNT_PLOT = SCRIPT_DIR / "VOL_island_count_vs_cutoff.png"
CLUSTER_SIZE_PLOT = SCRIPT_DIR / "VOL_cluster_size_vs_cutoff.png"
SIMILARITY_PLOT = SCRIPT_DIR / "VOL_similarity_vs_cutoff.png"
CALCULATIONS_PLOT = SCRIPT_DIR / "VOL_calculations_vs_cutoff.png"


def mean_value(values):
    return statistics.mean(values) if values else math.nan


def median_value(values):
    return statistics.median(values) if values else math.nan


def stdev_value(values):
    if len(values) > 1:
        return statistics.stdev(values)

    if len(values) == 1:
        return 0.0

    return math.nan


def min_value(values):
    return min(values) if values else math.nan


def max_value(values):
    return max(values) if values else math.nan


def discover_outputs():

    files = []

    for path in sorted(SCRIPT_DIR.glob("VOL_Island_*.out")):

        if "summary" in path.name.lower():
            continue

        if path.name.startswith("VOL_slurm_"):
            continue

        files.append(path)

    return files

def parse_output(path):

    run = {
        "source_file": str(path),
        "cutoff": math.nan,
        "molecules_clustered": 0,
        "total_islands": 0,
        "total_calculations": 0,
        "average_calculations_per_island": math.nan,
        "elapsed_seconds": math.nan,
        "molecules_per_second": math.nan,
        "virtual_memory_kb": 0,
        "physical_memory_kb": 0,
        "islands": [],
    }

    current = None

    cutoff_re = re.compile(
        r"^Similarity Cutoff:\s*([-+]?\d+(?:\.\d+)?)"
    )

    island_re = re.compile(
        r"^Volume_Overlap ISLAND\s+(\d+)"
    )

    head_re = re.compile(
        r"^Island Head:\s+(\S+)"
    )

    size_re = re.compile(
        r"^Island Size:\s+(\d+)"
    )

    member_re = re.compile(
        r"^\s*(\d+)\s+(\S+)\s+"
        r"([-+]?\d+(?:\.\d+)?)\s+"
        r"(HEAD|MEMBER)\s*$"
    )

    molecules_re = re.compile(
        r"^Molecules Clustered:\s*(\d+)"
    )

    total_islands_re = re.compile(
        r"^Total Islands:\s*(\d+)"
    )

    total_calculations_re = re.compile(
        r"^Total .*Calculations:\s*(\d+)"
    )

    average_calculations_re = re.compile(
        r"^Average .*Calculations per Island:\s*"
        r"([-+]?\d+(?:\.\d+)?)"
    )

    elapsed_re = re.compile(
        r"^Total elapsed time:\s*"
        r"([-+]?\d+(?:\.\d+)?)\s+seconds"
    )

    rate_re = re.compile(
        r"^Number of molecules per second:\s*"
        r"([-+]?\d+(?:\.\d+)?)"
    )

    virtual_re = re.compile(
        r"^Virtual memory used for this process:\s*(\d+)"
    )

    physical_re = re.compile(
        r"^Physical memory used for this process:\s*(\d+)"
    )

    with path.open("r", errors="replace") as handle:

        for raw in handle:

            line = raw.strip()

            match = cutoff_re.match(line)

            if match:
                run["cutoff"] = float(match.group(1))
                continue

            match = island_re.match(line)

            if match:

                current = {
                    "island_id": int(match.group(1)),
                    "head": "",
                    "size": 0,
                    "members": [],
                }

                run["islands"].append(current)
                continue

            if current is not None:

                match = head_re.match(line)

                if match:
                    current["head"] = match.group(1)
                    continue

                match = size_re.match(line)

                if match:
                    current["size"] = int(match.group(1))
                    continue

                match = member_re.match(line)

                if match:

                    rank = int(match.group(1))
                    molecule = match.group(2)
                    score = float(match.group(3))
                    role = match.group(4)

                    if role == "MEMBER":

                        current["members"].append(
                            {
                                "rank": rank,
                                "molecule": molecule,
                                "volume_overlap": score,
                            }
                        )

                    continue

            if line.startswith(
                "Similarity Island clustering complete"
            ):
                current = None
                continue

            match = molecules_re.match(line)

            if match:
                run["molecules_clustered"] = int(match.group(1))
                continue

            match = total_islands_re.match(line)

            if match:
                run["total_islands"] = int(match.group(1))
                continue

            match = total_calculations_re.match(line)

            if match:
                run["total_calculations"] = int(match.group(1))
                continue

            match = average_calculations_re.match(line)

            if match:
                run["average_calculations_per_island"] = float(
                    match.group(1)
                )
                continue

            match = elapsed_re.match(line)

            if match:
                run["elapsed_seconds"] = float(match.group(1))
                continue

            match = rate_re.match(line)

            if match:
                run["molecules_per_second"] = float(match.group(1))
                continue

            match = virtual_re.match(line)

            if match:
                run["virtual_memory_kb"] = int(match.group(1))
                continue

            match = physical_re.match(line)

            if match:
                run["physical_memory_kb"] = int(match.group(1))
                continue

    if run["total_islands"] == 0:
        run["total_islands"] = len(run["islands"])

    return run


def calculate_statistics(run):

    island_sizes = [
        island["size"]
        for island in run["islands"]
        if island["size"] > 0
    ]

    all_scores = []

    for island in run["islands"]:

        scores = [
            member["volume_overlap"]
            for member in island["members"]
        ]

        all_scores.extend(scores)

        island["member_count"] = len(scores)
        island["mean_overlap"] = mean_value(scores)
        island["median_overlap"] = median_value(scores)
        island["std_overlap"] = stdev_value(scores)
        island["min_overlap"] = min_value(scores)
        island["max_overlap"] = max_value(scores)

    run["largest_island"] = max(island_sizes) if island_sizes else 0
    run["smallest_island"] = min(island_sizes) if island_sizes else 0

    run["mean_island_size"] = mean_value(island_sizes)
    run["median_island_size"] = median_value(island_sizes)
    run["std_island_size"] = stdev_value(island_sizes)

    run["singleton_islands"] = sum(
        size == 1
        for size in island_sizes
    )

    run["member_overlap_count"] = len(all_scores)

    run["mean_member_overlap"] = mean_value(all_scores)
    run["median_member_overlap"] = median_value(all_scores)
    run["std_member_overlap"] = stdev_value(all_scores)
    run["min_member_overlap"] = min_value(all_scores)
    run["max_member_overlap"] = max_value(all_scores)


def write_summary(runs):

    fields = [
        "cutoff",
        "molecules_clustered",
        "total_islands",
        "largest_island",
        "smallest_island",
        "mean_island_size",
        "median_island_size",
        "std_island_size",
        "singleton_islands",
        "total_calculations",
        "average_calculations_per_island",
        "member_overlap_count",
        "mean_member_overlap",
        "median_member_overlap",
        "std_member_overlap",
        "min_member_overlap",
        "max_member_overlap",
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


def write_islands(runs):

    fields = [
        "cutoff",
        "island_id",
        "island_head",
        "island_size",
        "member_count",
        "mean_overlap",
        "median_overlap",
        "std_overlap",
        "min_overlap",
        "max_overlap",
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
                        "member_count": island["member_count"],
                        "mean_overlap": island["mean_overlap"],
                        "median_overlap": island["median_overlap"],
                        "std_overlap": island["std_overlap"],
                        "min_overlap": island["min_overlap"],
                        "max_overlap": island["max_overlap"],
                    }
                )


def write_members(runs):

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


def save_plot(path):

    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close()


def plot_runtime(runs):

    x = [
        run["cutoff"]
        for run in runs
        if not math.isnan(run["elapsed_seconds"])
    ]

    y = [
        run["elapsed_seconds"]
        for run in runs
        if not math.isnan(run["elapsed_seconds"])
    ]

    if not x:
        return

    plt.figure(figsize=(8, 6))
    plt.plot(x, y, marker="o")
    plt.xlabel("Volume Overlap cutoff")
    plt.ylabel("Elapsed time (seconds)")
    plt.title("Volume Overlap Runtime vs Cutoff")
    plt.grid(alpha=0.3)

    save_plot(RUNTIME_PLOT)


def plot_island_count(runs):

    x = [run["cutoff"] for run in runs]
    y = [run["total_islands"] for run in runs]

    plt.figure(figsize=(8, 6))
    plt.plot(x, y, marker="o")
    plt.xlabel("Volume Overlap cutoff")
    plt.ylabel("Number of islands")
    plt.title("Volume Overlap Island Count vs Cutoff")
    plt.grid(alpha=0.3)

    save_plot(ISLAND_COUNT_PLOT)


def plot_cluster_sizes(runs):

    plt.figure(figsize=(9, 6))

    for run in runs:

        sizes = [
            island["size"]
            for island in run["islands"]
            if island["size"] > 0
        ]

        if sizes:

            plt.scatter(
                [run["cutoff"]] * len(sizes),
                sizes,
                alpha=0.35,
                s=18,
            )

    x = [run["cutoff"] for run in runs]

    plt.plot(
        x,
        [run["mean_island_size"] for run in runs],
        marker="o",
        label="Mean island size",
    )

    plt.plot(
        x,
        [run["median_island_size"] for run in runs],
        marker="s",
        label="Median island size",
    )

    plt.plot(
        x,
        [run["largest_island"] for run in runs],
        marker="^",
        label="Largest island",
    )

    plt.xlabel("Volume Overlap cutoff")
    plt.ylabel("Molecules per island")
    plt.title("Volume Overlap Island Size vs Cutoff")
    plt.yscale("log")
    plt.grid(alpha=0.3)
    plt.legend()

    save_plot(CLUSTER_SIZE_PLOT)


def plot_similarity(runs):

    plt.figure(figsize=(9, 6))

    for run in runs:

        means = [
            island["mean_overlap"]
            for island in run["islands"]
            if not math.isnan(island["mean_overlap"])
        ]

        if means:

            plt.scatter(
                [run["cutoff"]] * len(means),
                means,
                alpha=0.35,
                s=18,
            )

    x = [run["cutoff"] for run in runs]

    plt.plot(
        x,
        [run["mean_member_overlap"] for run in runs],
        marker="o",
        label="Mean member overlap",
    )

    plt.plot(
        x,
        [run["median_member_overlap"] for run in runs],
        marker="s",
        label="Median member overlap",
    )

    plt.xlabel("Volume Overlap cutoff")
    plt.ylabel("Volume Overlap score")
    plt.title("Within-Island Volume Overlap vs Cutoff")

    plt.ylim(0.0, 1.0)

    plt.grid(alpha=0.3)
    plt.legend()

    save_plot(SIMILARITY_PLOT)


def plot_calculations(runs):

    valid = [
        run
        for run in runs
        if run["total_calculations"] > 0
    ]

    if not valid:
        return

    x = [run["cutoff"] for run in valid]
    y = [run["total_calculations"] for run in valid]

    plt.figure(figsize=(8, 6))
    plt.plot(x, y, marker="o")
    plt.xlabel("Volume Overlap cutoff")
    plt.ylabel("Total Volume Overlap calculations")
    plt.title("Volume Overlap Calculations vs Cutoff")
    plt.grid(alpha=0.3)

    save_plot(CALCULATIONS_PLOT)


def main():

    ANALYSIS_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    outputs = discover_outputs()

    if not outputs:
        print("ERROR: No Volume Overlap output files found.")
        sys.exit(1)

    runs = []

    print()
    print(
        f"Found {len(outputs)} Volume Overlap "
        f"output file(s)."
    )
    print()

    for output in outputs:

        print(f"Parsing: {output}")

        run = parse_output(output)

        if math.isnan(run["cutoff"]):
            print("WARNING: cutoff not found. Skipping.")
            continue

        if not run["islands"]:
            print("WARNING: no islands found. Skipping.")
            continue

        calculate_statistics(run)

        runs.append(run)

    if not runs:
        print("ERROR: No completed VOL runs could be analyzed.")
        sys.exit(1)

    runs.sort(
        key=lambda item: item["cutoff"]
    )

    write_summary(runs)
    write_islands(runs)
    write_members(runs)

    plot_runtime(runs)
    plot_island_count(runs)
    plot_cluster_sizes(runs)
    plot_similarity(runs)
    plot_calculations(runs)

    print()
    print("Volume Overlap analysis complete.")
    print()

    print("Cutoff   Islands   Largest   Mean size   Mean overlap")

    for run in runs:

        print(
            f"{run['cutoff']:>6.2f} "
            f"{run['total_islands']:>9d} "
            f"{run['largest_island']:>9d} "
            f"{run['mean_island_size']:>11.2f} "
            f"{run['mean_member_overlap']:>12.4f}"
        )

    print()
    print("Analysis directory:")
    print(ANALYSIS_DIR)
    print()


if __name__ == "__main__":
    main()
