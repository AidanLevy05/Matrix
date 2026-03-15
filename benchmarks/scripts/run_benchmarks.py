#!/usr/bin/env python3

import csv
import re
import statistics
import subprocess
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = ROOT / "benchmarks" / "data"
LOG_DIR = DATA_DIR / "logs"
FIG_DIR = ROOT / "benchmarks" / "figures"

SIZES = [250, 500, 1000, 1500]
PROCS = [1, 2, 4, 8, 16]
MAX_ATTEMPTS = 3

BENCHMARKS = {
    "multiply": {
        "binary": "multiplyTest",
        "parallel_pattern": re.compile(
            r"multiplyMatrix\(\) done in ([0-9]+(?:\.[0-9]+)?) seconds"
        ),
        "sequential_pattern": re.compile(
            r"Seq_multiplyMatrix\(\) done in ([0-9]+(?:\.[0-9]+)?) seconds"
        ),
    },
    "ref": {
        "binary": "refTest",
        "parallel_pattern": re.compile(r"ref\(\) done in ([0-9]+(?:\.[0-9]+)?) seconds"),
        "sequential_pattern": re.compile(
            r"Seq_ref\(\) done in ([0-9]+(?:\.[0-9]+)?) seconds"
        ),
    },
    "rref": {
        "binary": "rrefTest",
        "parallel_pattern": re.compile(
            r"rref\(\) done in ([0-9]+(?:\.[0-9]+)?) seconds"
        ),
        "sequential_pattern": re.compile(
            r"Seq_rref\(\) done in ([0-9]+(?:\.[0-9]+)?) seconds"
        ),
    },
    "lu": {
        "binary": "luTest",
        "parallel_pattern": re.compile(r"LU\(\) done in ([0-9]+(?:\.[0-9]+)?) seconds"),
        "sequential_pattern": re.compile(
            r"Sequential Seq_LU\(\) done in ([0-9]+(?:\.[0-9]+)?) seconds"
        ),
    },
}

def write_csv(path, rows, fieldnames):
    with path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_case(operation, size, procs, include_sequential):
    benchmark = BENCHMARKS[operation]
    binary = ROOT / benchmark["binary"]
    log_name = f"{operation}_n{size}_p{procs}_seq{int(include_sequential)}.log"
    log_path = LOG_DIR / log_name
    payload = f"{size}\n{1 if include_sequential else 0}\n"

    for attempt in range(1, MAX_ATTEMPTS + 1):
        started = time.perf_counter()
        completed = subprocess.run(
            ["mpirun", "--oversubscribe", "-np", str(procs), str(binary)],
            input=payload,
            cwd=ROOT,
            text=True,
            capture_output=True,
            check=False,
        )
        wall_seconds = time.perf_counter() - started

        output = completed.stdout + completed.stderr
        log_path.write_text(output, encoding="utf-8")

        parallel_match = benchmark["parallel_pattern"].search(output)
        sequential_match = benchmark["sequential_pattern"].search(output)

        if completed.returncode == 0 and parallel_match is not None:
            return {
                "operation": operation,
                "size": size,
                "procs": procs,
                "include_sequential": int(include_sequential),
                "attempt": attempt,
                "wall_seconds": wall_seconds,
                "parallel_seconds": float(parallel_match.group(1)),
                "sequential_seconds": (
                    float(sequential_match.group(1))
                    if sequential_match is not None
                    else None
                ),
                "log_path": log_path.relative_to(ROOT).as_posix(),
            }

        if attempt < MAX_ATTEMPTS:
            time.sleep(1.1)

    raise RuntimeError(
        f"Benchmark failed for {operation} size={size} procs={procs}; "
        f"see {log_path.relative_to(ROOT)}"
    )


def build_summary(raw_rows):
    grouped = {}
    for row in raw_rows:
        grouped.setdefault((row["operation"], row["size"]), []).append(row)

    summary_rows = []
    for key, rows in grouped.items():
        operation, size = key
        sequential_row = next(
            row for row in rows if row["mode"] == "sequential" and row["procs"] == 1
        )
        mpi_rows = [row for row in rows if row["mode"] == "mpi"]
        best_mpi = min(mpi_rows, key=lambda row: row["kernel_seconds"])
        summary_rows.append(
            {
                "operation": operation,
                "size": size,
                "sequential_seconds": sequential_row["kernel_seconds"],
                "best_mpi_procs": best_mpi["procs"],
                "best_mpi_seconds": best_mpi["kernel_seconds"],
                "best_speedup": sequential_row["kernel_seconds"]
                / best_mpi["kernel_seconds"],
            }
        )

    return sorted(summary_rows, key=lambda row: (row["operation"], row["size"]))


def plot_timings(raw_rows):
    plt.style.use("seaborn-v0_8-whitegrid")
    colors = {
        "sequential": "#222222",
        1: "#1b5e20",
        2: "#1565c0",
        4: "#ef6c00",
        8: "#6a1b9a",
        16: "#ad1457",
    }

    for operation in BENCHMARKS:
        fig, ax = plt.subplots(figsize=(8.5, 5.5))
        op_rows = [row for row in raw_rows if row["operation"] == operation]
        seq_rows = sorted(
            [row for row in op_rows if row["mode"] == "sequential"],
            key=lambda row: row["size"],
        )
        ax.plot(
            [row["size"] for row in seq_rows],
            [row["kernel_seconds"] for row in seq_rows],
            marker="o",
            linewidth=2.5,
            color=colors["sequential"],
            label="Sequential",
        )

        for procs in PROCS:
            mpi_rows = sorted(
                [
                    row
                    for row in op_rows
                    if row["mode"] == "mpi" and row["procs"] == procs
                ],
                key=lambda row: row["size"],
            )
            ax.plot(
                [row["size"] for row in mpi_rows],
                [row["kernel_seconds"] for row in mpi_rows],
                marker="o",
                linewidth=2,
                color=colors[procs],
                label=f"MPI np={procs}",
            )

        ax.set_title(f"{operation.upper()} Kernel Timing")
        ax.set_xlabel("Matrix size N (NxN)")
        ax.set_ylabel("Kernel time (seconds)")
        ax.legend(ncol=2, fontsize=9)
        fig.tight_layout()
        fig.savefig(FIG_DIR / f"{operation}_timing.png", dpi=220)
        plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True, sharey=True)
    for axis, operation in zip(axes.flat, BENCHMARKS):
        op_rows = [row for row in raw_rows if row["operation"] == operation]
        seq_by_size = {
            row["size"]: row["kernel_seconds"]
            for row in op_rows
            if row["mode"] == "sequential"
        }

        for procs in PROCS:
            mpi_rows = sorted(
                [
                    row
                    for row in op_rows
                    if row["mode"] == "mpi" and row["procs"] == procs
                ],
                key=lambda row: row["size"],
            )
            speedups = [
                seq_by_size[row["size"]] / row["kernel_seconds"] for row in mpi_rows
            ]
            axis.plot(
                [row["size"] for row in mpi_rows],
                speedups,
                marker="o",
                linewidth=2,
                color=colors[procs],
                label=f"np={procs}",
            )

        axis.set_title(operation.upper())
        axis.set_xlabel("Matrix size N")
        axis.set_ylabel("Speedup vs sequential")

    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=len(PROCS), frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(FIG_DIR / "speedup_summary.png", dpi=220)
    plt.close(fig)


def write_machine_info():
    memory_info = subprocess.run(
        ["free", "-h"],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=True,
    ).stdout.strip()
    cpu_count = subprocess.run(
        ["nproc"],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=True,
    ).stdout.strip()

    machine_info = [
        f"cpu_logical_cores={cpu_count}",
        "matrix_sizes=" + ",".join(str(size) for size in SIZES),
        "mpi_process_counts=" + ",".join(str(proc) for proc in PROCS),
        "timings=kernel times reported by the test executables",
        "",
        memory_info,
    ]
    (DATA_DIR / "machine_info.txt").write_text(
        "\n".join(machine_info) + "\n", encoding="ascii"
    )


def main():
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    subprocess.run(["make"], cwd=ROOT, check=True)

    captured_runs = []
    for operation in BENCHMARKS:
        for size in SIZES:
            print(f"[run] {operation} size={size} np=1 with sequential baseline")
            captured_runs.append(run_case(operation, size, 1, True))

            for procs in PROCS[1:]:
                print(f"[run] {operation} size={size} np={procs}")
                captured_runs.append(run_case(operation, size, procs, False))

    raw_rows = []
    for run in captured_runs:
        raw_rows.append(
            {
                "operation": run["operation"],
                "size": run["size"],
                "procs": run["procs"],
                "mode": "mpi",
                "kernel_seconds": run["parallel_seconds"],
                "wall_seconds": run["wall_seconds"],
                "attempt": run["attempt"],
                "log_path": run["log_path"],
            }
        )

        if run["sequential_seconds"] is not None:
            raw_rows.append(
                {
                    "operation": run["operation"],
                    "size": run["size"],
                    "procs": 1,
                    "mode": "sequential",
                    "kernel_seconds": run["sequential_seconds"],
                    "wall_seconds": "",
                    "attempt": run["attempt"],
                    "log_path": run["log_path"],
                }
            )

    raw_rows.sort(key=lambda row: (row["operation"], row["size"], row["mode"], row["procs"]))
    summary_rows = build_summary(raw_rows)

    write_csv(
        DATA_DIR / "timings.csv",
        raw_rows,
        [
            "operation",
            "size",
            "procs",
            "mode",
            "kernel_seconds",
            "wall_seconds",
            "attempt",
            "log_path",
        ],
    )
    write_csv(
        DATA_DIR / "summary.csv",
        summary_rows,
        [
            "operation",
            "size",
            "sequential_seconds",
            "best_mpi_procs",
            "best_mpi_seconds",
            "best_speedup",
        ],
    )
    plot_timings(raw_rows)
    write_machine_info()

    best_overall = {}
    for row in summary_rows:
        best_overall.setdefault(row["operation"], []).append(row["best_speedup"])

    print("\n[summary]")
    for operation in BENCHMARKS:
        operation_speedups = best_overall[operation]
        print(
            f"  {operation:8s} "
            f"median best speedup = {statistics.median(operation_speedups):.3f}x"
        )


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"benchmark failure: {exc}", file=sys.stderr)
        sys.exit(1)
