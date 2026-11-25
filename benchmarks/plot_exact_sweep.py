#!/usr/bin/env python3
"""Summarize exact TBC sweep and produce per-dataset plots.

Usage: python benchmarks/plot_exact_sweep.py results/exact_sweep/runs.csv

Only stdlib is required for summary CSV. Matplotlib is optional for PNG charts.
"""
from __future__ import annotations

import csv
import hashlib
import math
import re
import statistics
import sys
from collections import defaultdict
from pathlib import Path

LABELS = {
    "base": "C++ Base",
    "vertex": "C++ Vertex",
    "all": "C++ Vertex + Edges",
    "brunelli": "Julia / Brunelli",
}


def slugify(s: str) -> str:
    cleaned = re.sub(r"[^a-zA-Z0-9_.-]+", "_", s)[:65]
    return cleaned + "_" + hashlib.sha1(s.encode()).hexdigest()[:9]


def median_or_blank(numbers: list[float]):
    return statistics.median(numbers) if numbers else ""


def main() -> int:
    if len(sys.argv) != 2:
        print(__doc__)
        return 2
    csv_path = Path(sys.argv[1]).resolve()
    if not csv_path.is_file():
        raise FileNotFoundError(csv_path)

    grouped = defaultdict(list)
    with csv_path.open("r", newline="", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            if row.get("status") != "ok" or not row.get("compute_s"):
                continue
            key = (row["dataset"], row["algorithm"], int(row["threads"]))
            grouped[key].append(row)

    stats = {}
    for key, rows in grouped.items():
        times = [float(row["compute_s"]) for row in rows]
        rams = [float(row["peak_ram_mib"]) for row in rows if row.get("peak_ram_mib")]
        prep = [float(row["preprocess_s"]) for row in rows if row.get("preprocess_s")]
        stats[key] = {
            "n": len(times),
            "median_time_s": statistics.median(times),
            "mean_time_s": statistics.mean(times),
            "std_time_s": statistics.stdev(times) if len(times) > 1 else "",
            "median_peak_ram_mib": median_or_blank(rams),
            "median_preprocess_s": median_or_blank(prep),
        }

    output_dir = csv_path.parent
    summary_path = output_dir / "summary_median.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as f:
        fields = ["dataset", "algorithm", "threads", "successful_runs", "median_time_s",
                  "mean_time_s", "std_time_s", "median_peak_ram_mib", "median_preprocess_s",
                  "speedup_vs_1thread", "speed_ratio_to_brunelli_same_threads"]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for (dataset, algo, threads), result in sorted(stats.items()):
            baseline_1 = stats.get((dataset, algo, 1))
            brunelli = stats.get((dataset, "brunelli", threads))
            speedup = (baseline_1["median_time_s"] / result["median_time_s"]
                       if baseline_1 and result["median_time_s"] > 0 else "")
            ratio = (result["median_time_s"] / brunelli["median_time_s"]
                     if brunelli and brunelli["median_time_s"] > 0 else "")
            writer.writerow({
                "dataset": dataset, "algorithm": algo, "threads": threads,
                "successful_runs": result["n"],
                "median_time_s": result["median_time_s"],
                "mean_time_s": result["mean_time_s"],
                "std_time_s": result["std_time_s"],
                "median_peak_ram_mib": result["median_peak_ram_mib"],
                "median_preprocess_s": result["median_preprocess_s"],
                "speedup_vs_1thread": speedup,
                "speed_ratio_to_brunelli_same_threads": ratio,
            })

    print("Median summary:", summary_path)
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("Matplotlib not installed. Run 'python -m pip install matplotlib' for PNG charts.")
        return 0

    plots_dir = output_dir / "plots"
    plots_dir.mkdir(exist_ok=True)
    datasets = sorted({d for d, _a, _t in stats})
    plot_count = 0
    for dataset in datasets:
        slug = slugify(dataset)
        all_threads = sorted({p for (d, _a, p) in stats if d == dataset})
        # One separate figure for absolute compute time.
        fig, ax = plt.subplots(figsize=(9.0, 5.4))
        for algo, label in LABELS.items():
            points = [(p, stats[(dataset, algo, p)]["median_time_s"])
                      for p in all_threads if (dataset, algo, p) in stats]
            if points:
                ax.plot([p for p, _ in points], [t for _, t in points], marker="o", label=label)
        ax.set_title(f"Exact TBC / compute time\n{dataset}")
        ax.set_xlabel("Threads")
        ax.set_ylabel("Median compute time (s)")
        ax.set_xticks(all_threads)
        ax.grid(alpha=0.25)
        ax.legend()
        fig.tight_layout()
        fig.savefig(plots_dir / f"time_{slug}.png", dpi=175)
        plt.close(fig)
        plot_count += 1

        # A separate figure for speedup T(1)/T(p).
        fig, ax = plt.subplots(figsize=(9.0, 5.4))
        for algo, label in LABELS.items():
            one_thread = stats.get((dataset, algo, 1))
            if not one_thread:
                continue
            baseline = one_thread["median_time_s"]
            points = [(p, baseline / stats[(dataset, algo, p)]["median_time_s"])
                      for p in all_threads if (dataset, algo, p) in stats
                      and stats[(dataset, algo, p)]["median_time_s"] > 0]
            if points:
                ax.plot([p for p, _ in points], [s for _, s in points], marker="o", label=label)
        ax.plot(all_threads, all_threads, ls="--", lw=1, color="grey", label="Ideal")
        ax.set_title(f"Exact TBC / parallel speedup\n{dataset}")
        ax.set_xlabel("Threads")
        ax.set_ylabel("Speedup relative to 1 thread")
        ax.set_xticks(all_threads)
        ax.grid(alpha=0.25)
        ax.legend()
        fig.tight_layout()
        fig.savefig(plots_dir / f"speedup_{slug}.png", dpi=175)
        plt.close(fig)
        plot_count += 1

    print(f"Saved {plot_count} PNG figures in: {plots_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
