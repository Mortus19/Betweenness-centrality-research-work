#!/usr/bin/env python3
"""Exact TBC benchmark sweep: 3 C++ executables + TemporalGraphs.jl.

From the project root (Windows PowerShell):
    python benchmarks/run_exact_sweep.py --build --repeats 3 --threads 1-6

The runner does NOT download graph files. It recursively discovers every .patg
under the repository's Graphs/ directory and runs each experiment sequentially.
"""

from __future__ import annotations

import argparse
import csv
import fnmatch
import hashlib
import os
import platform
import random
import re
import shutil
import subprocess
import sys
import threading
import time
from datetime import datetime, timezone
from pathlib import Path

ALGOS = {
    "base": "TBC-base",
    "vertex": "TBC-optimized-only-vertex",
    "all": "TBC-optimized-all",
    "brunelli": None,
}
FIELDS = [
    "timestamp_utc", "git_commit", "dataset", "dataset_bytes", "algorithm",
    "threads", "repeat", "status", "compute_s", "preprocess_s",
    "total_without_startup_s", "wall_s", "peak_ram_mib", "peak_ram_source",
    "bc_sum_normalized", "return_code", "log_file", "error",
]


def parse_threads(spec: str) -> list[int]:
    result = []
    for chunk in spec.split(","):
        chunk = chunk.strip()
        if "-" in chunk:
            a, b = map(int, chunk.split("-", 1))
            if a > b:
                raise ValueError("Invalid threads range")
            result.extend(range(a, b + 1))
        else:
            result.append(int(chunk))
    if any(v < 1 for v in result):
        raise ValueError("Threads must be >= 1")
    return sorted(set(result))


def safe_name(s: str) -> str:
    cleaned = re.sub(r"[^a-zA-Z0-9_.-]+", "_", s)
    cleaned = cleaned[:68]
    return cleaned + "_" + hashlib.sha1(s.encode("utf-8")).hexdigest()[:9]


def repo_root() -> Path:
    return Path(__file__).resolve().parent.parent


def resolve_against(root: Path, path: str) -> Path:
    value = Path(path).expanduser()
    return (value if value.is_absolute() else root / value).resolve()


def get_git_commit(root: Path) -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"], cwd=root,
            stderr=subprocess.DEVNULL, text=True, timeout=5,
        ).strip()
    except (OSError, subprocess.SubprocessError):
        return "unknown"


def find_datasets(folder: Path, filters: list[str]) -> list[Path]:
    files = [p for p in folder.rglob("*") if p.is_file() and p.suffix.lower() == ".patg"]
    if filters:
        files = [p for p in files if any(
            fnmatch.fnmatch(p.name.lower(), pattern.lower()) or
            fnmatch.fnmatch(p.relative_to(folder).as_posix().lower(), pattern.lower())
            for pattern in filters
        )]
    return sorted(files, key=lambda f: (f.stat().st_size, f.as_posix()))


def check_build_type(build_dir: Path, allow_nonrelease: bool) -> None:
    cache = build_dir / "CMakeCache.txt"
    if not cache.is_file():
        raise RuntimeError(f"No CMakeCache.txt in {build_dir}. Configure Release profile in CLion first.")
    entries = {}
    for line in cache.read_text(encoding="utf-8", errors="replace").splitlines():
        if ":" in line and "=" in line:
            key = line.split(":", 1)[0]
            entries[key] = line.split("=", 1)[1]
    build_type = entries.get("CMAKE_BUILD_TYPE", "")
    multi = entries.get("CMAKE_CONFIGURATION_TYPES", "")
    if build_type.lower() != "release" and not multi and not allow_nonrelease:
        raise RuntimeError(
            f"CMake build type is {build_type or '(empty)'}, not Release. "
            "Switch to cmake-build-release or use --allow-nonrelease explicitly."
        )
    if build_type and build_type.lower() != "release":
        print(f"WARNING: CMake build type: {build_type}")


def parse_number(text: str, label_regex: str) -> float | None:
    match = re.search(r"^\s*" + label_regex + r"\s*:\s*([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)", text, re.MULTILINE | re.IGNORECASE)
    return float(match.group(1)) if match else None


def run_command(command: list[str], cwd: Path, timeout_s: float,
                env: dict[str, str]) -> tuple[int | None, str, float, float | None, str, str]:
    """Return (code, output, wall time, RSS MiB, RSS measurement kind, failure)."""
    try:
        process = subprocess.Popen(
            command, cwd=cwd, env=env,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, encoding="utf-8", errors="replace",
        )
    except OSError as exc:
        return None, "", 0.0, None, "unavailable", f"Cannot start process: {exc}"

    peak_memory = [0]
    kind = ["unavailable"]
    stop = threading.Event()

    def monitor() -> None:
        try:
            import psutil  # optional; pip install psutil
        except ImportError:
            return
        try:
            tracked = psutil.Process(process.pid)
        except psutil.Error:
            return
        while not stop.is_set():
            try:
                m = tracked.memory_info()
                if hasattr(m, "peak_wset"):
                    peak_memory[0] = max(peak_memory[0], m.peak_wset)
                    kind[0] = "windows_peak_working_set"
                else:
                    peak_memory[0] = max(peak_memory[0], m.rss)
                    kind[0] = "sampled_peak_rss"
            except psutil.Error:
                break
            stop.wait(0.05)

    sample_thread = threading.Thread(target=monitor, daemon=True)
    sample_thread.start()
    started = time.perf_counter()
    error = ""
    try:
        output, _ = process.communicate(timeout=timeout_s if timeout_s > 0 else None)
    except subprocess.TimeoutExpired:
        process.kill()
        output, _ = process.communicate()
        error = f"TIMEOUT after {timeout_s:g} s"
    finally:
        elapsed = time.perf_counter() - started
        stop.set()
        sample_thread.join(timeout=1)
    peak_mib = peak_memory[0] / (1024 * 1024) if peak_memory[0] else None
    return process.returncode, output or "", elapsed, peak_mib, kind[0], error


def load_done(csv_path: Path) -> set[tuple[str, str, str, str]]:
    if not csv_path.exists():
        return set()
    with csv_path.open("r", newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        return {
            (r["dataset"], r["algorithm"], r["threads"], r["repeat"])
            for r in reader if r.get("status") == "ok"
        }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--root", default=str(repo_root()), help="Repository root")
    parser.add_argument("--graphs-dir", default="Graphs", help="Directory recursively scanned for .patg")
    parser.add_argument("--build-dir", default="cmake-build-release")
    parser.add_argument("--output-dir", default="results/exact_sweep")
    parser.add_argument("--julia-script", default="benchmarks/benchmark_exact_julia.jl")
    parser.add_argument("--julia", default="julia", help="Path to Julia executable")
    parser.add_argument("--threads", default="1-6", help="1-6 or 1,2,4,6")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--algorithms", nargs="+", choices=list(ALGOS), default=list(ALGOS))
    parser.add_argument("--only", action="append", default=[], help="Filename/glob to select graphs; repeatable")
    parser.add_argument("--timeout", type=float, default=1800, help="Seconds per run, 0 disables timeout")
    parser.add_argument("--build", action="store_true", help="Build the three C++ Release targets first")
    parser.add_argument("--dry-run", action="store_true", help="List planned runs, without executing")
    parser.add_argument("--allow-nonrelease", action="store_true")
    parser.add_argument("--no-plot", action="store_true")
    parser.add_argument("--rerun", action="store_true", help="Ignore previously successful runs")
    args = parser.parse_args()

    if args.repeats < 1:
        parser.error("--repeats must be >= 1")
    try:
        thread_counts = parse_threads(args.threads)
    except ValueError as exc:
        parser.error(str(exc))

    root = Path(args.root).resolve()
    graphs_dir = resolve_against(root, args.graphs_dir)
    build_dir = resolve_against(root, args.build_dir)
    output_dir = resolve_against(root, args.output_dir)
    julia_script = resolve_against(root, args.julia_script)
    if not graphs_dir.is_dir():
        parser.error(f"Graphs directory not found: {graphs_dir}")
    datasets = find_datasets(graphs_dir, args.only)
    if not datasets:
        parser.error(f"No .patg graphs under: {graphs_dir}")
    print(f"Discovered {len(datasets)} temporal graph(s) under {graphs_dir}", flush=True)
    for graph in datasets:
        print(f"  {graph.relative_to(graphs_dir).as_posix():<55} {graph.stat().st_size / (1024 * 1024):>9.2f} MiB")

    algorithm_names = args.algorithms
    planned = len(datasets) * len(algorithm_names) * len(thread_counts) * args.repeats
    print(f"Plan: {len(datasets)} x {len(algorithm_names)} algorithms x {len(thread_counts)} thread counts x {args.repeats} repetitions = {planned} runs", flush=True)
    if args.dry_run:
        return 0

    if not julia_script.is_file() and "brunelli" in algorithm_names:
        parser.error(f"Julia benchmark script missing: {julia_script}")
    if "brunelli" in algorithm_names and not (shutil.which(args.julia) or Path(args.julia).is_file()):
        parser.error(f"Julia executable not found: {args.julia}")

    has_cpp = any(ALGOS[a] is not None for a in algorithm_names)
    if has_cpp:
        check_build_type(build_dir, args.allow_nonrelease)
        if args.build:
            targets = [ALGOS[a] for a in algorithm_names if ALGOS[a]]
            command = ["cmake", "--build", str(build_dir), "--config", "Release", "--target", *targets, "-j", "6"]
            print("Building:", subprocess.list2cmdline(command), flush=True)
            subprocess.run(command, cwd=root, check=True)

        for a in algorithm_names:
            if ALGOS[a] is None:
                continue
            exe = build_dir / (ALGOS[a] + (".exe" if os.name == "nt" else ""))
            if not exe.is_file():
                parser.error(f"Missing executable: {exe} (build Release targets first)")

    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / "runs.csv"
    logs = output_dir / "logs"
    values_dir = output_dir / "julia_values"
    logs.mkdir(exist_ok=True)
    values_dir.mkdir(exist_ok=True)
    already_done = set() if args.rerun else load_done(csv_path)
    append_header = not csv_path.exists() or csv_path.stat().st_size == 0
    commit = get_git_commit(root)
    print(f"Git commit: {commit} | machine: {platform.node()} | OS: {platform.platform()}")

    rng = random.Random(42)
    counter = 0
    successes = 0
    failures = 0
    with csv_path.open("a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS)
        if append_header:
            writer.writeheader()
            f.flush()

        for graph in datasets:
            dataset = graph.relative_to(root).as_posix() if graph.is_relative_to(root) else graph.name
            slug = safe_name(dataset)
            for rep in range(1, args.repeats + 1):
                tasks = [(a, p) for a in algorithm_names for p in thread_counts]
                rng.shuffle(tasks)  # reduce systematic thermal/order effects
                for algo, threads in tasks:
                    counter += 1
                    key = (dataset, algo, str(threads), str(rep))
                    if key in already_done:
                        print(f"[{counter}/{planned}] SKIP {dataset} / {algo} / p={threads} / r={rep}", flush=True)
                        continue
                    print(f"[{counter}/{planned}] {dataset} | {algo:8} | p={threads} | run={rep}", flush=True)
                    log_path = logs / f"{slug}__{algo}__p{threads}__r{rep}.log"
                    env = os.environ.copy()
                    env["OMP_NUM_THREADS"] = str(threads)
                    env["OMP_DYNAMIC"] = "FALSE"
                    env["JULIA_NUM_THREADS"] = str(threads)
                    if algo == "brunelli":
                        values = values_dir / f"{slug}.csv"
                        cmd = [args.julia, "--startup-file=no", f"--threads={threads}",
                               str(julia_script), str(graph), str(values)]
                    else:
                        exe = build_dir / (ALGOS[algo] + (".exe" if os.name == "nt" else ""))
                        cmd = [str(exe), str(threads), str(graph)]

                    code, output, wall, peak, peak_method, error = run_command(
                        cmd, root, args.timeout, env
                    )
                    log_path.write_text("COMMAND: " + subprocess.list2cmdline(cmd) + "\n\n" + output +
                                        ("\nERROR: " + error if error else ""), encoding="utf-8")
                    compute = parse_number(output, r"Compute time \(s\)" if algo == "brunelli" else r"Execution time \(s\)")
                    preprocessing = parse_number(output, r"Preprocessing \(s\)") if algo == "brunelli" else None
                    total = parse_number(output, r"Total \(s\)") if algo == "brunelli" else None
                    bc_sum = parse_number(output, r"BC sum normalized")
                    cpp_ram = parse_number(output, r"Peak RAM \(MiB\)") if algo != "brunelli" else None
                    if cpp_ram is not None:
                        peak, peak_method = cpp_ram, "cpp_peak_working_set"
                    ok = code == 0 and compute is not None and not error
                    if not ok and not error:
                        error = f"Exit code {code}, parsed compute time {compute}"
                    row = {
                        "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                        "git_commit": commit,
                        "dataset": dataset,
                        "dataset_bytes": graph.stat().st_size,
                        "algorithm": algo,
                        "threads": threads,
                        "repeat": rep,
                        "status": "ok" if ok else "failed",
                        "compute_s": compute if compute is not None else "",
                        "preprocess_s": preprocessing if preprocessing is not None else "",
                        "total_without_startup_s": total if total is not None else "",
                        "wall_s": round(wall, 6),
                        "peak_ram_mib": round(peak, 5) if peak is not None else "",
                        "peak_ram_source": peak_method,
                        "bc_sum_normalized": bc_sum if bc_sum is not None else "",
                        "return_code": code if code is not None else "",
                        "log_file": str(log_path.relative_to(output_dir)),
                        "error": error,
                    }
                    writer.writerow(row)
                    f.flush()
                    if ok:
                        successes += 1
                        print(f"       OK: compute={compute:.5f} s; peak RAM={row['peak_ram_mib'] or 'n/a'} MiB", flush=True)
                    else:
                        failures += 1
                        print(f"       FAILED: {error}; see {log_path}", flush=True)

    print(f"Complete: {successes} new successful runs, {failures} failed, results in {csv_path}", flush=True)
    if not args.no_plot:
        plotter = Path(__file__).with_name("plot_exact_sweep.py")
        if plotter.is_file():
            try:
                subprocess.run([sys.executable, str(plotter), str(csv_path)], cwd=root, check=True)
            except subprocess.CalledProcessError as exc:
                print(f"Plot step failed (results are still saved): {exc}", file=sys.stderr)
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
