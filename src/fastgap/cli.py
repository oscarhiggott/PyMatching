"""Command line tools: ``python -m fastgap {bench,build-index,landmark-report,serve}``."""

from __future__ import annotations

import argparse
import csv
import json
import os
import struct
import subprocess
import sys
import time
from typing import List, Optional

import numpy as np

from fastgap import _cpp_fastgap
from fastgap.decoder import GapDecoder
from fastgap.exact import exact_gap_batch
from fastgap.index import GapIndex

CSV_COLUMNS = [
    "shot",
    "num_defects",
    "w_star",
    "gap_exact",
    "gap_lb",
    "gap_ub",
    "t_decode_ns",
    "t_gap_ns",
    "t_exact_ns",
    "walk_simple",
    "all_hops_exact",
    "censored",
    "strategy_hash",
]


# ---------------------------------------------------------------------------------------------
# Inputs


def add_dem_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--dem", help="path to a graphlike .dem file")
    p.add_argument("--gen", help="generate a stim circuit instead, e.g. surface_code:rotated_memory_x")
    p.add_argument("-d", "--distance", type=int, default=5)
    p.add_argument("-r", "--rounds", type=int, default=None, help="default: distance")
    p.add_argument("-p", "--noise", type=float, default=1e-3, help="uniform circuit noise strength")
    p.add_argument("--observable", type=int, default=0)


def load_dem(args):
    import stim

    if args.dem:
        return stim.DetectorErrorModel.from_file(args.dem)
    if args.gen:
        p = args.noise
        circuit = stim.Circuit.generated(
            args.gen,
            distance=args.distance,
            rounds=args.rounds or args.distance,
            after_clifford_depolarization=p,
            before_round_data_depolarization=p,
            before_measure_flip_probability=p,
            after_reset_flip_probability=p,
        )
        return circuit.detector_error_model(decompose_errors=True)
    raise SystemExit("pass --dem FILE or --gen CODE_TASK")


def add_index_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--idx", help="prebuilt index file (from build-index)")
    p.add_argument("--strategy", default="auto")
    p.add_argument("--num-landmarks", type=int, default=16)


def get_index(args, dem) -> GapIndex:
    if getattr(args, "idx", None) and os.path.exists(args.idx):
        return GapIndex.load(args.idx)
    return GapIndex.from_dem(
        dem,
        observable=args.observable,
        strategy=args.strategy,
        num_landmarks=args.num_landmarks,
    )


def sample_shots(dem, n: int, seed: Optional[int]) -> np.ndarray:
    dets, _, _ = dem.compile_sampler(seed=seed).sample(n, bit_packed=True)
    return np.ascontiguousarray(dets)


def percentiles(x) -> dict:
    x = np.asarray(x, dtype=np.float64)
    if len(x) == 0:
        return {}
    return {
        "median": float(np.median(x)),
        "p99": float(np.percentile(x, 99)),
        "p99.9": float(np.percentile(x, 99.9)),
        "mean": float(np.mean(x)),
        "max": float(np.max(x)),
    }


# ---------------------------------------------------------------------------------------------
# build-index


def cmd_build_index(args) -> int:
    dem = load_dem(args)
    t0 = time.perf_counter()
    idx = GapIndex.from_dem(
        dem,
        observable=args.observable,
        strategy=args.strategy,
        num_landmarks=args.num_landmarks,
    )
    elapsed = time.perf_counter() - t0
    idx.save(args.out)
    report = idx.memory_report()
    report.update({"strategy": idx.strategy_name, "strategy_hash": idx.strategy_hash, "build_seconds": elapsed,
                   "landmarks": [l.label for l in idx.landmarks], "out": args.out})
    print(json.dumps(report, indent=2))
    return 0


# ---------------------------------------------------------------------------------------------
# serve: a minimal pipe protocol used to measure end-to-end latency (syndrome arriving from
# outside the process, answer leaving it).
#   request:  uint32 length, then `length` bytes of one bit-packed shot (length 0 = quit)
#   response: uint64 prediction mask, int64 gap_lb (internal units, -1 = inf), uint8 censored

RESPONSE = struct.Struct("<Qqb")


def cmd_serve(args) -> int:
    idx = GapIndex.load(args.idx)
    dec = GapDecoder(idx, num_threads=args.threads)
    if args.pin is not None:
        _cpp_fastgap.pin_current_thread(args.pin)
    inp = sys.stdin.buffer
    out = sys.stdout.buffer
    core = dec._core
    threshold = -1 if args.threshold is None else int(np.ceil(args.threshold * idx.normalising_constant))
    while True:
        header = inp.read(4)
        if len(header) < 4:
            return 0
        (n,) = struct.unpack("<I", header)
        if n == 0:
            return 0
        shot = np.frombuffer(inp.read(n), dtype=np.uint8).reshape(1, n)
        r = core.decode_batch(shot, threshold=threshold, upper_bound=False, pricing="index", check=False)
        out.write(RESPONSE.pack(int(r["prediction_mask"][0]), int(r["gap_lb"][0]), int(r["censored"][0])))
        out.flush()


def measure_end_to_end(idx_path: str, shots: np.ndarray, threshold: Optional[float], pin: Optional[int],
                       threads: int = 1) -> np.ndarray:
    cmd = [sys.executable, "-m", "fastgap", "serve", "--idx", idx_path, "--threads", str(threads)]
    if threshold is not None:
        cmd += ["--threshold", str(threshold)]
    if pin is not None:
        cmd += ["--pin", str(pin)]
    proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, bufsize=0)
    times = np.zeros(len(shots), dtype=np.int64)
    try:
        for i, row in enumerate(shots):
            payload = struct.pack("<I", len(row)) + row.tobytes()
            t0 = time.perf_counter_ns()
            proc.stdin.write(payload)
            got = b""
            while len(got) < RESPONSE.size:
                chunk = proc.stdout.read(RESPONSE.size - len(got))
                if not chunk:
                    raise RuntimeError("fastgap serve exited early")
                got += chunk
            times[i] = time.perf_counter_ns() - t0
        proc.stdin.write(struct.pack("<I", 0))
    finally:
        proc.stdin.close()
        proc.wait(timeout=30)
    return times


# ---------------------------------------------------------------------------------------------
# bench


def cmd_bench(args) -> int:
    dem = load_dem(args)
    idx = get_index(args, dem)
    if args.pin is not None:
        _cpp_fastgap.pin_current_thread(args.pin)
    shots = sample_shots(dem, args.shots, args.seed)
    warm = sample_shots(dem, args.warmup, None if args.seed is None else args.seed + 1)

    # Latency at each thread count. Shots are always decoded one after another; the threads of a
    # decoder cooperate on each shot's gap search, so more threads lower per-shot latency.
    pin = args.pin is not None
    thread_counts = sorted(set([1] + list(args.threads)))
    by_threads = {}
    for t in thread_counts:
        dec = GapDecoder(idx, num_threads=t, pin_threads=pin)
        dec.decode_batch(warm, threshold=args.threshold, upper_bound=args.upper_bound)
        by_threads[t] = dec.decode_batch(shots, threshold=args.threshold, upper_bound=args.upper_bound)
    res = by_threads[1]
    for t, res_t in by_threads.items():
        same = (np.array_equal(res.gap_lb_int, res_t.gap_lb_int) and np.array_equal(res.gap_ub_int, res_t.gap_ub_int)
                and np.array_equal(res.prediction, res_t.prediction) and np.array_equal(res.walk_simple, res_t.walk_simple))
        if not same:
            raise AssertionError(f"{t}-thread results differ from single-threaded results")

    exact = None
    violations = 0
    if not args.no_exact:
        exact_gap_batch(idx, warm)
        exact = exact_gap_batch(idx, shots)
        bad = res.gap_lb_int > exact.gap_int
        violations = int(np.sum(bad))
        ub_bad = int(np.sum(res.walk_simple & (res.gap_ub_int < exact.gap_int)))
        if violations or ub_bad:
            print(f"SOUNDNESS VIOLATION: {violations} shots with gap_lb > gap_exact, {ub_bad} with gap_ub < gap_exact",
                  file=sys.stderr)

    e2e = None
    if args.e2e:
        path = args.idx
        if not path:
            path = os.path.join(os.path.dirname(os.path.abspath(args.out or ".")) or ".", ".fastgap_bench.idx")
            idx.save(path)
        e2e = {}
        for t in thread_counts:
            measure_end_to_end(path, warm[: min(len(warm), 200)], args.threshold, args.pin, t)
            e2e[t] = measure_end_to_end(path, shots[: args.e2e_shots], args.threshold, args.pin, t)

    nc = idx.normalising_constant
    if args.out:
        with open(args.out, "w", newline="") as f:
            w = csv.writer(f)
            extra = [t for t in thread_counts if t != 1]
            w.writerow(CSV_COLUMNS + [f"t_gap_ns_t{t}" for t in extra] + ["t_exact_mstar_ns", "t_exact_comp_ns"])
            for i in range(len(res)):
                w.writerow([
                    i,
                    int(res.num_defects[i]),
                    repr(float(res.weight[i])),
                    "" if exact is None else repr(float(exact.gap[i])),
                    repr(float(res.gap_lb[i])),
                    repr(float(res.gap_ub[i])),
                    int(res.t_decode_ns[i]),
                    int(res.t_gap_ns[i]),
                    "" if exact is None else int(exact.t_exact_ns[i]),
                    int(res.walk_simple[i]),
                    int(res.all_hops_exact[i]),
                    int(res.censored[i]),
                    idx.strategy_hash,
                ] + [int(by_threads[t].t_gap_ns[i]) for t in extra] + (
                    ["", ""] if exact is None else [int(exact.t_mstar_ns[i]), int(exact.t_comp_ns[i])]))

    n = len(res)
    summary = {
        "num_shots": n,
        "num_detectors": idx.num_detectors,
        "strategy": idx.strategy_name,
        "strategy_hash": idx.strategy_hash,
        "num_landmarks": len(idx.landmarks),
        "normalising_constant": nc,
        "threads": thread_counts,
        "pinned": args.pin is not None and sys.platform.startswith("linux"),
        # Compute-only per-shot latency at each thread count (stage timers inside decode_batch).
        "latency_ns_by_threads": {
            str(t): {
                "decode": percentiles(r.t_decode_ns),
                "gap": percentiles(r.t_gap_ns),
                "total": percentiles(r.t_decode_ns + r.t_gap_ns),
                "gap_median_speedup": float(np.median(res.t_gap_ns) / max(1.0, np.median(r.t_gap_ns))),
            }
            for t, r in by_threads.items()
        },
        "gap_overhead_ratio_median": float(np.median(res.t_gap_ns) / max(1.0, np.median(res.t_decode_ns))),
        "walk_simple_fraction": float(np.mean(res.walk_simple)),
        "all_hops_exact_fraction": float(np.mean(res.all_hops_exact)),
        "censored_fraction": float(np.mean(res.censored)),
        "blossom_shot_fraction": float(np.mean(res.num_blossoms > 0)),
        "mean_defects": float(np.mean(res.num_defects)),
    }
    if exact is not None:
        summary["exact_ns"] = percentiles(exact.t_exact_ns)
        # The baseline's two decodes separately: the one in the class of M* plays the role of the
        # ordinary decode; the complementary-class decode is what the baseline pays for the gap.
        summary["exact_mstar_decode_ns"] = percentiles(exact.t_mstar_ns)
        summary["exact_comp_decode_ns"] = percentiles(exact.t_comp_ns)
        # Gap computation only: baseline complementary decode latency / fastgap gap-stage latency.
        summary["gap_speedup_vs_exact_by_threads"] = {
            str(t): {
                q: float(summary["exact_comp_decode_ns"][q] / max(1.0, summary["latency_ns_by_threads"][str(t)]["gap"][q]))
                for q in ("median", "p99", "p99.9")
            }
            for t in by_threads
        }
        summary["soundness_violations"] = violations
        summary["gap_lb_equals_exact_fraction"] = float(np.mean(res.gap_lb_int == exact.gap_int))
        summary["mean_gap_shortfall_db"] = float(np.mean(GapIndex.to_db(np.where(
            np.isinf(exact.gap), 0.0, exact.gap - np.where(np.isinf(res.gap_lb), exact.gap, res.gap_lb)))))
        if args.threshold is not None:
            cls = res.classify(args.threshold)
            summary["uncertain_fraction_at_threshold"] = float(np.mean(cls == 0))
    if e2e is not None:
        summary["end_to_end_ns_by_threads"] = {str(t): percentiles(v) for t, v in e2e.items()}
    text = json.dumps(summary, indent=2)
    print(text)
    if args.summary:
        with open(args.summary, "w") as f:
            f.write(text + "\n")
    return 1 if violations else 0


# ---------------------------------------------------------------------------------------------
# landmark-report (CLAUDE.md §7.5)


def cmd_landmark_report(args) -> int:
    dem = load_dem(args)
    idx = get_index(args, dem)
    shots = sample_shots(dem, args.shots, args.seed)
    dec = GapDecoder(idx)
    t0 = time.perf_counter()
    res = dec.decode_batch(shots)
    wall = time.perf_counter() - t0
    sigma = dec.decode_batch(shots, pricing="exact")
    exact = exact_gap_batch(idx, shots)
    core = idx._core

    # LB/dist on the hop pairs of real walks.
    ratios: List[float] = []
    dense = np.unpackbits(shots, axis=1, bitorder="little")[:, : idx.num_detectors]
    for row in dense[: args.hop_shots]:
        r = dec.decode(row, keep_walk=True)
        for a, b, _, _exact in r.walk:
            if a < 0 or b < 0:
                continue
            d = core.exact_distance(a, b)
            if d > 0:
                ratios.append(core.lower_bound(a, b) / d)

    finite = np.isfinite(exact.gap_int) & np.isfinite(res.gap_lb_int) & np.isfinite(sigma.gap_lb_int)
    nc = idx.normalising_constant
    lb_part = GapIndex.to_db((sigma.gap_lb_int - res.gap_lb_int)[finite] / nc)
    alt_part = GapIndex.to_db((exact.gap_int - sigma.gap_lb_int)[finite] / nc)
    total = GapIndex.to_db((exact.gap_int - res.gap_lb_int)[finite] / nc)
    edges = np.array([0, 1e-9, 0.5, 1, 2, 3, 5, 10, 20, np.inf])

    def hist(x):
        counts, _ = np.histogram(x, bins=edges)
        return {f"[{edges[i]:g},{edges[i + 1]:g})": int(c) for i, c in enumerate(counts)}

    ratios_arr = np.asarray(ratios)
    report = {
        "dem_detectors": idx.num_detectors,
        "strategy": idx.strategy_name,
        "strategy_hash": idx.strategy_hash,
        "landmarks": [l.label for l in idx.landmarks],
        "memory": idx.memory_report(),
        "num_shots": len(res),
        "lb_over_dist_on_hops": {
            "count": int(len(ratios_arr)),
            "mean": float(ratios_arr.mean()) if len(ratios_arr) else None,
            "fraction_exact": float(np.mean(ratios_arr >= 1.0)) if len(ratios_arr) else None,
            "quantiles": {q: float(np.quantile(ratios_arr, q)) for q in (0.01, 0.1, 0.5, 0.9)} if len(ratios_arr) else {},
        },
        "fraction_gap_lb_equals_delta_sigma": float(np.mean(res.gap_lb_int == sigma.gap_lb_int)),
        "fraction_gap_lb_equals_gap_exact": float(np.mean(res.gap_lb_int == exact.gap_int)),
        "shortfall_db_histogram": {"total": hist(total), "lb_part": hist(lb_part), "alternation_blossom_part": hist(alt_part)},
        "us_per_shot": 1e6 * wall / len(res),
        "soundness_violations": int(np.sum(res.gap_lb_int > exact.gap_int)),
    }
    text = json.dumps(report, indent=2)
    print(text)
    if args.out:
        with open(args.out, "w") as f:
            f.write(text + "\n")
    return 1 if report["soundness_violations"] else 0


# ---------------------------------------------------------------------------------------------


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m fastgap", description=__doc__)
    sub = parser.add_subparsers(dest="cmd", required=True)

    b = sub.add_parser("bench", help="per-shot CSV and latency summary")
    add_dem_args(b)
    add_index_args(b)
    b.add_argument("--shots", type=int, default=10000)
    b.add_argument("--warmup", type=int, default=1000)
    b.add_argument("--threads", type=int, nargs="+", default=[1],
                   help="thread counts per shot to measure latency at (1 is always included), e.g. --threads 4 10")
    b.add_argument("--seed", type=int, default=None)
    b.add_argument("--threshold", type=float, default=None, help="post-selection threshold in nats")
    b.add_argument("--upper-bound", action="store_true", default=True)
    b.add_argument("--no-upper-bound", dest="upper_bound", action="store_false")
    b.add_argument("--no-exact", action="store_true", help="skip the exact two-decode baseline")
    b.add_argument("--e2e", action="store_true", help="also measure end-to-end latency through a pipe")
    b.add_argument("--e2e-shots", type=int, default=2000)
    b.add_argument("--pin", type=int, default=None, help="pin threads starting at this CPU (Linux)")
    b.add_argument("--out", help="per-shot CSV path")
    b.add_argument("--summary", help="summary JSON path")
    b.set_defaults(func=cmd_bench)

    bi = sub.add_parser("build-index", help="build and save a GapIndex")
    add_dem_args(bi)
    bi.add_argument("--strategy", default="auto")
    bi.add_argument("--num-landmarks", type=int, default=16)
    bi.add_argument("--out", required=True)
    bi.set_defaults(func=cmd_build_index)

    lr = sub.add_parser("landmark-report", help="tightness report for a DEM and landmark strategy")
    add_dem_args(lr)
    add_index_args(lr)
    lr.add_argument("--shots", type=int, default=10000)
    lr.add_argument("--hop-shots", type=int, default=2000, help="shots whose walks are used for LB/dist")
    lr.add_argument("--seed", type=int, default=None)
    lr.add_argument("--out", help="report JSON path")
    lr.set_defaults(func=cmd_landmark_report)

    s = sub.add_parser("serve", help="decode shots from stdin (used by bench --e2e)")
    s.add_argument("--idx", required=True)
    s.add_argument("--threshold", type=float, default=None)
    s.add_argument("--pin", type=int, default=None)
    s.add_argument("--threads", type=int, default=1, help="threads per shot")
    s.set_defaults(func=cmd_serve)

    args = parser.parse_args(argv)
    return int(args.func(args) or 0)
