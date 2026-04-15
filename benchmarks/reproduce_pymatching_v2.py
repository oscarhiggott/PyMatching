import time
import csv
import numpy as np
import matplotlib.pyplot as plt
import stim
import pymatching

DISTANCES = [5, 7, 9, 13, 17]
NOISE_RATES = [0.003, 0.005, 0.007]
NUM_SHOTS = 100000

def benchmark_surface_code(distance: int, p: float, num_shots: int):
    circuit = stim.Circuit.generated(
        "surface_code:rotated_memory_x",
        distance=distance,
        rounds=distance,
        after_clifford_depolarization=p,
        before_round_data_depolarization=p,
        before_measure_flip_probability=p,
        after_reset_flip_probability=p,
    )

    dem = circuit.detector_error_model(decompose_errors=True)
    matching = pymatching.Matching.from_detector_error_model(dem)

    sampler = circuit.compile_detector_sampler()
    shots, actual_observables = sampler.sample(shots=num_shots, separate_observables=True)

    matching.decode_batch(shots[0:1, :])  # warmup

    t0 = time.time()
    predictions = matching.decode_batch(shots)
    t1 = time.time()

    num_errors = int(np.sum(np.any(predictions != actual_observables, axis=1)))
    logical_error_rate = num_errors / num_shots

    return {
        "distance": distance,
        "p": p,
        "num_shots": num_shots,
        "num_detectors": circuit.num_detectors,
        "detectors_per_round": circuit.num_detectors / distance,
        "seconds_per_shot": (t1 - t0) / num_shots,
        "seconds_per_round": (t1 - t0) / num_shots / distance,
        "num_errors": num_errors,
        "logical_error_rate": logical_error_rate,
    }

# --- Run benchmarks ---
rows = []
for p in NOISE_RATES:
    for d in DISTANCES:
        row = benchmark_surface_code(d, p, NUM_SHOTS)
        rows.append(row)
        print(f"p={p:.3f} d={d:2d}  err_rate={row['logical_error_rate']:.5f}  "
              f"time/round={row['seconds_per_round']*1e6:.2f}μs  "
              f"errors={row['num_errors']}/{NUM_SHOTS}")

filename_base = "benchmarks/pymatching_v2_baseline"

with open(f"{filename_base}.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)

# --- Fidelity plot ---
plt.figure(figsize=(6, 4))
for p in NOISE_RATES:
    subset = [r for r in rows if r["p"] == p]
    x = [r["distance"] for r in subset]
    y = [r["logical_error_rate"] for r in subset]
    plt.semilogy(x, y, marker="s", label=f"p={p}")

plt.xlabel("Distance")
plt.ylabel("Logical error rate per shot")
plt.title(f"PyMatching v2 (K=1) fidelity baseline, {NUM_SHOTS} shots")
plt.legend()
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.tight_layout()
plt.savefig(f"{filename_base}_fidelity.png", dpi=200)
plt.close()

# --- Runtime plot (use p=0.005 as representative) ---
subset = [r for r in rows if r["p"] == 0.005]
x = np.array([r["distance"] for r in subset], dtype=float)
y = np.array([r["seconds_per_round"] for r in subset], dtype=float)

plt.figure(figsize=(6, 4))
plt.loglog(x, y, marker="o", label="PyMatching v2 (K=1)")
plt.xlabel("Distance")
plt.ylabel("Decoding time per round (s)")
plt.title("PyMatching v2 runtime, p=0.005")
plt.legend()
plt.tight_layout()
plt.savefig(f"{filename_base}_runtime.png", dpi=200)
plt.close()

print(f"\nSaved: {filename_base}_fidelity.png, {filename_base}_runtime.png")
