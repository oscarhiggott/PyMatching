import time
import csv
import numpy as np
import matplotlib.pyplot as plt
import stim
import pymatching

DISTANCES = [5, 7, 9, 13, 17, 23, 29, 39, 50]
P = 0.001
NUM_SHOTS = 10000

def time_surface_code_circuit(distance: int, p: float, num_shots: int = 10000):
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
    shots, actual_observables = sampler.sample(
        shots=num_shots,
        separate_observables=True,
    )

    # warmup so internal structures are cached
    matching.decode_batch(shots[0:1, :])

    t0 = time.time()
    predictions = matching.decode_batch(shots)
    t1 = time.time()

    seconds_per_shot = (t1 - t0) / num_shots
    seconds_per_round = seconds_per_shot / distance

    num_detectors = circuit.num_detectors
    detectors_per_round = num_detectors / distance

    return {
        "distance": distance,
        "p": p,
        "num_shots": num_shots,
        "num_detectors": num_detectors,
        "detectors_per_round": detectors_per_round,
        "seconds_per_shot": seconds_per_shot,
        "seconds_per_round": seconds_per_round,
    }

rows = []
for d in DISTANCES:
    row = time_surface_code_circuit(d, P, NUM_SHOTS)
    rows.append(row)
    print(row)

filename_base = f"benchmarks/pymatching_v2_reproduced_p_{P}"

with open(f"{filename_base}.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)

# Re-extract rows from the saved CSV file
with open(f"{filename_base}.csv", "r", newline="") as f:
    reader = csv.DictReader(f)
    rows = []
    for row in reader:
        rows.append({
            "distance": int(row["distance"]),
            "p": float(row["p"]),
            "num_shots": int(row["num_shots"]),
            "num_detectors": int(row["num_detectors"]),
            "detectors_per_round": float(row["detectors_per_round"]),
            "seconds_per_shot": float(row["seconds_per_shot"]),
            "seconds_per_round": float(row["seconds_per_round"]),
        })

# Plot against distance
x = np.array([r["distance"] for r in rows], dtype=float)
y = np.array([r["seconds_per_round"] for r in rows], dtype=float)

plt.figure(figsize=(6, 4))
plt.loglog(x, y, marker="o", label="PyMatching v2")

# Fit using detectors_per_round, but still plot the fit versus distance
fit_mask = np.array([r["distance"] > 10 for r in rows])
fit_detectors = np.array([r["detectors_per_round"] for r in rows], dtype=float)[fit_mask]
fit_times = y[fit_mask]

coef = np.polyfit(np.log(fit_detectors), np.log(fit_times), 1)
slope, intercept = coef

all_detectors = np.array([r["detectors_per_round"] for r in rows], dtype=float)
y_fit = np.exp(intercept) * all_detectors**slope

# sort by distance so the fitted line is drawn cleanly
sort_idx = np.argsort(x)
plt.loglog(
    x[sort_idx],
    y_fit[sort_idx],
    "--",
    label=fr"fit: T $\propto N^{{{slope:.2f}}}$"
)

plt.xlabel("Distance")
plt.ylabel("decoding time per round (s)")
plt.title(f"PyMatching v2 benchmark, p={P}")
plt.legend()
plt.tight_layout()
plt.savefig(f"{filename_base}.png", dpi=200)
plt.show()