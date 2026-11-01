# fastgap

This fork adds **fastgap**, a Python package shipped in the same distribution as `pymatching`.
It decodes each shot with the ordinary sparse-blossom decode and returns a
**certified lower bound on the complementary gap** `Δ = W^c − W*` (the weight difference between
the best matching in the other logical class and the minimum-weight matching). The bound comes
from the decoder's own final dual state. fastgap never runs a second decode and never sweeps the
empty spacetime volume.

Why: the complementary gap is the standard per-shot confidence signal for post-selection in magic-state cultivation
and decoder confidence in yoked surface code and other code concatenation schemes.
The usual way to compute it is to decode a second time, forced into the other logical class.
That second decode costs as much as the first, or more, and grows with the spacetime volume.
fastgap replaces it with a small search over the defects alone.
The search is sound by construction: `gap_lb ≤ Δ` on every shot, checked against the exact
two-decode baseline on 10⁵ shots per test configuration. Predictions and
`W*` stay bit-identical to upstream PyMatching.

```python
import stim, fastgap

circuit = stim.Circuit.generated("surface_code:rotated_memory_x", distance=7, rounds=7,
                                 after_clifford_depolarization=1e-3)
dem = circuit.detector_error_model(decompose_errors=True)

idx = fastgap.GapIndex.from_dem(dem, observable=0, strategy="auto", num_landmarks=16)
idx.save("rmx_d7.idx")                      # cache the index per DEM
dec = fastgap.GapDecoder(idx, num_threads=4)  # 4 threads share each shot's gap search

shots = circuit.compile_detector_sampler().sample(10_000, bit_packed=True)
res = dec.decode_batch(shots, threshold=None, upper_bound=True)
res.prediction    # identical to pymatching.Matching.decode_batch
res.weight        # W*, identical to upstream
res.gap_lb        # certified lower bound on the gap (nats); res.gap_lb_db in dB
res.gap_ub        # upper bound when the returned walk is simple, else inf

prediction, w_star, gap = fastgap.exact_gap_batch(idx, shots)   # exact two-decode baseline
assert (res.gap_lb_int <= fastgap.exact_gap_batch(idx, shots).gap_int).all()
```

For post-selection at threshold `τ`, pass `threshold=τ`: the search stops as soon as the bound
reaches `τ` (`censored=True`), and `res.classify(τ)` sorts shots into keep, discard and uncertain.

How it works:

1. Observable `k` is gauge-fixed onto the boundary, which splits the boundary into sides `L` and `R`.
2. The radii `y_v`, `y_S` are read before shattering the blossoms.
3. A dense Dijkstra runs over alternating "hop-relay-hop-relay" walks from `L` to `R`,
    priced by reduced costs `σ(u,v) = dist(u,v) − ρ_u − ρ_v + 2β(u,v)`.
    Each hop distance is a lower bound taken from landmark potentials.
4. `num_threads = T` threads cooperate on *each* shot to lower its latency; shots are decoded one
   after another. The calling thread runs the Dijkstra, and the other `T − 1` threads price, in
   advance, the rows of the states it is about to settle. Results are bit-identical for every
   `T`. Use at most the size of the fastest core cluster (e.g. 4–5 on an Apple M5 Pro).

### Gap-computation latency

The baseline runs two PyMatching decodes with the observable as an extra detector, toggled off and
on. Both methods run one ordinary decode, which gives the prediction (min weight matching).
The baseline then pays a second, complementary-class decode for the gap, while fastgap pays its gap search,
(without `gap_ub`). The table compares only those two costs:
speed-up = baseline complementary decode / fastgap gap stage.

Setup: rotated surface-code memory X, `d` rounds, uniform circuit noise `p`, 2×10⁴ shots per DEM,
compute-only per-shot timers (median / p99), Apple M5 Pro. There were 0 soundness violations.

| DEM | mean defects | baseline gap (complementary decode) | fastgap gap, T=1 | fastgap gap, T=4 | speed-up T=1 (med / p99) | speed-up T=4 (med / p99) |
|---|---|---|---|---|---|---|
| d=5, p=5e-3 | 8.4 | 6.92 / 14.0 µs | 0.62 / 1.7 µs | 0.62 / 2.3 µs | 11.1× / 8.2× | 11.1× / 6.0× |
| d=7, p=5e-3 | 25.3 | 23.92 / 42.9 µs | 1.67 / 4.0 µs | 2.42 / 4.9 µs | 14.3× / 10.7× | 9.9× / 8.8× |
| d=5, p=1e-3 | 1.8 | 3.50 / 9.5 µs | 0.12 / 0.8 µs | 0.12 / 0.7 µs | 28.0× / 12.6× | 28.0× / 13.3× |
| d=9, p=1e-3 | 12.0 | 37.62 / 60.7 µs | 0.83 / 2.2 µs | 0.92 / 3.0 µs | 45.1× / 27.0× | 41.1× / 20.0× |
| d=13, p=1e-3 | 38.3 | 154.46 / 226.9 µs | 3.04 / 6.7 µs | 3.62 / 6.8 µs | 50.8× / 33.8× | 42.6× / 33.6× |
| d=15, p=1e-3 | 59.6 | 286.17 / 421.0 µs | 5.58 / 11.5 µs | 5.46 / 10.2 µs | 51.3× / 36.5× | 52.4× / 41.2× |
| d=17, p=1e-3 | 87.8 | 518.92 / 787.9 µs | 10.12 / 19.6 µs | 8.54 / 15.5 µs | 51.3× / 40.1× | 60.7× / 51.0× |

The shared ordinary decode costs about the same in both methods (e.g. 7.8 µs for fastgap vs
8.4 µs for the baseline at d=17), so it is left out. Raw per-shot CSVs and summaries are in
`benchmarks/fastgap/results/`; more detail is in `NOTES.md`.

Landmark placement is pluggable (`fastgap.landmarks`). It changes only how tight the bound is,
never whether it is valid.

Command line:

```
python -m fastgap build-index --dem X.dem --strategy auto --out X.idx
python -m fastgap bench --dem X.dem --idx X.idx --shots 100000 --threads 4 --out X.csv --e2e
python -m fastgap landmark-report --dem X.dem --idx X.idx --shots 10000
```

Supported in v1: graphlike DEMs with boundaries, non-negative weights, and one observable per
index. Not supported: correlated matching, and codes with a logical that never touches a
boundary (toric codes raise `fastgap.BoundarylessLogicalError`). The design document is
`CLAUDE.md`. Implementation notes and measured results are in `NOTES.md`, and benchmarks are in
`benchmarks/fastgap/`.

## Attribution

Cooked with Claude Opus 5.5. Paper for `fastgap` is coming soon.

Please also cite the following as our method is based on it.

```
@article{Higgott2025sparseblossom,
  doi = {10.22331/q-2025-01-20-1600},
  url = {https://doi.org/10.22331/q-2025-01-20-1600},
  title = {Sparse {B}lossom: correcting a million errors per core second with minimum-weight matching},
  author = {Higgott, Oscar and Gidney, Craig},
  journal = {{Quantum}},
  issn = {2521-327X},
  publisher = {{Verein zur F{\"{o}}rderung des Open Access Publizierens in den Quantenwissenschaften}},
  volume = {9},
  pages = {1600},
  month = jan,
  year = {2025}
}
```
