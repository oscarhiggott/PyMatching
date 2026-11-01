# fastgap implementation notes

Fork point: upstream PyMatching `6f63b2b9474ba0fa7e511fe52bffdce858a06984`. The file and function
names in CLAUDE.md §5.1 were checked against this commit and are unchanged.

Environment used for the benchmark numbers below: Apple M5 Pro (arm64; 5 "Super" cores sharing
one L2, plus 10 performance cores in two L2 clusters), macOS 26.5, Apple clang 21, Python 3.14,
stim 1.16. Older sections (speed history, soundness tables before the re-run) were measured on an
earlier Apple M-series laptop with Apple clang 16 and Python 3.12. The timer is `std::chrono::steady_clock`, which ticks in
steps of ~41.7 ns on this machine, so sub-microsecond medians are quantised. macOS has no hard
thread-affinity API, so `--pin` is a no-op here; pinning only takes effect on Linux.

## Layout

| Path | What |
|---|---|
| `src/pymatching/fastgap/canonical.{h,cc}` | §2.3 gauge fixing on the raw flattened DEM, the parallel-edge check, canonical and baseline DEMs/graphs |
| `src/pymatching/fastgap/interior_graph.{h,cc}` | CSR copy of the decoder's own `MatchingGraph` (integer weights, sides L/R), Dijkstra, A* |
| `src/pymatching/fastgap/dual_state.{h,cc}` | §5.1 `extract_dual_state`: dual radii before shattering, then the upstream shatter for `mate` |
| `src/pymatching/fastgap/gap_index.{h,cc}` | §5.3 index: `h_L`, `h_R`, landmark potentials, exact local tables, save/load |
| `src/pymatching/fastgap/gap_decoder.{h,cc}` | §3.3–3.6 dense Dijkstra over hop/relay states, the upper bound, diagnostics, batch decoding |
| `src/pymatching/fastgap/exact.{h,cc}` | §5.2 two-decode baseline and the brute-force subset DP |
| `src/pymatching/fastgap/thread_pool.h` | `SpinTeam` (the threads that share one shot's search) and `ThreadPool` (shot-parallel exact baseline) |
| `src/pymatching/fastgap/fastgap.pybind.cc` | `fastgap._cpp_fastgap` (internal integer units only) |
| `src/fastgap/` | Python API (`GapIndex`, `GapDecoder`, `exact_gap_batch`), landmark strategies, CLI |
| `tests/fastgap/` | §9 tests in Python; C++ tests are the `*.test.cc` files next to each source |
| `benchmarks/fastgap/` | DEM generator, index config, run script, raw CSVs and summaries |

The only upstream edit is the pair of declarations added to `driver/mwpm_decoding.h`. Build
plumbing changed too: `CMakeLists.txt` adds `FASTGAP_*` source lists, the `_cpp_fastgap` module
and `-DFASTGAP_CHECK` on the test target; `setup.py` builds both extensions from one CMake tree;
`BUILD` adds `-DFASTGAP_CHECK` to `cc_test`. Upstream tests still pass: 116 upstream C++ tests,
92 upstream Python tests (13 skipped by upstream for optional dependencies).

Note for macOS (case-insensitive file systems): `cmake -B build` collides with Bazel's `BUILD`
file. Use another directory name, e.g. `cmake -B cmake-build` (already git-ignored).

## Local tables removed (branch `remove-local-tables`, 2026-11-01)

The exact local distance tables are gone. The index now holds only the boundary potentials and
the landmark potentials, and hop distances are priced by the landmark bound alone:
`d̂(u, v) = LB(u, v)`. The sections above that mention tables describe the earlier design.

- `GapIndex::set_landmarks(sources, config)` no longer takes a radius. `table_radius`,
  `table_bytes`, `table_lookup` and the `--table-radius` CLI flag are removed. The index file
  format is now version 2.
- `compute_row` prices a row from the potential-major landmark block only. No node→defect map is
  needed any more.
- `exact_distance` (used only for `gap_ub`) is always an A* with the landmark heuristic.
- `all_hops_exact` is now true only when the walk has no defect–defect hop.
- Tests: the table tests became `DHatIsALowerBound` and `EveryNodeALandmarkGivesExactSigma`
  (every node a landmark ⇒ LB = dist ⇒ index pricing equals exact pricing). The phenomenological
  tightness test now expects > 60% exact shots instead of > 90% (measured: 67%). Every soundness
  test is unchanged and passes: 30 C++ fastgap tests and 75 Python tests.
- `GapIndexCore.dual_state(events)` (new binding) returns one shot's final dual state: `y_leaf`,
  `rho`, `y_blossom`, ancestor blossoms and `mate`. It is used for the paper's figures.
- Benchmarks were re-run with `--no-upper-bound` (the paper does not report `gap_ub`), with
  `rmx_d15_p1e-3` added. All 7 DEMs: 0 soundness violations. Gap-stage medians at T = 1 fell
  (d = 17: 18.9 → 10.1 µs) because there is no table scan and no `gap_ub`. Since rows are now
  cheaper, T = 4 helps only at d = 17 (10.1 → 8.5 µs); at d = 7 and d = 13 it is slower than
  T = 1. Both numbers are reported as measured.

Decode that both methods share, after the re-run (median / p99):

| DEM | mean defects | fastgap decode (PyMatching) | baseline decode in the class of M* |
|---|---|---|---|
| rmx_d5_p5e-3 | 8.4 | 0.83 / 5.0 µs | 0.88 / 4.7 µs |
| rmx_d7_p5e-3 | 25.3 | 3.21 / 12.8 µs | 3.33 / 12.0 µs |
| rmx_d5_p1e-3 | 1.8 | 0.08 / 1.2 µs | 0.12 / 1.2 µs |
| rmx_d9_p1e-3 | 12.0 | 0.96 / 3.8 µs | 1.04 / 4.0 µs |
| rmx_d13_p1e-3 | 38.3 | 3.17 / 7.8 µs | 3.50 / 8.7 µs |
| rmx_d15_p1e-3 | 59.6 | 5.12 / 10.8 µs | 5.50 / 11.8 µs |
| rmx_d17_p1e-3 | 87.8 | 7.79 / 14.4 µs | 8.42 / 16.3 µs |

## Phase 1: dual state (§5.1)

- `extract_dual_state` validates and indexes the events first, then calls the upstream
  `process_timeline_until_completion`. Before shattering, it walks every defect's leaf →
  `blossom_parent` chain. It checks that each region is frozen, that
  `get_distance_at_time(now) == y_intercept()`, and that `region_that_arrived_top` agrees with
  the chain. Each defect stores its ancestors top-first with prefix sums, so `β(u,v)` is a
  common-prefix walk; when the top blossoms differ it returns 0 immediately.
- Leaf identification: the leaf is read from `nodes[d].region_that_arrived` *after* flooding,
  without hooking `create_detection_event`. The code checks on every shot that this region has
  no children and that `reached_from_source == &node`. The checks are cheap and have not fired on
  any of the several million shots decoded by the test suite: the source node of a detection
  event is never released by its own leaf region.
- `W* = Σ y` is read before shattering. Tests compare it bit-for-bit with upstream's weight on
  random grids (blossom-heavy) and on surface codes.
- Upstream corner cases: events on user-declared boundary nodes are skipped exactly as upstream
  does. Negative weights are rejected (p > 0.5 at canonicalisation; any negative-weight state in
  the Mwpm). "No perfect matching" propagates, and the Mwpm stays usable afterwards (tested). Bad
  shots (duplicate or out-of-range events) are rejected before flooding, so the Mwpm is never
  left half-processed.

## Phase 2: references (§5.2)

- `Pricing::EXACT` runs the same search with exact `dist`, computing one Dijkstra row per settled
  state. This gives `Δ_σ`.
- `ExactBaseline` runs two upstream decodes on the baseline graph, with the observable detector
  off and on. The baseline graph is built from the same canonical spec. Its normalising constant
  is asserted equal to the decode graph's, both when the index is built and when the baseline is
  constructed.
- `brute_force_gap`: a subset DP over ≤ 20 defects in the L/R model, including the bare L–R
  chain. On random small DEMs it equals the two-decode baseline in both `W_even` and `W_odd`.

## Phase 3: index and priced search (§5.3)

- Canonicalisation works on `dem.flattened()`, splitting components at separators, as
  PyMatching's own iterator does. The canonical DEM keeps every instruction in the same order and
  changes only the observable-k labels of graphlike components, so PyMatching builds a
  structurally identical graph. Merges and the first-edge observables are the same, and so are
  `M*` and the tie-breaking. Upstream PyMatching itself parses `str(dem)`, which is not an exact
  round trip of the probabilities, so fastgap also builds from `str(dem)`. The index stores that
  text and re-canonicalises it on load.
- The distance graph is read from the decoder's own `MatchingGraph`, so its edges, weights and
  units match the decoder's by construction. A 64-bit fingerprint of the graph is stored in the
  index file and checked on load.
- Potentials: exact int64 in the index and the file. Pricing uses a saturated int32 copy (cap
  2³¹−1, §2.1), which halves the work and lets the bound vectorise; saturation keeps the bound
  valid.
- Tables: a truncated Dijkstra from every detector, stored as CSR sorted by neighbour id. The
  default radius is 4 × the median edge weight (§7.4).
- Search: dense Dijkstra over `|D|+2` states, as in §3.3. Implementation details that matter for
  speed (prices are identical to the plain per-pair formula; tests compare against it):
  - Labels live per *hop target*, not per state. Every target `u` (a defect not matched to `L`)
    leads by its relay to exactly one place, `mate(u)` or `R`, so the label of state `mate(u)`
    can sit at `u`'s position. Relaxing a row is then one contiguous, branch-free `min` pass, and
    the argmin for the next iteration is fused into it. Targets that lead to states are sorted by
    state index, so the first minimum is the serial tie-break. The start states (`L`, and defects
    matched to `L`) have label 0, are never relaxed, and are kept in a separate sorted list.
  - Each settled row is priced in one pass. The landmark bound is computed over a
    potential-major block of the targets. The row's exact table is then scanned once through a
    node→defect map, instead of one binary search per pair.
  - Ties between finished walks are broken by (iteration, defect index), so the walk, and with it
    `gap_ub`, never depends on the thread count.
- **Intra-shot parallelism.** `num_threads = T` threads cooperate on *each* shot; shots are always
  decoded one after another (`decode_batch` is a plain loop). The aim is lower per-shot latency,
  not throughput.
  - Thread 0 runs the Dijkstra. The other `T − 1` threads *prefetch price rows*: each one
    repeatedly claims (with a CAS) the free state whose tentative label is smallest and prices
    its row into a shared row cache. Those are the states thread 0 settles next. Thread 0
    publishes a label to the prefetchers only when it drops.
  - When thread 0 settles `x`, it uses row `x` if it is ready. Otherwise it prices the row itself,
    or, if a prefetcher holds the claim, waits for it. A row's content does not depend on who
    priced it, so every thread count gives bit-identical results (`gap_lb`, `gap_ub`, walk,
    censoring). C++ and Python tests check this at T ∈ {2, 3, 4, 7}.
  - Shared rows are int32, saturated at `cap = d_LR + 1`. This is exact, not approximate. After
    `L` settles, `best ≤ d_LR`, so no state with a label ≥ `cap` is ever settled or on the walk,
    and no candidate ≥ `cap` is ever the answer. Before `L` settles, only label-0 states settle.
    When `d_LR + 1` does not fit in int32, the shot runs single-threaded.
  - The threads are a `SpinTeam`. A condvar wake-up costs several µs, more than a whole search, so
    idle workers spin for 2 ms after their last job before they sleep. On macOS they request the
    user-interactive QoS class (performance cores).
  - **Rejected: one barrier per settled state.** In that design each thread owned a slice of the
    targets and all threads joined an all-reduce per iteration. On this M5 Pro an all-reduce costs
    70–150 ns at T = 2–4, about half of one row's work, so d=17 only reached 1.3×.
  - **Rejected: gathering the potentials in parallel** at the start of a shot. The extra barrier
    costs more than the ~0.5 µs it saves.
  - Scaling is bounded by the core cluster. On this machine the 5 "Super" cores share one L2.
    T ≤ 5 stays inside that cluster; at T > 5, prefetched rows cross clusters, and reading them
    costs thread 0 more than pricing them itself. In microbenchmarks T = 4 and T = 5 were
    within ~5% of each other.
- Speed history at d=17, p=1e-3 (88 defects per shot on average; median gap-stage time per shot,
  no upper bound). Rows above the line were measured on the earlier machine, rows below it on
  the M5 Pro:

  | step | time |
  |---|---|
  | first version: per-pair binary search plus int64 LB | 58.5 µs |
  | row scan of the table, gathered potentials | 37.6 µs |
  | int32 saturated potentials | 31.8 µs |
  | active-target compaction | 27.9 µs |
  | — M5 Pro: same code — | 21.0 µs |
  | labels per target position, branch-free relax, fused argmin (T = 1) | 15.5 µs |
  | row prefetching, T = 4 | 9.3 µs |

  For comparison, the decode is about 9.5 µs (earlier machine) and the exact baseline about
  600 µs.
- **Rejected idea, recorded so nobody retries it:** A* with the potential `h_R − ρ`. The intended
  consistency argument needs an *upper* bound on `h_R(u)` across the relay `u ⇒ m`, but tightness
  only gives `h_R(u) ≤ h_R(m) + ρ_u + ρ_m`, which leaves a slack of 2(ρ_u + ρ_m). Relays are free
  in reduced-cost space, so there is no consistent goal potential of this form. The §9 soundness
  test (`gap_lb ≤ gap_exact`) failed immediately on surface codes, and the change was reverted.

## Phase 4: upper bound and diagnostics (§5.4)

- Exact hop distances use an A* on the interior graph with the landmark bound as heuristic.
  The heuristic is consistent because a max of 1-Lipschitz potentials is 1-Lipschitz, saturated
  or not. It is verified against brute-force Dijkstra in `interior_graph.test.cc`. Table hits are
  reused, and L/R hops use `h_L`/`h_R`. These are the only targeted distance queries.
- Diagnostics per shot: `walk_simple`, `all_hops_exact`, `num_hops`, `hop_slack = Σ(dist − d̂)`
  (only with `upper_bound`), `num_blossoms`, and `blossom_pass` (a run inside a blossom entered
  and left by hops, i.e. `n_S > 0`).
- Debug invariants (`check=True`; always on in the C++ tests through `FASTGAP_CHECK`):
  - boundary reduced costs are ≥ 0;
  - relays to L and R are tight;
  - for every defect-defect relay, `ρ_u + ρ_v − 2β == dist` (brute-force Dijkstra);
  - all-pairs `σ ≥ 0` (all pairs when there are ≤ 64 defects; for larger shots, all pairs among
    the first 64);
  - `d̂ ≤ dist`, with equality on table hits;
  - `gap_ub ≥ gap_lb`.

## Phase 5: bindings, CLI, benchmark (§5.5)

- `python -m fastgap bench | build-index | landmark-report | serve`. The bench CSV has exactly
  the columns of §5.5 plus a trailing `strategy_hash` column, so every row carries the hash (§7.1).
- Timing:
  - compute-only per-stage timers in C++;
  - 1000 warm-up shots;
  - median, p99 and p99.9;
  - per-shot latency at T = 1 and at every T in `--threads` (threads per shot; shots are decoded
    one after another, so no throughput or core-seconds figures are reported);
  - end-to-end latency (`--e2e`): a `python -m fastgap serve` subprocess reads bit-packed shots
    from a pipe and writes answers back. The bench measures the round trip per shot, so this
    number includes process-boundary I/O and the Python wrapper, as Craig Gidney asked for.

## Results

### Soundness (§9.3, §9.4, §9.5): 10⁵ shots per configuration

Every configuration below passed all of these checks:
- `gap_lb ≤ gap_exact` in internal integer units on every shot;
- `gap_ub ≥ gap_exact` on every simple walk;
- predictions and weights bit-identical to `pymatching.Matching.decode_batch` on the original DEM;
- the in-decoder invariants on every shot for d ≤ 7, and on the first 2000 shots for d ≥ 9.

The raw per-configuration stats are in `benchmarks/fastgap/results/soundness_1e5.jsonl`. Strategy
`auto` (= `faces_corners` here), K ≤ 16, default table radius. Times are the per-shot
compute-only stage timers on the M5 Pro. Shots were decoded one after another with T = 4 threads
per shot and `parallel_grain = 1`, so even tiny shots took the parallel path; that costs ~0.2 µs
on d = 3 shots, which the default grain avoids. The test also checks every shot bit-for-bit
against a single-threaded decoder. The gap time includes `gap_ub`.

| basis | d | p | blossom shots | blossom passes | walk simple | gap_lb == exact | t_decode (median) | t_gap (median / p99) | t_exact (median) |
|---|---|---|---|---|---|---|---|---|---|
| X | 3 | 0.005 | 1.4% | 0.01% | 100.00% | 100.0% | 0.08 µs | 0.38 / 1.1 µs | 0.9 µs |
| Z | 3 | 0.005 | 0.9% | 0.02% | 100.00% | 100.0% | 0.08 µs | 0.33 / 1.1 µs | 0.9 µs |
| X | 5 | 0.005 | 9.7% | 0.22% | 99.81% | 99.6% | 0.88 µs | 1.17 / 2.6 µs | 8.8 µs |
| Z | 5 | 0.005 | 8.2% | 0.22% | 99.88% | 99.6% | 0.83 µs | 1.12 / 2.5 µs | 8.0 µs |
| X | 7 | 0.005 | 29.1% | 0.93% | 98.95% | 95.4% | 3.42 µs | 2.67 / 5.7 µs | 29.7 µs |
| Z | 7 | 0.005 | 26.4% | 0.88% | 99.19% | 95.1% | 3.38 µs | 2.50 / 5.3 µs | 27.9 µs |
| X | 3 | 0.01 | 3.6% | 0.07% | 99.98% | 99.9% | 0.17 µs | 0.50 / 1.4 µs | 1.2 µs |
| Z | 3 | 0.01 | 3.0% | 0.04% | 100.00% | 100.0% | 0.17 µs | 0.46 / 1.3 µs | 1.1 µs |
| X | 5 | 0.01 | 24.5% | 0.65% | 99.12% | 98.7% | 2.50 µs | 1.75 / 3.8 µs | 12.0 µs |
| Z | 5 | 0.01 | 23.3% | 0.77% | 99.38% | 98.9% | 2.42 µs | 1.67 / 3.6 µs | 11.3 µs |
| X | 7 | 0.01 | 66.3% | 3.01% | 95.77% | 94.0% | 10.08 µs | 4.58 / 9.6 µs | 39.2 µs |
| Z | 7 | 0.01 | 64.9% | 2.96% | 96.76% | 94.9% | 9.75 µs | 4.12 / 8.4 µs | 37.4 µs |
| X | 9 | 0.0005 | 0.7% | 0.02% | 99.99% | 88.9% | 0.46 µs | 1.00 / 5.0 µs | 33.6 µs |
| Z | 9 | 0.0005 | 0.7% | 0.01% | 100.00% | 91.3% | 0.46 µs | 0.96 / 4.8 µs | 30.3 µs |
| X | 13 | 0.0005 | 1.8% | 0.04% | 99.97% | 45.3% | 1.54 µs | 3.79 / 13.1 µs | 133.8 µs |
| Z | 13 | 0.0005 | 1.9% | 0.06% | 99.97% | 53.5% | 1.58 µs | 3.50 / 12.0 µs | 123.2 µs |
| X | 17 | 0.0005 | 3.7% | 0.08% | 99.92% | 17.8% | 3.71 µs | 7.88 / 22.7 µs | 397.0 µs |
| Z | 17 | 0.0005 | 3.9% | 0.10% | 99.94% | 25.5% | 3.71 µs | 7.29 / 19.5 µs | 362.8 µs |
| X | 9 | 0.001 | 2.7% | 0.07% | 99.95% | 78.5% | 1.00 µs | 1.75 / 6.0 µs | 42.0 µs |
| Z | 9 | 0.001 | 2.9% | 0.08% | 99.96% | 82.0% | 1.00 µs | 1.58 / 5.7 µs | 38.4 µs |
| X | 13 | 0.001 | 7.7% | 0.20% | 99.83% | 40.0% | 3.33 µs | 5.50 / 12.9 µs | 170.3 µs |
| Z | 13 | 0.001 | 8.0% | 0.20% | 99.87% | 46.0% | 3.33 µs | 5.12 / 11.6 µs | 156.2 µs |
| X | 17 | 0.001 | 15.7% | 0.39% | 99.65% | 24.4% | 8.00 µs | 11.96 / 22.0 µs | 556.8 µs |
| Z | 17 | 0.001 | 16.7% | 0.37% | 99.72% | 31.0% | 7.96 µs | 10.92 / 19.5 µs | 502.1 µs |

Reading the table:
- Blossoms are common at d=7, p=1e-2 (65% of shots), but the bound almost never loses anything
  to them. The `Δ − Δ_σ` part of the shortfall is ~0 on every DEM measured.
- All of the looseness at large d comes from the landmark bound on long hops (`Δ_σ − Δ_σ,H`).

### Tightness vs index parameters (d=13, p=1e-3, rotated memory X, 5000 shots)

| strategy | K | table radius | gap_lb == exact | mean shortfall | p90 shortfall | gap stage (median, before speed-ups) | table memory |
|---|---|---|---|---|---|---|---|
| boundary_only | 0 | 4w | 16.2% | 16.9 dB | 37.7 dB | 12.2 µs | 3.9 MB |
| faces_corners | 12 | 4w | 40.0% | 5.9 dB | 18.1 dB | 13.0 µs | 3.9 MB |
| fps | 16 | 4w | 43.3% | 5.4 dB | 16.5 dB | 12.6 µs | 3.9 MB |
| faces_corners + fps | 32 | 4w | 47.8% | 4.8 dB | 14.8 dB | 14.0 µs | 3.9 MB |
| faces_corners + fps | 64 | 4w | 60.7% | 3.4 dB | 12.5 dB | 17.2 µs | 3.9 MB |
| greedy_cover | 16 | 4w | 58.2% | 2.8 dB | 10.1 dB | 14.2 µs | 3.9 MB |
| faces_corners + fps | 16 | **8w** | **92.7%** | **0.55 dB** | **0.0 dB** | 18.0 µs | 16.3 MB |
| faces_corners + fps | 16 | 0 | 12.8% | 10.4 dB | 25.5 dB | 7.3 µs | 0 |

(`w` = median edge weight. The landmark-report histogram puts the entire shortfall in the LB
part: `Δ_σ == Δ` on > 99.9% of shots.)

Take-aways for future experiments:
1. **The table radius is the main tightness lever at large distance.** The default of 4 edge
   weights (§7.4) is enough for d ≤ 7. For d ≥ 13, memory of about 16 MB at R = 8w buys exact
   bounds on more than 90% of shots.
2. **Among landmark-only options, `greedy_cover` is the best use of a fixed K.** Plain FPS and
   face/corner placements need about 4× more landmarks for the same tightness.
3. Post-selection decisions (`threshold`) are cheap regardless. Early termination stops at τ, and
   `classify(τ)` is never wrong, which is tested against the exact gap.

### Benchmarks (latency only)

Every DEM in `benchmarks/fastgap/index_config.json` was benchmarked with 2×10⁴ shots, 1000
warm-up shots, the `auto` strategy and the default table radius. Shots are sampled from the DEM
and decoded one after another. `T` is the number of threads that cooperate on each shot's gap
search, *including* the calling thread (thread 0 plus T − 1 prefetchers). T = 1 is the baseline
and T = 4 comes from the config. The raw per-shot CSVs (`t_decode_ns`, `t_gap_ns` for T = 1,
`t_gap_ns_t4`, `t_exact_ns`) and summaries are in `benchmarks/fastgap/results/`.

**Gap computation: fastgap vs the exact PyMatching baseline.** The baseline is the method of
Gidney's cultivation and yoking code: two PyMatching decodes of the canonical DEM, with the
observable as an extra detector toggled off and on. The two decodes are timed separately. The
one in the class of M* yields W* and the prediction, just like an ordinary decode. The one in
the complementary class yields W^c, and is what the baseline pays for the gap. fastgap also runs
one ordinary decode, then computes the gap from its dual state (`t_gap`, including `gap_ub`).
The decode is common to both methods and is reported separately, so the speed-ups below are for
the gap computation only:

    speed-up = baseline complementary-class decode latency / fastgap gap-stage latency

All values are compute-only per-shot stage timers (median / p99).

| DEM | mean defects | baseline gap (complementary decode) | fastgap gap, T=1 | fastgap gap, T=4 | speed-up T=1 (med / p99) | speed-up T=4 (med / p99) | violations |
|---|---|---|---|---|---|---|---|
| rmx_d5_p5e-3 | 8.4 | 7.17 / 14.6 µs | 0.75 / 2.5 µs | 0.75 / 2.7 µs | 9.6× / 5.7× | 9.6× / 5.5× | 0 |
| rmx_d7_p5e-3 | 25.3 | 24.67 / 44.5 µs | 2.75 / 7.1 µs | 2.62 / 5.6 µs | 9.0× / 6.3× | 9.4× / 7.9× | 0 |
| rmx_d5_p1e-3 | 1.8 | 3.67 / 10.0 µs | 0.12 / 0.9 µs | 0.12 / 0.9 µs | 29.3× / 10.9× | 29.3× / 10.9× | 0 |
| rmx_d9_p1e-3 | 12.0 | 38.50 / 62.7 µs | 1.50 / 5.9 µs | 1.54 / 5.7 µs | 25.7× / 10.7× | 25.0× / 11.0× | 0 |
| rmx_d13_p1e-3 | 38.3 | 158.71 / 232.8 µs | 6.71 / 14.8 µs | 5.50 / 13.0 µs | 23.7× / 15.8× | 28.9× / 17.9× | 0 |
| rmx_d17_p1e-3 | 87.8 | 535.88 / 808.5 µs | 18.88 / 34.2 µs | 11.92 / 21.9 µs | 28.4× / 23.6× | 45.0× / 37.0× | 0 |

The decode that both methods share (not part of the speed-ups):

| DEM | mean defects | fastgap decode (PyMatching) | baseline decode in the class of M* |
|---|---|---|---|
| rmx_d5_p5e-3 | 8.4 | 0.83 / 5.2 µs | 0.92 / 4.9 µs |
| rmx_d7_p5e-3 | 25.3 | 3.29 / 13.3 µs | 3.46 / 12.5 µs |
| rmx_d5_p1e-3 | 1.8 | 0.12 / 1.2 µs | 0.12 / 1.3 µs |
| rmx_d9_p1e-3 | 12.0 | 1.00 / 3.9 µs | 1.08 / 4.2 µs |
| rmx_d13_p1e-3 | 38.3 | 3.33 / 7.9 µs | 3.62 / 8.8 µs |
| rmx_d17_p1e-3 | 87.8 | 7.83 / 14.1 µs | 8.71 / 16.0 µs |

Other per-DEM numbers:

| DEM | detectors | fastgap gap speed-up T=4 vs T=1 (median) | end-to-end med, T=1 / T=4 (µs) | gap_lb == exact | LB/dist on hops (median) |
|---|---|---|---|---|---|
| rmx_d5_p5e-3 | 120 | 1.00× | 15.2 / 13.7 | 99.6% | 1.000 |
| rmx_d7_p5e-3 | 336 | 1.05× | 21.3 / 19.1 | 95.5% | 1.000 |
| rmx_d5_p1e-3 | 120 | 1.00× | 13.6 / 13.5 | 100.0% | 1.000 |
| rmx_d9_p1e-3 | 720 | 0.97× | 16.0 / 14.0 | 78.8% | 0.996 |
| rmx_d13_p1e-3 | 2184 | 1.22× | 24.7 / 20.0 | 40.0% | 0.996 |
| rmx_d17_p1e-3 | 4896 | 1.58× | 44.2 / 31.0 | 24.3% | 1.000 |

Notes:
- The baseline's decode in the class of M* runs on a graph with one extra node (the
  observable detector), so it is slightly slower than fastgap's ordinary decode.
- End-to-end is the round trip of one shot through the `serve` pipe protocol, with `serve
  --threads T`: syndrome in from another process, answer out. Roughly 10 µs of it is the fixed
  cost of the pipe and the Python wrapper.
- No throughput figures: shot-level parallelism is gone, and every thread works on the shot in
  hand.

Reading the tables:
- fastgap's gap computation is 9–29× faster than the baseline's complementary decode at
  T = 1, and up to 45× at T = 4 (d = 17). The complementary decode has to grow a region
  across the whole spacetime volume, while the gap search touches only the defects.
- Four threads help only on large shots. At d = 17 the gap stage drops from 18.9 to 11.9 µs, and
  the speed-up rises from 28× to 45×; at d = 13, from 24× to 29×. On shots with fewer than
  `parallel_grain` = 16 defects a single thread does the search, and at d = 9 T = 4 is ~3%
  slower.
- The single-threaded gap stage is faster than before this change (d = 17: 29.5 → 18.9 µs
  median, including `gap_ub`), thanks to the target-indexed, branch-free relaxation.

## Known gaps / not done

- **The real cultivation escape-stage DEM has not been tested.** The DEM comes from Gidney's
  `magic-state-cultivation` repository, which is not available in this environment. What is in
  place:
  - the `grafted` strategy, tested with virtual-node annotations supplied through `meta`;
  - `GapIndex.from_gap_dem`, which builds fastgap's graph from exactly the edges PyMatching keeps
    and asserts that edges, weights and `normalising_constant` match the reference graph;
  - tests of `from_gap_dem` on synthetic gap DEMs that contain 3-detector errors PyMatching drops.

  Run `from_gap_dem` on the real `gap_dem` before using fastgap for cultivation, and check which
  coordinate the virtual pair nodes carry. The `grafted` default assumes a 5th coordinate.
- **Yoked patches** are approximated by two independent patches whose logicals are XORed into one
  observable. It passes canonicalisation and soundness, but it is not Gidney's actual yoked
  circuit.
- **Bazel** is not installed here. The `BUILD` change is a one-line `copts` addition, and the
  globs already exclude `*.test.cc`, `*.test.h` and `*.pybind.cc` from the library target.
- **CPU pinning** is implemented for Linux only.
