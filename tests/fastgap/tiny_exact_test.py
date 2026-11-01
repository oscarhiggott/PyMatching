"""CLAUDE.md §9.1 (tiny exact) and §9.2 (gauge)."""

import numpy as np
import pytest
import stim

import fastgap
from fastgap_helpers import random_grid_dem, toric_dem, yoked_dem, phenomenological_dem, surface_dem


def _shots(dem, n, seed):
    dets, _, _ = dem.compile_sampler(seed=seed).sample(n)
    return dets


@pytest.mark.parametrize("seed", range(12))
@pytest.mark.parametrize("num_landmarks", [0, 3])
def test_tiny_exact_chain(seed, num_landmarks):
    """gap_lb <= Delta_sigma <= Delta == brute force == two-decode, gap_ub >= Delta when simple,
    with the in-decoder invariants (all-pairs sigma >= 0, relays tight, Sigma y = W*, d_hat <= dist)
    switched on."""
    rng = np.random.default_rng(seed)
    dem = random_grid_dem(rng, 3 + seed % 4, 3 + seed % 3)
    idx = fastgap.GapIndex.from_dem(dem, strategy="fps", num_landmarks=num_landmarks)
    dec = fastgap.GapDecoder(idx, check=True)
    dets = _shots(dem, 150, seed)
    dets = dets[dets.sum(axis=1) <= 10]
    res = dec.decode_batch(dets, upper_bound=True, bit_packed_shots=False)
    sigma = dec.decode_batch(dets, upper_bound=True, bit_packed_shots=False, pricing="exact")
    exact = fastgap.exact_gap_batch(idx, dets, bit_packed_shots=False)
    for i, row in enumerate(dets):
        brute = fastgap.brute_force_gap(idx, np.flatnonzero(row))
        assert brute["gap"] == exact.gap_int[i]
        assert brute["w_star"] == res.weight_int[i] == exact.w_star_int[i]
    assert np.all(res.gap_lb_int <= sigma.gap_lb_int)
    assert np.all(sigma.gap_lb_int <= exact.gap_int)
    assert np.all(~res.walk_simple | (res.gap_ub_int >= exact.gap_int))
    assert np.all(~sigma.walk_simple | (sigma.gap_ub_int >= exact.gap_int))
    # Sigma y = W*: the weight matches upstream PyMatching on the original DEM.
    import pymatching

    _, w = pymatching.Matching.from_detector_error_model(dem).decode_batch(dets, return_weights=True)
    assert np.array_equal(res.weight, w)


def test_hand_computed_repetition():
    dem = stim.DetectorErrorModel("""
        error(0.1) D0 L0
        error(0.1) D0 D1
        error(0.1) D1 D2
        error(0.1) D2
    """)
    idx = fastgap.GapIndex.from_dem(dem, strategy="boundary_only")
    dec = fastgap.GapDecoder(idx, check=True)
    w = np.log(0.9 / 0.1)
    r = dec.decode([1, 0, 0], upper_bound=True)
    assert r.prediction.tolist() == [1]
    assert r.weight == pytest.approx(w)
    assert r.gap_lb == pytest.approx(2 * w)
    assert r.gap_ub == pytest.approx(2 * w)
    assert r.gap_lb_db == pytest.approx(10 * np.log10(np.e) * 2 * w)
    r0 = dec.decode([0, 0, 0], upper_bound=True)
    assert r0.gap_lb == pytest.approx(idx.d_lr)
    assert idx.d_lr == pytest.approx(4 * w)


def test_gauge_shift_and_invariance():
    """Different gauges of the same code: identical gaps, predictions differing exactly by the
    per-shot shift; predictions in original labels always match upstream."""
    import pymatching

    rng = np.random.default_rng(3)
    dem = random_grid_dem(rng, 5, 4)
    idx = fastgap.GapIndex.from_dem(dem, strategy="fps", num_landmarks=4)
    gauge = idx.gauge
    assert np.any(gauge), "test DEM should need a non-trivial gauge"
    dets = _shots(dem, 2000, 1)
    res = fastgap.GapDecoder(idx).decode_batch(dets, bit_packed_shots=False)
    pred = pymatching.Matching.from_detector_error_model(dem).decode_batch(dets)
    assert np.array_equal(res.prediction, pred)

    # Canonical DEM: all-zero gauge, predictions shifted by XOR_{x in D} s(x).
    canon = idx.canonical_dem()
    idx2 = fastgap.GapIndex.from_dem(canon, strategy="fps", num_landmarks=4)
    assert not np.any(idx2.gauge)
    res2 = fastgap.GapDecoder(idx2).decode_batch(dets, bit_packed_shots=False)
    shift = (dets.astype(np.int64) @ gauge.astype(np.int64)) % 2
    assert np.array_equal(res2.prediction[:, 0], res.prediction[:, 0] ^ shift.astype(np.uint8))
    e1 = fastgap.exact_gap_batch(idx, dets, bit_packed_shots=False)
    e2 = fastgap.exact_gap_batch(idx2, dets, bit_packed_shots=False)
    assert np.array_equal(e1.gap_int, e2.gap_int)
    assert np.array_equal(res.gap_lb_int, res2.gap_lb_int)


def test_toroidal_dem_raises():
    with pytest.raises(fastgap.BoundarylessLogicalError):
        fastgap.GapIndex.from_dem(toric_dem(4))
    with pytest.raises(ValueError):
        fastgap.GapIndex.from_dem(toric_dem(5))


def test_parallel_edge_conflict_raises():
    dem = stim.DetectorErrorModel("""
        error(0.1) D0 L0
        error(0.1) D0
        error(0.1) D0 D1
        error(0.1) D1
    """)
    with pytest.raises(fastgap.ParallelEdgeConflictError):
        fastgap.GapIndex.from_dem(dem)


@pytest.mark.parametrize("name,dem", [
    ("surface_x", surface_dem("rotated_memory_x", 3, 0.01)),
    ("surface_z", surface_dem("rotated_memory_z", 3, 0.01)),
    ("yoked", yoked_dem(3, 0.01)),
    ("phenomenological", phenomenological_dem(5, 0.02)),
])
def test_experiment_dems_canonicalise(name, dem):
    """Every DEM family we plan to use passes §2.3 step 2, and the canonical DEM has the
    observable only on boundary edges."""
    idx = fastgap.GapIndex.from_dem(dem)
    canon = idx.canonical_dem()
    for inst in canon.flattened():
        if inst.type != "error":
            continue
        dets, obs = 0, False
        for t in inst.targets_copy() + [stim.target_separator()]:
            if t.is_separator():
                assert not (dets == 2 and obs)
                dets, obs = 0, False
            elif t.is_relative_detector_id():
                dets += 1
            elif t.is_logical_observable_id() and t.val == 0:
                obs = not obs
    dets = _shots(dem, 500, 0)
    res = fastgap.GapDecoder(idx, check=True).decode_batch(dets, bit_packed_shots=False, upper_bound=True)
    exact = fastgap.exact_gap_batch(idx, dets, bit_packed_shots=False)
    assert np.all(res.gap_lb_int <= exact.gap_int)


def test_threshold_censoring_user_units():
    dem = surface_dem("rotated_memory_x", 5, 0.01)
    idx = fastgap.GapIndex.from_dem(dem)
    dec = fastgap.GapDecoder(idx)
    dets = _shots(dem, 2000, 5)
    full = dec.decode_batch(dets, bit_packed_shots=False)
    for tau in (0.5, 2.0, 5.0, 20.0):
        r = dec.decode_batch(dets, threshold=tau, bit_packed_shots=False)
        assert np.array_equal(r.gap_lb >= tau, full.gap_lb >= tau)
        assert np.all(r.gap_lb[r.censored] >= tau)
        assert np.all(r.gap_lb_int[~r.censored] == full.gap_lb_int[~r.censored])
        cls = dec.decode_batch(dets, threshold=tau, upper_bound=True, bit_packed_shots=False).classify(tau)
        exact = fastgap.exact_gap_batch(idx, dets, bit_packed_shots=False)
        assert np.all(exact.gap[cls == 1] >= tau)
        assert np.all(exact.gap[cls == -1] < tau)
