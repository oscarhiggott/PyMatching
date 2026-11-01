"""Index round-trip (§9.6), landmark strategies (§7), the reference gap-DEM path (§7.3), CLI."""

import json
import subprocess
import sys

import numpy as np
import pytest

import fastgap
from fastgap import landmarks
from fastgap_helpers import phenomenological_dem, random_grid_dem, sample, surface_dem, with_obs_detector


def test_index_round_trip_bit_identical(tmp_path):
    dem = surface_dem("rotated_memory_z", 5, 0.01)
    idx = fastgap.GapIndex.from_dem(dem)
    path = tmp_path / "x.idx"
    idx.save(path)
    loaded = fastgap.GapIndex.load(path)
    assert loaded.strategy_hash == idx.strategy_hash
    assert loaded.strategy_config == idx.strategy_config
    assert [l.label for l in loaded.landmarks] == [l.label for l in idx.landmarks]
    dets, _ = sample(dem, 5000, 1)
    a = fastgap.GapDecoder(idx).decode_batch(dets, upper_bound=True)
    b = fastgap.GapDecoder(loaded).decode_batch(dets, upper_bound=True)
    for field in ("prediction", "weight_int", "gap_lb_int", "gap_ub_int", "censored", "walk_simple", "all_hops_exact"):
        assert np.array_equal(getattr(a, field), getattr(b, field)), field


def test_index_rejects_corrupt_file(tmp_path):
    dem = surface_dem("rotated_memory_z", 3, 0.01)
    idx = fastgap.GapIndex.from_dem(dem)
    path = tmp_path / "x.idx"
    idx.save(path)
    data = bytearray(path.read_bytes())
    data[-3] ^= 0xFF
    data[:8] = b"NOTFGAP!"
    path.write_bytes(bytes(data))
    with pytest.raises(RuntimeError):
        fastgap.GapIndex.load(path)


@pytest.mark.parametrize("strategy", ["boundary_only", "faces", "faces_corners", "fps", "greedy_cover", "phenomenological", "grafted", "auto"])
def test_strategies_are_sound(strategy):
    dem = surface_dem("rotated_memory_x", 5, 0.01)
    idx = fastgap.GapIndex.from_dem(dem, strategy=strategy, num_landmarks=8)
    assert len(idx.landmarks) <= 8
    assert idx.strategy_config["num_landmarks"] == 8
    g = idx.graph
    rng = np.random.default_rng(0)
    for u in rng.integers(0, g.num_nodes, size=10):
        row = g.dijkstra([u])
        for v in range(g.num_nodes):
            assert idx._core.lower_bound(int(u), v) <= row[v]
    dets, _ = sample(dem, 3000, 2)
    res = fastgap.GapDecoder(idx).decode_batch(dets)
    exact = fastgap.exact_gap_batch(idx, dets)
    assert np.all(res.gap_lb_int <= exact.gap_int)


def test_more_landmarks_never_hurt_much():
    dem = surface_dem("rotated_memory_x", 7, 0.01)
    dets, _ = sample(dem, 3000, 3)
    exact = fastgap.exact_gap_batch(fastgap.GapIndex.from_dem(dem, strategy="boundary_only"), dets)
    fracs = {}
    for s in ("boundary_only", "faces_corners"):
        idx = fastgap.GapIndex.from_dem(dem, strategy=s)
        res = fastgap.GapDecoder(idx).decode_batch(dets)
        fracs[s] = np.mean(res.gap_lb_int == exact.gap_int)
    assert fracs["faces_corners"] >= fracs["boundary_only"]


def test_auto_picks_by_dem():
    g = fastgap.GapIndex.from_dem(surface_dem("rotated_memory_x", 3, 0.01))
    assert g.strategy_name == "faces_corners"
    no_coords = "\n".join(l for l in str(surface_dem("rotated_memory_x", 3, 0.01)).splitlines() if not l.strip().startswith("detector"))
    g2 = fastgap.GapIndex.from_dem(no_coords)
    assert g2.strategy_name == "fps"


def test_custom_strategy_registration():
    class FirstNodes:
        name = "first_nodes"

        def params(self):
            return {"n": 2}

        def select(self, graph, coords, meta, k):
            return [landmarks.LandmarkSource(frozenset([0, 1]), "first")]

    landmarks.register("first_nodes", FirstNodes)
    idx = fastgap.GapIndex.from_dem(surface_dem("rotated_memory_x", 3, 0.01), strategy="first_nodes")
    assert idx.strategy_name == "first_nodes"
    assert idx.landmarks[0].nodes == frozenset([0, 1])
    assert isinstance(FirstNodes(), fastgap.LandmarkStrategy)


def test_grafted_on_synthetic_virtual_nodes():
    # Pretend a patch region is irregular: tag some detectors as virtual pair nodes via meta.
    dem = surface_dem("rotated_memory_x", 5, 0.01)
    idx = fastgap.GapIndex.from_dem(
        dem, strategy="grafted", num_landmarks=10,
        meta={"virtual_nodes": {"r": [10, 11, 12], "g": [20, 21]}, "switch_time": 2},
    )
    labels = [l.label for l in idx.landmarks]
    assert labels[:3] == ["virtual:all", "virtual:g", "virtual:r"]
    assert "slice:switch" in labels
    assert len(labels) == 10


def test_from_gap_dem_matches_reference():
    """§7.3: fastgap's graph equals the edges PyMatching keeps for a gap_dem, even when the gap
    DEM has 3-detector errors (observable on interior edges) that PyMatching drops."""
    rng = np.random.default_rng(8)
    raw = random_grid_dem(rng, 6, 5, gauge_density=0.1)  # observable on some interior edges too
    gap_dem = with_obs_detector(raw)
    idx = fastgap.GapIndex.from_gap_dem(gap_dem, strategy="fps", num_landmarks=4)
    assert idx.num_detectors == raw.num_detectors + 1
    assert any(len([t for t in inst.targets_copy() if t.is_relative_detector_id()]) == 3
               for inst in gap_dem.flattened() if inst.type == "error"), "want 3-detector errors in gap_dem"
    # Sample from the edges the reference graph keeps (errors PyMatching dropped cannot occur in
    # the reference decoder's model either).
    import pymatching

    m = pymatching.Matching.from_detector_error_model(gap_dem)
    pymatching.set_seed(4)
    full = np.array([m.add_noise()[1] for _ in range(2000)], dtype=np.uint8)
    full[:, -1] = 0
    res = fastgap.GapDecoder(idx, check=True).decode_batch(full, bit_packed_shots=False, upper_bound=True)
    exact = fastgap.exact_gap_batch(idx, full, bit_packed_shots=False)
    assert np.all(res.gap_lb_int <= exact.gap_int)
    assert np.all(~res.walk_simple | (res.gap_ub_int >= exact.gap_int))

    # Reference two-decode (Gidney's method) directly on gap_dem with pymatching. Dropping the
    # 3-detector errors can leave pieces whose only "boundary" is the observable detector; the
    # reference decoder cannot match those shots in one class, so they are skipped here.
    compared = 0
    for i, row in enumerate(full):
        on = row.copy()
        on[-1] = 1
        try:
            _, w_off = m.decode(row, return_weight=True)
            _, w_on = m.decode(on, return_weight=True)
        except ValueError:
            continue
        compared += 1
        assert abs(abs(w_on - w_off) - exact.gap[i]) < 1e-9
    assert compared > 0.5 * len(full)


def test_from_gap_dem_of_canonical_surface_code_equals_dem_path():
    dem = surface_dem("rotated_memory_x", 5, 0.01)
    idx = fastgap.GapIndex.from_dem(dem)
    gap_dem = with_obs_detector(idx.canonical_dem())
    ref = fastgap.GapIndex.from_gap_dem(gap_dem)
    assert ref.normalising_constant == idx.normalising_constant
    assert ref.d_lr_int == idx.d_lr_int
    dets, _, _ = dem.compile_sampler(seed=0).sample(2000)
    full = np.zeros((len(dets), ref.num_detectors), dtype=np.uint8)
    full[:, : dem.num_detectors] = dets
    a = fastgap.exact_gap_batch(idx, dets, bit_packed_shots=False)
    b = fastgap.exact_gap_batch(ref, full, bit_packed_shots=False)
    assert np.array_equal(a.gap_int, b.gap_int)


def test_phenomenological_faces_corners_tight():
    dem = phenomenological_dem(7, 0.03)
    idx = fastgap.GapIndex.from_dem(dem, strategy="phenomenological")
    dets, _ = sample(dem, 3000, 0)
    res = fastgap.GapDecoder(idx).decode_batch(dets)
    exact = fastgap.exact_gap_batch(idx, dets)
    assert np.all(res.gap_lb_int <= exact.gap_int)
    # Landmarks only (no local tables): about two thirds of the shots are exact here.
    assert np.mean(res.gap_lb_int == exact.gap_int) > 0.6


def test_cli_build_index_bench_and_report(tmp_path):
    idx_path = tmp_path / "s.idx"
    run = lambda *a: subprocess.run([sys.executable, "-m", "fastgap", *a], check=True, capture_output=True, text=True)
    out = run("build-index", "--gen", "surface_code:rotated_memory_x", "-d", "3", "-p", "0.01", "--out", str(idx_path))
    assert json.loads(out.stdout)["strategy"] == "faces_corners"
    csv_path = tmp_path / "b.csv"
    out = run("bench", "--gen", "surface_code:rotated_memory_x", "-d", "3", "-p", "0.01", "--idx", str(idx_path),
              "--shots", "500", "--warmup", "100", "--threads", "2", "3", "--seed", "1", "--out", str(csv_path),
              "--e2e", "--e2e-shots", "50")
    summary = json.loads(out.stdout)
    assert summary["soundness_violations"] == 0
    assert summary["threads"] == [1, 2, 3]
    assert set(summary["latency_ns_by_threads"]) == {"1", "2", "3"}
    assert set(summary["end_to_end_ns_by_threads"]) == {"1", "2", "3"}
    assert "throughput" not in summary
    header = csv_path.read_text().splitlines()[0].split(",")
    assert header[:12] == ["shot", "num_defects", "w_star", "gap_exact", "gap_lb", "gap_ub", "t_decode_ns",
                           "t_gap_ns", "t_exact_ns", "walk_simple", "all_hops_exact", "censored"]
    assert header[-4:] == ["t_gap_ns_t2", "t_gap_ns_t3", "t_exact_mstar_ns", "t_exact_comp_ns"]
    assert set(summary["gap_speedup_vs_exact_by_threads"]) == {"1", "2", "3"}
    out = run("landmark-report", "--gen", "surface_code:rotated_memory_x", "-d", "3", "-p", "0.01",
              "--idx", str(idx_path), "--shots", "300", "--hop-shots", "100", "--seed", "2")
    report = json.loads(out.stdout)
    assert report["soundness_violations"] == 0
    assert 0 <= report["fraction_gap_lb_equals_gap_exact"] <= 1
