# Copyright 2022 PyMatching Contributors

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#      http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os
from pathlib import Path

import numpy as np
from scipy.sparse import csc_matrix
import pytest
import networkx as nx

import pymatching
from pymatching import Matching


def repetition_code(n: int):
    row_ind, col_ind = zip(*((i, j) for i in range(n) for j in (i, (i + 1) % n)))
    data = np.ones(2 * n, dtype=np.uint8)
    return csc_matrix((data, (row_ind, col_ind)))


weight_fixtures = [10, 15, 20, 100]


@pytest.mark.parametrize("n", weight_fixtures)
def test_matching_weight(n: int):
    p = 0.4
    H = repetition_code(n)
    noise = np.random.rand(n) < p
    weights = np.random.rand(n)
    s = H @ noise % 2
    m = Matching(H, spacelike_weights=weights)
    corr, weight = m.decode(s, return_weight=True)
    expected_weight = np.sum(weights[corr == 1])
    assert expected_weight == pytest.approx(weight, rel=0.001)


def test_negative_weight_repetition_code():
    m = Matching()
    m.add_edge(0, 1, 0, -1)
    m.add_edge(1, 2, 1, -1)
    m.add_edge(2, 3, 2, -1)
    m.add_edge(3, 4, 3, -1)
    m.add_edge(4, 5, 4, -1)
    m.add_edge(5, 0, 5, -1)
    c, w = m.decode([0, 1, 1, 0, 0, 0], return_weight=True)
    assert np.array_equal(c, np.array([1, 0, 1, 1, 1, 1]))
    assert w == -5


def test_isolated_negative_weight():
    m = Matching()
    m.add_edge(0, 1, 0, 1)
    m.add_edge(1, 2, 1, -10)
    m.add_edge(2, 3, 2, 1)
    m.add_edge(3, 0, 3, 1)
    c, w = m.decode([0, 1, 1, 0], return_weight=True)
    assert np.array_equal(c, np.array([0, 1, 0, 0]))
    assert w == -10


def test_negative_and_positive_in_matching():
    g = nx.Graph()
    g.add_edge(0, 1, fault_ids=0, weight=1)
    g.add_edge(1, 2, fault_ids=1, weight=-10)
    g.add_edge(2, 3, fault_ids=2, weight=1)
    g.add_edge(3, 0, fault_ids=3, weight=1)
    m = Matching(g)
    c, w = m.decode([0, 1, 0, 1], return_weight=True)
    assert np.array_equal(c, np.array([0, 1, 1, 0]))
    assert w == pytest.approx(-9, rel=0.001)


def test_decode_to_matched_detection_events():
    num_nodes = 20
    m = Matching()
    m.add_boundary_edge(0)
    for i in range(num_nodes):
        m.add_edge(i, i + 1)
    m.add_boundary_edge(num_nodes)

    dets = np.array([2, 10, 12, 18])
    syndrome = np.zeros(m.num_detectors, dtype=np.uint8)
    syndrome[dets] = 1

    arr = m.decode_to_matched_dets_array(syndrome)
    assert np.array_equal(arr, np.array([[2, -1], [10, 12], [18, -1]]))

    d = m.decode_to_matched_dets_dict(syndrome)
    assert d == {2: None, 10: 12, 12: 10, 18: None}


def test_decode_to_matched_detection_events_with_negative_weights_raises_value_error():
    m = Matching()
    m.add_edge(0, 1, weight=-1)
    with pytest.raises(ValueError):
        m.decode_to_matched_dets_array([0, 0])

    with pytest.raises(ValueError):
        m.decode_to_matched_dets_dict([0, 0])


def test_matching_solution_integral_weights():
    m = Matching()
    m.add_boundary_edge(0, weight=3)
    m.add_edge(0, 1, weight=-1)
    m.add_edge(1, 2, weight=100)
    m.add_edge(2, 3, weight=12)
    corr, tot_weight = m.decode([0, 0, 0, 1], return_weight=True)
    assert tot_weight == 114
    m.add_edge(3, 4, weight=16777215)
    corr, tot_weight = m.decode([0, 0, 0, 0, 1], return_weight=True)
    assert tot_weight == 114 + 16777215
    m.add_edge(4, 5, weight=-16777215)
    corr, tot_weight = m.decode([0, 0, 0, 0, 0, 1], return_weight=True)
    assert tot_weight == 114


def get_full_data_path(filename: str) -> str:
    for data_dir in ("./PyMatching/data/", "./data/", "../data/", "../../data/"):
        fullpath = os.path.join(data_dir, filename)
        if os.path.isfile(fullpath):
            return fullpath
    raise ValueError(f"No data directory found inside {os.getcwd()}")


def test_surface_code_solution_weights(data_dir: Path):
    stim = pytest.importorskip("stim")
    dem = stim.DetectorErrorModel.from_file(
        data_dir / "surface_code_rotated_memory_x_13_0.01.dem"
    )
    m = Matching.from_detector_error_model(dem)
    shots = stim.read_shot_data_file(
        path=data_dir / "surface_code_rotated_memory_x_13_0.01_1000_shots.b8",
        format="b8",
        num_detectors=m.num_detectors,
        num_observables=m.num_fault_ids,
    )
    with open(
        data_dir
        / "surface_code_rotated_memory_x_13_0.01_1000_shots_no_buckets_weights_pymatchingv0.7_exact.txt",
        "r",
        encoding="utf-8",
    ) as f:
        expected_weights = [float(w) for w in f.readlines()]
    with open(
        data_dir
        / "surface_code_rotated_memory_x_13_0.01_1000_shots_no_buckets_predictions_pymatchingv0.7_exact.txt",
        "r",
        encoding="utf-8",
    ) as f:
        expected_observables = [int(w) for w in f.readlines()]
    assert shots.shape == (1000, m.num_detectors + m.num_fault_ids)
    weights = []
    predicted_observables = []
    for i in range(min(shots.shape[0], 1000)):
        prediction, weight = m.decode(
            shots[i, 0: -m.num_fault_ids], return_weight=True
        )
        weights.append(weight)
        predicted_observables.append(prediction)
    for weight, expected_weight in zip(weights, expected_weights):
        assert weight == pytest.approx(expected_weight, rel=1e-8)
    assert predicted_observables == expected_observables[0: len(predicted_observables)]

    expected_observables_arr = np.zeros((shots.shape[0], 1), dtype=np.uint8)
    expected_observables_arr[:, 0] = np.array(expected_observables)

    sampler = dem.compile_sampler()
    temp_shots, _, _ = sampler.sample(shots=10, bit_packed=True)
    assert temp_shots.shape[1] == np.ceil(dem.num_detectors // 8)

    batch_predictions = m.decode_batch(shots[:, 0: -m.num_fault_ids])
    assert np.array_equal(batch_predictions, expected_observables_arr)

    batch_predictions, batch_weights = m.decode_batch(
        shots[:, 0: -m.num_fault_ids], return_weights=True
    )
    assert np.array_equal(batch_predictions, expected_observables_arr)
    assert np.allclose(batch_weights, expected_weights, rtol=1e-8)

    bitpacked_shots = np.packbits(
        shots[:, 0: dem.num_detectors], bitorder="little", axis=1
    )
    batch_predictions_from_bitpacked, bitpacked_batch_weights = m.decode_batch(
        bitpacked_shots, return_weights=True, bit_packed_shots=True
    )
    assert np.array_equal(batch_predictions_from_bitpacked, expected_observables_arr)
    assert np.allclose(bitpacked_batch_weights, expected_weights, rtol=1e-8)

    bitpacked_batch_predictions_from_bitpacked, bitpacked_batch_weights = (
        m.decode_batch(
            bitpacked_shots,
            return_weights=True,
            bit_packed_shots=True,
            bit_packed_predictions=True,
        )
    )
    assert np.array_equal(
        bitpacked_batch_predictions_from_bitpacked, expected_observables_arr
    )
    assert np.allclose(bitpacked_batch_weights, expected_weights, rtol=1e-8)


def test_surface_code_solution_weights_with_correlations(data_dir: Path):
    stim = pytest.importorskip("stim")
    dem = stim.DetectorErrorModel.from_file(
        data_dir / "surface_code_rotated_memory_x_13_0.01.dem"
    )
    m = Matching.from_detector_error_model(dem)
    shots = stim.read_shot_data_file(
        path=data_dir / "surface_code_rotated_memory_x_13_0.01_1000_shots.b8",
        format="b8",
        num_detectors=m.num_detectors,
        num_observables=m.num_fault_ids,
    )
    # Test correlated decoding
    corr_weights_path = data_dir / "surface_code_rotated_memory_x_13_0.01_1000_shots_no_buckets_weights_pymatching_correlated.txt"
    corr_predictions_path = data_dir / "surface_code_rotated_memory_x_13_0.01_1000_shots_no_buckets_predictions_pymatching_correlated.txt"
    with open(
        corr_weights_path,
        "r",
        encoding="utf-8",
    ) as f:
        expected_corr_weights = [float(w) for w in f.readlines()]
    with open(
        corr_predictions_path,
        "r",
        encoding="utf-8",
    ) as f:
        expected_corr_observables = [int(w) for w in f.readlines()]

    expected_corr_observables_arr = np.zeros((shots.shape[0], 1), dtype=np.uint8)
    expected_corr_observables_arr[:, 0] = np.array(expected_corr_observables)

    m_corr = Matching.from_detector_error_model(dem, enable_correlations=True)
    corr_predictions, corr_weights = m_corr.decode_batch(
        shots[:, 0: -m.num_fault_ids],
        return_weights=True,
        enable_correlations=True,
    )
    # # Uncomment to update data files
    # with open(
    #     corr_weights_path,
    #     "w",
    #     encoding="utf-8",
    # ) as f:
    #     for w in corr_weights:
    #         print(w, file=f)
    # with open(
    #     corr_predictions_path,
    #     "w",
    #     encoding="utf-8",
    # ) as f:
    #     for pred in corr_predictions:
    #         print(pred[0], file=f)
    assert np.array_equal(corr_predictions, expected_corr_observables_arr)
    assert np.allclose(corr_weights, expected_corr_weights, rtol=1e-8)


def test_decode_batch_to_bitpacked_predictions():
    m = pymatching.Matching()
    m.add_edge(0, 1, fault_ids={0})
    m.add_edge(1, 2, fault_ids={10})
    m.add_edge(2, 3, fault_ids={3, 5})
    m.add_edge(3, 4, fault_ids={20, 16})

    predictions = m.decode_batch(
        np.array([[1, 0, 1, 0, 0], [0, 0, 1, 0, 1]], dtype=np.uint8),
        bit_packed_predictions=True,
    )
    assert np.array_equal(
        predictions, np.array([[1, 4, 0], [40, 0, 17]], dtype=np.uint8)
    )

    predictions = m.decode_batch(
        np.array([[5], [20]], dtype=np.uint8),
        bit_packed_shots=True,
        bit_packed_predictions=True,
    )
    assert np.array_equal(
        predictions, np.array([[1, 4, 0], [40, 0, 17]], dtype=np.uint8)
    )
    with pytest.raises(ValueError):
        m.decode_batch(np.array([[]], dtype=np.uint8))


def test_detection_event_too_large_raises_value_error():
    m = pymatching.Matching()
    m.add_edge(0, 1)
    with pytest.raises(ValueError):
        m.decode([1, 0, 1])
    with pytest.raises(ValueError):
        m.decode_batch([[8]], bit_packed_shots=True)


def test_decode_self_loops():
    m = Matching()
    m.add_boundary_edge(0, fault_ids={0}, weight=1)
    m.add_edge(0, 1, fault_ids={1}, weight=3)
    m.add_edge(1, 2, fault_ids={2}, weight=3)
    m.add_edge(2, 2, fault_ids={3}, weight=-100)
    m.add_edge(3, 3, fault_ids={4}, weight=-200)
    m.add_edge(4, 4, fault_ids={5}, weight=4)
    corr, weight = m.decode([0, 0, 1, 0, 0], return_weight=True)
    assert np.array_equal(corr, np.array([1, 1, 1, 1, 1, 0], dtype=np.uint8))
    assert weight == -293
    with pytest.raises(ValueError):
        m.decode([0, 0, 1, 1, 0])


def test_decode_wrong_syndrome_type_raises_type_error():
    with pytest.raises(ValueError):
        m = Matching()
        m.add_edge(0, 1)
        m.decode([0, "A"])


def test_syndrome_on_boundary_nodes():
    m = Matching()
    m.add_edge(0, 1, fault_ids={0})
    m.add_edge(1, 2, fault_ids={1})
    m.add_edge(2, 3, fault_ids={2})
    m.add_edge(3, 4, fault_ids={3})
    m.set_boundary_nodes({3, 4})
    m.decode([0, 0, 0, 1, 0])
    m.decode([0, 0, 0, 0, 1])
    m.decode([0, 0, 0, 1, 1])
    m.decode([1, 0, 1, 0, 1])


def test_decode_to_edges():
    m = Matching()
    m.add_boundary_edge(0)
    for i in range(10):
        m.add_edge(i, i + 1)
    edges = m.decode_to_edges_array([0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0])
    assert np.array_equal(
        edges,
        np.array([[9, 8], [5, 6], [4, 3], [5, 4], [0, 1], [0, -1]], dtype=np.int64),
    )
    edges = m.decode_to_edges_array(
        [False, True, False, True, False, False, True, False, True, True, False]
    )
    assert np.array_equal(
        edges,
        np.array([[9, 8], [5, 6], [4, 3], [5, 4], [0, 1], [0, -1]], dtype=np.int64),
    )


def test_parallel_boundary_edges_decoding():
    m = Matching()
    m.set_boundary_nodes({0, 2})
    m.add_edge(0, 1, fault_ids=0, weight=3.5)
    m.add_edge(1, 2, fault_ids=1, weight=2.5)
    assert np.array_equal(m.decode([0, 1]), np.array([0, 1], dtype=np.uint8))
    m.add_boundary_edge(1, fault_ids=100, weight=100)
    # Test pm::SearchGraph
    assert np.array_equal(np.nonzero(m.decode([0, 1]))[0], np.array([1], dtype=int))

    m = Matching()
    m.add_edge(0, 1, fault_ids=0, weight=-1)
    m.add_edge(0, 2, fault_ids=1, weight=3)
    m.add_boundary_edge(0, fault_ids=2, weight=-0.5)
    m.add_edge(0, 3, fault_ids=3, weight=-3)
    m.add_edge(0, 4, fault_ids=4, weight=-2)
    assert np.array_equal(
        m.decode([1, 0, 0, 0, 0]), np.array([0, 0, 1, 0, 0], dtype=np.uint8)
    )
    m.set_boundary_nodes({1, 2, 3, 4})
    assert np.array_equal(
        m.decode([1, 0, 0, 0, 0]), np.array([0, 0, 0, 1, 0], dtype=np.uint8)
    )


def test_decode_to_edges_with_correlations():
    stim = pytest.importorskip("stim")
    dem = stim.DetectorErrorModel("""
        error(0.1) D0 D1
        error(0.1) D2 D3
        error(0.3) D0 D1 ^ D2 D3
    """)
    m = Matching.from_detector_error_model(dem, enable_correlations=True)
    syndrome = np.array([1, 1, 1, 1])
    edges = m.decode_to_edges_array(syndrome, enable_correlations=True)
    edges.sort(axis=1)
    edges = edges[np.lexsort((edges[:, 1], edges[:, 0]))]
    expected_edges = np.array([[0, 1], [2, 3]])
    assert np.array_equal(edges, expected_edges)


def test_correlated_matching_handles_single_detector_components():
    stim = pytest.importorskip("stim")
    p = 0.1
    circuit = stim.Circuit.generated(
        code_task="surface_code:rotated_memory_x",
        distance=5,
        rounds=5,
        before_round_data_depolarization=p,
    )
    circ_str = str(circuit).replace(
        f"DEPOLARIZE1({p})", f"PAULI_CHANNEL_1(0, {p}, 0)"
    )
    noisy_circuit = stim.Circuit(circ_str)
    dem = noisy_circuit.detector_error_model(
        decompose_errors=True, approximate_disjoint_errors=True
    )
    m = Matching.from_detector_error_model(dem, enable_correlations=True)
    assert m.num_detectors > 0


def test_load_from_circuit_with_correlations():
    stim = pytest.importorskip("stim")
    circuit = stim.Circuit.generated(
        code_task="surface_code:rotated_memory_x",
        distance=3,
        rounds=3,
        after_clifford_depolarization=0.001
    )
    shots = circuit.compile_detector_sampler().sample(shots=10)
    matching_1 = pymatching.Matching(circuit, enable_correlations=True)
    matching_2 = pymatching.Matching.from_stim_circuit(circuit=circuit, enable_correlations=True)
    for m in (matching_1, matching_2):
        predictions, weights = m.decode_batch(shots=shots, return_weights=True, enable_correlations=True)


def test_use_correlations_with_uncorrelated_dem_load_raises_value_error(tmp_path):
    stim = pytest.importorskip("stim")
    d = 3
    p = 0.001
    circuit = stim.Circuit.generated(
        code_task="surface_code:rotated_memory_x",
        distance=d,
        rounds=d,
        after_clifford_depolarization=p
    )
    dem = circuit.detector_error_model(decompose_errors=True)
    shots = circuit.compile_detector_sampler().sample(shots=10)
    matching_1 = pymatching.Matching(circuit, enable_correlations=False)
    matching_2 = pymatching.Matching.from_stim_circuit(circuit=circuit, enable_correlations=False)
    matching_3 = pymatching.Matching.from_detector_error_model(
        model=dem,
        enable_correlations=False
    )
    fn = f"surface_code_x_d{d}_r{d}_p{p}"
    stim_file = tmp_path / f"{fn}.stim"
    circuit.to_file(stim_file)
    matching_4 = pymatching.Matching.from_stim_circuit_file(stim_file, enable_correlations=False)
    dem_file = tmp_path / f"{fn}.dem"
    dem.to_file(dem_file)
    matching_5 = pymatching.Matching.from_detector_error_model_file(dem_file, enable_correlations=False)
    for m in (matching_1, matching_2, matching_3, matching_4, matching_5):
        with pytest.raises(ValueError):
            m.decode_batch(shots=shots, return_weights=True, enable_correlations=True)
        with pytest.raises(ValueError):
            m.decode_batch(shots=shots, return_weights=False, enable_correlations=True)
        with pytest.raises(ValueError):
            m.decode_to_edges_array(shots[0], enable_correlations=True)
        with pytest.raises(ValueError):
            m.decode(shots[0], enable_correlations=True)


def test_use_correlations_without_decompose_errors_raises_value_error(tmp_path):
    stim = pytest.importorskip("stim")
    d = 3
    p = 0.001
    circuit = stim.Circuit.generated(
        code_task="surface_code:rotated_memory_x",
        distance=d,
        rounds=d,
        after_clifford_depolarization=p
    )
    dem = circuit.detector_error_model(decompose_errors=False)
    dem_file = tmp_path / "surface_code.dem"
    dem.to_file(dem_file)
    with pytest.raises(ValueError):
        pymatching.Matching.from_detector_error_model(dem, enable_correlations=True)
    with pytest.raises(ValueError):
        pymatching.Matching(dem, enable_correlations=True)
    with pytest.raises(ValueError):
        pymatching.Matching.from_detector_error_model_file(dem_file, enable_correlations=True)
