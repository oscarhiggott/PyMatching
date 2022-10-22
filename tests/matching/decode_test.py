# Copyright 2020 Oscar Higgott

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#      http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import pathlib
import numpy as np
import stim
from scipy.sparse import csc_matrix
import pytest
import networkx as nx
import os

from pymatching import Matching

THIS_DIR = pathlib.Path(__file__).parent.resolve()
DATA_DIR = os.path.join(pathlib.Path(THIS_DIR).parent.parent.absolute(), "data")


def repetition_code(n):
    row_ind, col_ind = zip(*((i, j) for i in range(n) for j in (i, (i + 1) % n)))
    data = np.ones(2 * n, dtype=np.uint8)
    return csc_matrix((data, (row_ind, col_ind)))


weight_fixtures = [
    10, 15, 20, 100
]


@pytest.mark.parametrize("n", weight_fixtures)
def test_matching_weight(n):
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
    assert d == {
        2: None,
        10: 12,
        12: 10,
        18: None
    }


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


def test_surface_code_solution_weights():
    dem = stim.DetectorErrorModel.from_file(os.path.join(DATA_DIR, "surface_code_rotated_memory_x_13_0.01.dem"))
    m = Matching.from_detector_error_model(dem)
    shots = stim.read_shot_data_file(path=os.path.join(DATA_DIR, "surface_code_rotated_memory_x_13_0.01_1000_shots.dets"),
                                     format="dets", num_detectors=m.num_detectors,
                                     num_observables=m.num_fault_ids)
    with open(os.path.join(
            DATA_DIR,
            "surface_code_rotated_memory_x_13_0.01_1000_shots_no_buckets_weights_pymatchingv0.7_exact.txt"),
            "r") as f:
        expected_weights = [float(w) for w in f.readlines()]
    with open(os.path.join(
            DATA_DIR,
            "surface_code_rotated_memory_x_13_0.01_1000_shots_no_buckets_predictions_pymatchingv0.7_exact.txt"),
            "r") as f:
        expected_observables = [float(w) for w in f.readlines()]
    assert shots.shape == (1000, m.num_detectors + m.num_fault_ids)
    weights = []
    predicted_observables = []
    for i in range(min(shots.shape[0], 1000)):
        prediction, weight = m.decode(shots[i, 0:-m.num_fault_ids], return_weight=True)
        weights.append(weight)
        predicted_observables.append(prediction)
    for i in range(len(weights)):
        assert weights[i] == pytest.approx(expected_weights[i], rel=1e-8)
    assert predicted_observables == expected_observables[0:len(predicted_observables)]


def test_deprecated_position_arguments_raise_deprecation_warning():
    m = Matching()
    m.add_edge(0, 1, fault_ids={0}, weight=4)
    m.add_edge(1, 2, fault_ids={1}, weight=9)
    with pytest.warns(DeprecationWarning):
        correction = m.decode([0, 1, 1], 3)
        assert np.array_equal(correction, np.array([0, 1]))
    with pytest.warns(DeprecationWarning):
        correction, weight = m.decode([0, 1, 1], 3, True)
        assert np.array_equal(correction, np.array([0, 1]))
        assert weight == 9.0
