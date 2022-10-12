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

import numpy as np
import pymatching
from scipy.sparse import csc_matrix, csr_matrix
import pytest
import networkx as nx
import retworkx as rx
import matplotlib.pyplot as plt
import stim
import doctest

from pymatching._cpp_pymatching import MatchingGraph
from pymatching import Matching


def repetition_code(n):
    row_ind, col_ind = zip(*((i, j) for i in range(n) for j in (i, (i + 1) % n)))
    data = np.ones(2 * n, dtype=np.uint8)
    return csr_matrix((data, (row_ind, col_ind)))


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

    arr = m.decode_to_matched_detection_events_array(syndrome)
    assert np.array_equal(arr, np.array([[2, -1], [10, 12], [18, -1]]))

    d = m.decode_to_matched_detection_events_dict(syndrome)
    assert d == {
        2: None,
        10: 12,
        12: 10,
        18: None
    }
