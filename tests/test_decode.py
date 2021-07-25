# Copyright 2021 Oscar Higgott

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

import pytest
from unittest.mock import patch
import numpy as np
from scipy.sparse import load_npz, csr_matrix
import pytest
import networkx as nx

from pymatching._cpp_mwpm import (BlossomFailureException, decode_match_neighbourhood,
                                  decode)
from pymatching import Matching
from pymatching.matching import _local_matching

TEST_DIR = dir_path = os.path.dirname(os.path.realpath(__file__))


def test_mwpm_decode_method():
    fn = "css_2D-toric_(4,4)_[[18,2,3]]_Hx.npz"
    H = load_npz(os.path.join(TEST_DIR, 'data', fn))
    m = Matching(H)
    n = np.zeros(H.shape[1], dtype=int)
    n[5] = 1
    n[10] = 1
    z = H.dot(n) % 2
    c = m.decode(z)
    assert(np.array_equal(c,n))


noisy_fixtures = [
    (
        np.array([
            [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        ]),
        np.array([
            [0,0,0,0,0,0,0,1,0],
            [0,0,0,0,0,0,1,0,0],
            [0,0,0,0,0,0,0,0,0]
        ]).T,
        np.array([0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0])
    ),
    (
        np.array([
            [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        ]),
        np.array([
            [0,0,0,0,0,1,0,0,0],
            [0,0,0,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0]
        ]).T,
        np.array([0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0])
    )
]


@pytest.mark.parametrize("n,z_err,c_expected", noisy_fixtures)
def test_mwpm_noisy_decode(n, z_err, c_expected):
    fn = "css_2D-toric_(4,4)_[[18,2,3]]_Hx.npz"
    H = load_npz(os.path.join(TEST_DIR, 'data', fn))
    m = Matching(H, repetitions=z_err.shape[1])
    n_all = np.cumsum(n, 0) % 2
    z_noiseless = H.dot(n_all.T) % 2
    z_noisy = (z_noiseless + z_err) % 2
    z_noisy[:, 1:] = (z_noisy[:, 1:] - z_noisy[:, :-1]) % 2
    c = m.decode(z_noisy)
    assert(np.array_equal(c, c_expected))


def test_precompute_shortest_paths():
    g = nx.Graph()
    g.add_edge(0, 1, qubit_id=0)
    g.add_edge(1, 2, qubit_id=1)
    m = Matching(g)
    assert not m.stabiliser_graph.has_computed_all_pairs_shortest_paths()
    m2 = Matching(g, precompute_shortest_paths=True)
    assert m2.stabiliser_graph.has_computed_all_pairs_shortest_paths()


def test_decode_all_neighbours():
    g = nx.Graph()
    g.add_edge(0, 1, qubit_id=0)
    g.add_edge(1, 2, qubit_id=1)
    m = Matching(g)
    noise = m.decode([1,0,1], num_neighbours=None)
    assert np.array_equal(noise, np.array([1,1]))


def test_bad_syndrome_raises_value_error():
    g = nx.Graph()
    g.add_edge(0, 1, qubit_id=0)
    g.add_edge(1, 2, qubit_id=1)
    m = Matching(g)
    with pytest.raises(TypeError):
        m.decode('test')


distance_fixtures = [
    (2,11,1),
    (3,13,2),
    (2,8,1),
    (2,98, 11)
]


@pytest.mark.parametrize("node1,node2,expected", distance_fixtures)
def test_spacetime_distance(node1, node2, expected):
    fn = "css_2D-toric_(4,4)_[[18,2,3]]_Hx.npz"
    H = load_npz(os.path.join(TEST_DIR, 'data', fn))
    m = Matching(H)
    d = m.stabiliser_graph.space_time_distance(node1, node2)
    assert(d == expected)


spacetime_path_fixtures = [
    (1,184,[1,4]),
    (2,62,[2,8])
]


@pytest.mark.parametrize("node1,node2,expected", spacetime_path_fixtures)
def test_spacetime_shortest_path(node1, node2, expected):
    fn = "css_2D-toric_(4,4)_[[18,2,3]]_Hx.npz"
    H = load_npz(os.path.join(TEST_DIR, 'data', fn))
    m = Matching(H)
    path = m.stabiliser_graph.space_time_shortest_path(node1, node2)
    assert(path == expected)


def test_3d_syndrome_raises_value_error_when_repetitions_not_set():
    H = csr_matrix(np.array([[1, 1, 0], [0, 1, 1]]))
    m = Matching(H)
    with pytest.raises(ValueError):
        m.decode(np.array([[1, 0], [0, 0]]))


def test_double_weight_matching():
    w = nx.Graph()
    w.add_edge(0, 1, qubit_id=0, weight=0.97)
    w.add_edge(2, 3, qubit_id=1, weight=1.98)
    w.add_edge(0, 2, qubit_id=2, weight=1.1)
    w.add_edge(1, 3, qubit_id=3, weight=1.2)
    m = Matching(w)
    assert(
        list(m.decode(np.array([1,1,1,1]))) == list(np.array([0,0,1,1]))
        )


def test_matching_correct():
    g = nx.Graph()
    g.add_edge(0, 1, weight=1.24, qubit_id=0)
    g.add_edge(1, 2, weight=1.31, qubit_id=1)
    g.add_edge(2, 3, weight=1.41, qubit_id=2)
    g.add_edge(0, 4, weight=1.51, qubit_id=3)
    g.add_edge(1, 5, weight=1.65, qubit_id=4)
    g.add_edge(2, 6, weight=1.15, qubit_id=5)
    g.add_edge(3, 7, weight=1.44, qubit_id=6)
    g.add_edge(4, 5, weight=1.70, qubit_id=7)
    g.add_edge(5, 6, weight=1.9, qubit_id=8)
    g.add_edge(6, 7, weight=1.12, qubit_id=9)
    g.add_edge(4, 8, weight=1.87, qubit_id=10)
    g.add_edge(5, 9, weight=1.91, qubit_id=11)
    g.add_edge(6, 10, weight=1.09, qubit_id=12)
    g.add_edge(7, 11, weight=1.21, qubit_id=13)
    g.add_edge(8, 9, weight=1.99, qubit_id=14)
    g.add_edge(9, 10, weight=1.01, qubit_id=15)
    g.add_edge(10, 11, weight=1.06, qubit_id=16)
    g.add_edge(8, 12, weight=1.16, qubit_id=17)
    g.add_edge(9, 13, weight=1.38, qubit_id=18)
    g.add_edge(10, 14, weight=1.66, qubit_id=19)
    g.add_edge(11, 15, weight=1.58, qubit_id=20)
    g.add_edge(12, 13, weight=1.12, qubit_id=21)
    g.add_edge(13, 14, weight=1.50, qubit_id=22)
    g.add_edge(14, 15, weight=1.00, qubit_id=23)

    m = Matching(g)
    assert sum(m.decode([0]*16, num_neighbours=20)) == 0
    assert sum(m.decode([0]*16, num_neighbours=None)) == 0
    z = np.zeros(16, dtype=np.uint8)
    z[0] = 1
    z[5] = 1
    z[6] = 1
    z[11] = 1
    z[14] = 1
    z[15] = 1
    assert np.array_equal(m.decode(z, num_neighbours=20).nonzero()[0], np.array([0,4,12,16,23]))
    assert np.array_equal(m.decode(z, num_neighbours=None).nonzero()[0], np.array([0,4,12,16,23]))


@pytest.mark.parametrize("cluster_size", range(3, 10, 2))
def test_local_matching_clusters(cluster_size):
    g = nx.Graph()
    qid = 0
    for i in range(cluster_size):
        g.add_edge(i, i+1, weight=1.0, qubit_id=qid)
        qid += 1
    g.add_edge(cluster_size, cluster_size+1, weight=2*cluster_size, qubit_id=qid)
    qid += 1
    for i in range(cluster_size+1, 2*cluster_size + 1):
        g.add_edge(i, i+1, weight=1.0, qubit_id=qid)
        qid += 1
    m = Matching(g)
    m.decode([1]*(cluster_size+1)*2, num_neighbours=cluster_size)


G = nx.Graph()
n = 100
for i in range(n):
    G.add_edge(i, (i+1) % n, qubit_id=i)
M = Matching(G)
defects = np.array(list(range(n)))


def test_local_matching_raises_value_error():
    with pytest.raises(ValueError):
        for x in (-10, -5, 0):
            _local_matching(M.stabiliser_graph, defects, x, False)


@pytest.mark.parametrize("num_neighbours", [2, 5, 20, 50])
def test_local_matching_raises_blossom_error(num_neighbours):
    with patch('pymatching.matching._py_decode_match_neighbourhood') as mock_decode:
        mock_decode.side_effect = BlossomFailureException
        with pytest.raises(BlossomFailureException):
            _local_matching(M.stabiliser_graph, defects, num_neighbours, False)
        assert mock_decode.call_count == 1+int(np.ceil(np.log2(n)-np.log2(num_neighbours)))


def test_local_matching_catches_blossom_errors():
    with patch('pymatching.matching._py_decode_match_neighbourhood') as mock_decode:
        mock_decode.side_effect = [BlossomFailureException]*3 + [None]
        _local_matching(M.stabiliser_graph, defects, 2, False)
        assert mock_decode.call_count == 4


def test_decoding_large_defect_id_raises_value_error():
    g = nx.Graph()
    g.add_edge(0, 1)
    g.add_edge(1, 2)
    m = Matching(g)
    with pytest.raises(ValueError):
        decode_match_neighbourhood(m.stabiliser_graph, np.array([1, 4]))
        decode(m.stabiliser_graph, np.array([1, 4]))


def test_decode_with_odd_number_of_defects():
    g = nx.Graph()
    g.add_edge(0, 1)
    g.add_edge(1, 2)
    g.add_edge(2, 0)
    m = Matching(g)
    with pytest.raises(ValueError):
        decode_match_neighbourhood(m.stabiliser_graph, np.array([1]))
    with pytest.raises(ValueError):
        decode(m.stabiliser_graph, np.array([1]))
    g.nodes[2]['is_boundary'] = True
    m2 = Matching(g)
    decode_match_neighbourhood(m2.stabiliser_graph, np.array([1]))
    decode(m2.stabiliser_graph, np.array([1]))
    g.nodes[2]['is_boundary'] = False
    g.nodes[1]['is_boundary'] = True
    m3 = Matching(g)
    decode_match_neighbourhood(m3.stabiliser_graph, np.array([1]))
    decode(m3.stabiliser_graph, np.array([1]))


def test_decode_with_multiple_components():
    g = nx.Graph()
    g.add_edge(0, 1)
    g.add_edge(1, 2)
    g.add_edge(2, 0)
    g.add_edge(3, 4)
    g.add_edge(4, 5)
    g.add_edge(3, 5)

    m = Matching(g)
    for z in (np.array([0]), np.arange(6)):
        with pytest.raises(ValueError):
            decode_match_neighbourhood(m.stabiliser_graph, z)
        with pytest.raises(ValueError):
            decode(m.stabiliser_graph, z)

    g.nodes[0]['is_boundary'] = True
    m2 = Matching(g)
    decode_match_neighbourhood(m2.stabiliser_graph, np.array([1]))
    decode(m2.stabiliser_graph, np.array([1]))
    for z in (np.arange(6), np.array([3])):
        with pytest.raises(ValueError):
            decode_match_neighbourhood(m2.stabiliser_graph, z)
        with pytest.raises(ValueError):
            decode(m2.stabiliser_graph, z)

    g.nodes[4]['is_boundary'] = True
    m3 = Matching(g)
    for z in (np.array([0]), np.arange(6), np.array([3]), np.array([1,3])):
        decode_match_neighbourhood(m3.stabiliser_graph, z)
        decode(m3.stabiliser_graph, z)
