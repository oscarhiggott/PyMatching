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

import os

import pytest
from unittest.mock import patch
import numpy as np
from scipy.sparse import csc_matrix, load_npz, csr_matrix
import pytest
import networkx as nx
import matplotlib.pyplot as plt

from pymatching._cpp_mwpm import WeightedStabiliserGraph, BlossomFailureException
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
    z_noisy[:,1:] = (z_noisy[:,1:] - z_noisy[:,:-1]) % 2
    c = m.decode(z_noisy)
    assert(np.array_equal(c, c_expected))


def test_bad_qubit_id_raises_value_error():
    g = nx.Graph()
    g.add_edge(0,1, qubit_id='test')
    with pytest.raises(ValueError):
        m = Matching(g)
    g = nx.Graph()
    g.add_edge(0,1, qubit_id=[[1],[2]])
    with pytest.raises(ValueError):
        m = Matching(g)


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
    with pytest.raises(ValueError):
        noise = m.decode('test')


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


def test_boundary_from_check_matrix():
    H = csr_matrix(np.array([[1,1,0,0,0],[0,1,1,0,0],
                             [0,0,1,1,0],[0,0,0,1,1]]))
    m = Matching(H)
    assert m.boundary == [4]
    assert np.array_equal(m.decode(np.array([1,0,0,0])), np.array([1,0,0,0,0]))
    assert np.array_equal(m.decode(np.array([0,1,0,0])), np.array([1,1,0,0,0]))
    assert np.array_equal(m.decode(np.array([0,1,1,0])), np.array([0,0,1,0,0]))
    assert np.array_equal(m.decode(np.array([0,0,1,0])), np.array([0,0,0,1,1]))


def test_boundary_from_networkx():
    g = nx.Graph()
    g.add_edge(4,0, qubit_id=0)
    g.add_edge(0,1, qubit_id=1)
    g.add_edge(1,2, qubit_id=2)
    g.add_edge(2,3, qubit_id=3)
    g.add_edge(3,4, qubit_id=4)
    g.nodes()[4]['is_boundary'] = True
    m = Matching(g)
    assert m.boundary == [4]
    assert np.array_equal(m.decode(np.array([1,0,0,0])), np.array([1,0,0,0,0]))
    assert np.array_equal(m.decode(np.array([0,1,0,0])), np.array([1,1,0,0,0]))
    assert np.array_equal(m.decode(np.array([0,1,1,0])), np.array([0,0,1,0,0]))
    assert np.array_equal(m.decode(np.array([0,0,1,0])), np.array([0,0,0,1,1]))


def test_boundaries_from_networkx():
    g = nx.Graph()
    g.add_edge(0,1, qubit_id=0)
    g.add_edge(1,2, qubit_id=1)
    g.add_edge(2,3, qubit_id=2)
    g.add_edge(3,4, qubit_id=3)
    g.add_edge(4,5, qubit_id=4)
    g.add_edge(0,5, qubit_id=-1, weight=0.0)
    g.nodes()[0]['is_boundary'] = True
    g.nodes()[5]['is_boundary'] = True
    m = Matching(g)
    assert m.boundary == [0,5]
    assert np.array_equal(m.decode(np.array([0,1,0,0,0,0])), np.array([1,0,0,0,0]))
    assert np.array_equal(m.decode(np.array([0,0,1,0,0])), np.array([1,1,0,0,0]))
    assert np.array_equal(m.decode(np.array([0,0,1,1,0])), np.array([0,0,1,0,0]))
    assert np.array_equal(m.decode(np.array([0,0,0,1,0])), np.array([0,0,0,1,1]))


def test_nonzero_matrix_elements_not_one_raises_value_error():
    H = csr_matrix(np.array([[0,1.01,1.01],[1.01,1.01,0]]))
    with pytest.raises(ValueError):
        Matching(H)


def test_too_many_checks_per_qubit_raises_value_error():
    H = csr_matrix(np.array([[1,1,0,0],[1,0,1,0],[1,0,0,1]]))
    with pytest.raises(ValueError):
        Matching(H)


def test_negative_weight_raises_value_error():
    g = nx.Graph()
    g.add_edge(0,1,weight=-1)
    with pytest.raises(ValueError):
        Matching(g)
    with pytest.raises(ValueError):
        Matching(csr_matrix([[1,1,0],[0,1,1]]), spacelike_weights=np.array([1,1,-1]))


def test_odd_3d_syndrome_raises_value_error():
    H = csr_matrix(np.array([[1,1,0],[0,1,1]]))
    m = Matching(H)
    with pytest.raises(ValueError):
        m.decode(np.array([[1,0],[0,0]]))


def test_add_noise_to_unweighted_returns_none():
    m = Matching(csr_matrix(np.array([[1,1,0],[0,1,1]])))
    assert m.add_noise() == None
    m = Matching(csr_matrix(np.array([[1,1,0],[0,1,1]])), 
             error_probabilities=np.array([0.5,0.7,-0.1]))
    assert m.add_noise() == None


def test_error_probability_from_array():
    H = csr_matrix(np.array([[1,1,0,0,0],[0,1,1,0,0],
                             [0,0,1,1,0],[0,0,0,1,1]]))
    m = Matching(H, error_probabilities=np.array([0.,0.,0.,0.,1.]))
    assert np.array_equal(m.add_noise()[0], np.array([0,0,0,0,1]))
    assert np.array_equal(m.add_noise()[1], np.array([0,0,0,1,1]))
    m = Matching(H, error_probabilities=np.array([0.,0.,0.,0.,0.]))
    assert np.array_equal(m.add_noise()[0], np.array([0,0,0,0,0]))
    assert np.array_equal(m.add_noise()[1], np.array([0,0,0,0,0]))
    m = Matching(H, error_probabilities=0.0)
    assert np.array_equal(m.add_noise()[0], np.array([0,0,0,0,0]))
    assert np.array_equal(m.add_noise()[1], np.array([0,0,0,0,0]))
    m = Matching(H, error_probabilities=1.0)
    assert np.array_equal(m.add_noise()[0], np.array([1,1,1,1,1]))
    assert np.array_equal(m.add_noise()[1], np.array([0,0,0,0,0]))


def test_weighted_mwpm_from_array():
    H = csc_matrix([[1,0],[1,1],[0,1]])
    m = Matching(H, spacelike_weights=np.array([1., 2.]))
    assert m.stabiliser_graph.distance(0, 1) == 1.
    assert m.stabiliser_graph.distance(1, 2) == 2.
    with pytest.raises(ValueError):
        m = Matching(H, spacelike_weights=np.array([1.]))
    with pytest.raises(ValueError):
        m = Matching(H, spacelike_weights=np.array([1., -2.]))


def test_unweighted_stabiliser_graph_from_networkx():
    w = nx.Graph()
    w.add_edge(0, 1, qubit_id=0, weight=7.0)
    w.add_edge(0, 5, qubit_id=1, weight=14.0)
    w.add_edge(0, 2, qubit_id=2, weight=9.0)
    w.add_edge(1, 2, qubit_id=-1, weight=10.0)
    w.add_edge(1, 3, qubit_id=3, weight=15.0)
    w.add_edge(2, 5, qubit_id=4, weight=2.0)
    w.add_edge(2, 3, qubit_id=-1, weight=11.0)
    w.add_edge(3, 4, qubit_id=5, weight=6.0)
    w.add_edge(4, 5, qubit_id=6, weight=9.0)
    m = Matching(w)
    assert(m.num_qubits == 7)
    assert(m.num_stabilisers == 6)
    assert(m.stabiliser_graph.shortest_path(3, 5) == [3, 2, 5])
    assert(m.stabiliser_graph.distance(5, 0) == pytest.approx(11.0))
    assert(np.array_equal(
        m.decode(np.array([1,0,1,0,0,0])),
        np.array([0,0,1,0,0,0,0]))
    )
    with pytest.raises(ValueError):
        m.decode(np.array([1,1,0]))
    with pytest.raises(ValueError):
        m.decode(np.array([1,1,1,0,0,0]))
    assert(np.array_equal(
        m.decode(np.array([1,0,0,0,0,1])),
        np.array([0,0,1,0,1,0,0]))
    )
    assert(np.array_equal(
        m.decode(np.array([0,1,0,0,0,1])),
        np.array([0,0,0,0,1,0,0]))
    )


def test_mwmpm_from_networkx():
    g = nx.Graph()
    g.add_edge(0, 1, qubit_id=0)
    g.add_edge(0, 2, qubit_id=1)
    g.add_edge(1, 2, qubit_id=2)
    m = Matching(g)
    assert(isinstance(m.stabiliser_graph, WeightedStabiliserGraph))
    assert(m.num_stabilisers == 3)
    assert(m.num_qubits == 3)
    assert(m.stabiliser_graph.distance(0,2) == 1)
    assert(m.stabiliser_graph.shortest_path(0,2) == [0,2])

    g = nx.Graph()
    g.add_edge(0, 1)
    g.add_edge(0, 2)
    g.add_edge(1, 2)
    m = Matching(g)
    assert(isinstance(m.stabiliser_graph, WeightedStabiliserGraph))
    assert(m.num_stabilisers == 3)
    assert(m.num_qubits == 0)
    assert(m.stabiliser_graph.distance(0,2) == 1)
    assert(m.stabiliser_graph.shortest_path(0,2) == [0,2])

    g = nx.Graph()
    g.add_edge(0, 1, weight=1.5)
    g.add_edge(0, 2, weight=1.7)
    g.add_edge(1, 2, weight=1.2)
    m = Matching(g)
    assert(isinstance(m.stabiliser_graph, WeightedStabiliserGraph))
    assert(m.num_stabilisers == 3)
    assert(m.num_qubits == 0)
    assert(m.stabiliser_graph.distance(0,2) == pytest.approx(1.7))
    assert(m.stabiliser_graph.shortest_path(0,2) == [0,2])


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


def test_repr():
    g = nx.Graph()
    g.add_edge(0, 1, qubit_id=0)
    g.add_edge(1, 2, qubit_id=1)
    g.add_edge(2, 3, qubit_id=2)
    g.nodes[0]['is_boundary'] = True
    g.nodes[3]['is_boundary'] = True
    g.add_edge(0, 3, weight=0.0)
    m = Matching(g)
    assert m.__repr__() == ("<pymatching.Matching object with 3 qubits, "
                            "2 stabilisers, 2 boundary nodes, and 4 edges>")


def test_wrong_connected_components_raises_value_error():
    g = nx.Graph()
    g.add_edge(0, 1, qubit_id=0)
    g.add_edge(1, 2, qubit_id=1)
    g.add_edge(2, 0, qubit_id=2)
    g.add_edge(3, 4, qubit_id=3)
    g.add_edge(4, 5, qubit_id=4)
    g.add_edge(5, 3, qubit_id=5)
    with pytest.raises(ValueError):
        Matching(g)
    g = nx.Graph()
    g.add_edge(0, 1, qubit_id=0)
    g.add_edge(1, 2, qubit_id=1)
    g.add_edge(2, 0, qubit_id=2)
    m = Matching(g)
    assert m.stabiliser_graph.get_num_connected_components() == 1


def test_high_qubit_id_raises_value_error():
    g = nx.Graph()
    g.add_edge(0,1,qubit_id=1)
    with pytest.raises(ValueError):
        Matching(g)


def test_high_node_id_raises_value_error():
    g = nx.Graph()
    g.add_edge(1, 2)
    with pytest.raises(ValueError):
        Matching(g)


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


@pytest.mark.parametrize("cluster_size", range(3,10,2))
def test_local_matching_connected(cluster_size):
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


def test_matching_edges():
    g = nx.Graph()
    g.add_edge(0, 1, qubit_id=0, weight=1.1, error_probability=0.1)
    g.add_edge(1, 2, qubit_id=1, weight=2.1, error_probability=0.2)
    g.add_edge(2, 3, qubit_id={2,3}, weight=0.9, error_probability=0.3)
    g.nodes[0]['is_boundary'] = True
    g.nodes[3]['is_boundary'] = True
    g.add_edge(0, 3, weight=0.0)
    m = Matching(g)
    es = list(m.edges())
    expected_edges = [
        (0,1,{'qubit_id': {0}, 'weight': 1.1, 'error_probability': 0.1}),
        (0,3,{'qubit_id': set(), 'weight': 0.0, 'error_probability': -1.0}),
        (1,2,{'qubit_id': {1}, 'weight': 2.1, 'error_probability': 0.2}),
        (2,3,{'qubit_id': {2,3}, 'weight': 0.9, 'error_probability': 0.3})
        
    ]
    assert es == expected_edges


def test_matching_to_networkx():
    g = nx.Graph()
    g.add_edge(0, 1, qubit_id={0}, weight=1.1, error_probability=0.1)
    g.add_edge(1, 2, qubit_id={1}, weight=2.1, error_probability=0.2)
    g.add_edge(2, 3, qubit_id={2,3}, weight=0.9, error_probability=0.3)
    g.nodes[0]['is_boundary'] = True
    g.nodes[3]['is_boundary'] = True
    g.add_edge(0, 3, weight=0.0)
    m = Matching(g)

    g.edges[(0,3)]['qubit_id'] = set()
    g.edges[(0,3)]['error_probability'] = -1.0
    g.nodes[1]['is_boundary'] = False
    g.nodes[2]['is_boundary'] = False
    
    g2 = m.to_networkx()

    assert g.nodes(data=True) == g2.nodes(data=True)
    gedges = [({s,t},d) for (s, t, d) in g.edges(data=True)]
    g2edges = [({s,t},d) for (s, t, d) in g2.edges(data=True)]
    assert sorted(gedges) == sorted(g2edges)


def test_draw_matching():
    g = nx.Graph()
    g.add_edge(0, 1, qubit_id={0}, weight=1.1, error_probability=0.1)
    g.add_edge(1, 2, qubit_id={1}, weight=2.1, error_probability=0.2)
    g.add_edge(2, 3, qubit_id={2,3}, weight=0.9, error_probability=0.3)
    g.nodes[0]['is_boundary'] = True
    g.nodes[3]['is_boundary'] = True
    g.add_edge(0, 3, weight=0.0)
    m = Matching(g)
    plt.figure()
    m.draw()


G = nx.Graph()
n = 100
for i in range(n):
    G.add_edge(i, (i+1) % n, qubit_id=i)
M = Matching(G)
defects = np.array(list(range(n)))


def test_local_matching_raises_value_error():
    with pytest.raises(ValueError):
        for x in (-10, -5, 0):
            _local_matching(M.stabiliser_graph, defects, 20, False, x)


@pytest.mark.parametrize("num_attempts", [1,3,5])
def test_local_matching_raises_blossom_error(num_attempts):
    with patch('pymatching.matching._py_decode_match_neighbourhood') as mock_decode:
        mock_decode.side_effect = BlossomFailureException
        with pytest.raises(BlossomFailureException):
            _local_matching(M.stabiliser_graph, defects, 20, False, num_attempts)
        assert mock_decode.call_count == num_attempts


def test_local_matching_catches_blossom_errors():
    with patch('pymatching.matching._py_decode_match_neighbourhood') as mock_decode:
        mock_decode.side_effect = [BlossomFailureException]*3 + [0]
        _local_matching(M.stabiliser_graph, defects, 20, False, 5)
        assert mock_decode.call_count == 4
