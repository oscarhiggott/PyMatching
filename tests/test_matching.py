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
import pytest
from scipy.sparse import csc_matrix, csr_matrix
import pytest
import networkx as nx
import matplotlib.pyplot as plt

from pymatching._cpp_mwpm import WeightedStabiliserGraph
from pymatching import Matching


def test_bad_qubit_id_raises_value_error():
    g = nx.Graph()
    g.add_edge(0,1, qubit_id='test')
    with pytest.raises(ValueError):
        m = Matching(g)
    g = nx.Graph()
    g.add_edge(0,1, qubit_id=[[1],[2]])
    with pytest.raises(ValueError):
        m = Matching(g)


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
    assert(m.num_detectors == 6)
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


def test_mwpm_from_networkx():
    g = nx.Graph()
    g.add_edge(0, 1, qubit_id=0)
    g.add_edge(0, 2, qubit_id=1)
    g.add_edge(1, 2, qubit_id=2)
    m = Matching(g)
    assert(isinstance(m.stabiliser_graph, WeightedStabiliserGraph))
    assert(m.num_detectors == 3)
    assert(m.num_qubits == 3)
    assert(m.stabiliser_graph.distance(0,2) == 1)
    assert(m.stabiliser_graph.shortest_path(0,2) == [0,2])

    g = nx.Graph()
    g.add_edge(0, 1)
    g.add_edge(0, 2)
    g.add_edge(1, 2)
    m = Matching(g)
    assert(isinstance(m.stabiliser_graph, WeightedStabiliserGraph))
    assert(m.num_detectors == 3)
    assert(m.num_qubits == 0)
    assert(m.stabiliser_graph.distance(0,2) == 1)
    assert(m.stabiliser_graph.shortest_path(0,2) == [0,2])

    g = nx.Graph()
    g.add_edge(0, 1, weight=1.5)
    g.add_edge(0, 2, weight=1.7)
    g.add_edge(1, 2, weight=1.2)
    m = Matching(g)
    assert(isinstance(m.stabiliser_graph, WeightedStabiliserGraph))
    assert(m.num_detectors == 3)
    assert(m.num_qubits == 0)
    assert(m.stabiliser_graph.distance(0,2) == pytest.approx(1.7))
    assert(m.stabiliser_graph.shortest_path(0,2) == [0,2])


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


def test_load_matching_from_dense_array():
    H = np.array([[1, 1, 0], [0, 1, 1]])
    m = Matching()
    m.load_from_check_matrix(H)
