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

from pymatching._cpp_mwpm import MatchingGraph
from pymatching import Matching


def test_bad_fault_ids_raises_value_error():
    g = nx.Graph()
    g.add_edge(0,1, fault_ids='test')
    with pytest.raises(ValueError):
        m = Matching(g)
    g = nx.Graph()
    g.add_edge(0,1, fault_ids=[[1],[2]])
    with pytest.raises(ValueError):
        m = Matching(g)


def test_boundary_from_check_matrix():
    H = csr_matrix(np.array([[1,1,0,0,0],[0,1,1,0,0],
                             [0,0,1,1,0],[0,0,0,1,1]]))
    m = Matching(H)
    assert m.boundary == {4}
    assert np.array_equal(m.decode(np.array([1,0,0,0])), np.array([1,0,0,0,0]))
    assert np.array_equal(m.decode(np.array([0,1,0,0])), np.array([1,1,0,0,0]))
    assert np.array_equal(m.decode(np.array([0,1,1,0])), np.array([0,0,1,0,0]))
    assert np.array_equal(m.decode(np.array([0,0,1,0])), np.array([0,0,0,1,1]))


def test_boundary_from_networkx():
    g = nx.Graph()
    g.add_edge(4,0, fault_ids=0)
    g.add_edge(0,1, fault_ids=1)
    g.add_edge(1,2, fault_ids=2)
    g.add_edge(2,3, fault_ids=3)
    g.add_edge(3,4, fault_ids=4)
    g.nodes()[4]['is_boundary'] = True
    m = Matching(g)
    assert m.boundary == {4}
    assert np.array_equal(m.decode(np.array([1,0,0,0])), np.array([1,0,0,0,0]))
    assert np.array_equal(m.decode(np.array([0,1,0,0])), np.array([1,1,0,0,0]))
    assert np.array_equal(m.decode(np.array([0,1,1,0])), np.array([0,0,1,0,0]))
    assert np.array_equal(m.decode(np.array([0,0,1,0])), np.array([0,0,0,1,1]))


def test_boundaries_from_networkx():
    g = nx.Graph()
    g.add_edge(0,1, fault_ids=0)
    g.add_edge(1,2, fault_ids=1)
    g.add_edge(2,3, fault_ids=2)
    g.add_edge(3,4, fault_ids=3)
    g.add_edge(4,5, fault_ids=4)
    g.add_edge(0,5, fault_ids=-1, weight=0.0)
    g.nodes()[0]['is_boundary'] = True
    g.nodes()[5]['is_boundary'] = True
    m = Matching(g)
    assert m.boundary == {0, 5}
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


def test_wrong_check_matrix_type_raises_type_error():
    with pytest.raises(TypeError):
        Matching("test")
    m = Matching()
    with pytest.raises(TypeError):
        m.load_from_check_matrix("test")


def test_wrong_networkx_graph_type_raises_type_error():
    m = Matching()
    with pytest.raises(TypeError):
        m.load_from_networkx("test")


def test_error_probability_from_array():
    H = csr_matrix(np.array([[1, 1, 0, 0, 0],[0, 1, 1, 0, 0],
                             [0, 0, 1, 1, 0],[0, 0, 0, 1, 1]]))
    m = Matching(H, error_probabilities=np.array([0., 0., 0., 0., 1.]))
    assert np.array_equal(m.add_noise()[0], np.array([0, 0, 0, 0, 1]))
    assert np.array_equal(m.add_noise()[1], np.array([0, 0, 0, 1, 0]))
    m = Matching(H, error_probabilities=np.array([0., 0., 0., 0., 0.]))
    assert np.array_equal(m.add_noise()[0], np.array([0, 0, 0, 0, 0]))
    assert np.array_equal(m.add_noise()[1], np.array([0, 0, 0, 0, 0]))
    m = Matching(H, error_probabilities=0.0)
    assert np.array_equal(m.add_noise()[0], np.array([0, 0, 0, 0, 0]))
    assert np.array_equal(m.add_noise()[1], np.array([0, 0, 0, 0, 0]))
    m = Matching(H, error_probabilities=1.0)
    assert np.array_equal(m.add_noise()[0], np.array([1, 1, 1, 1, 1]))
    assert np.array_equal(m.add_noise()[1], np.array([0, 0, 0, 0, 0]))


def test_weighted_mwpm_from_array():
    H = csc_matrix([[1,0],[1,1],[0,1]])
    m = Matching(H, spacelike_weights=np.array([1., 2.]))
    assert m.matching_graph.distance(0, 1) == 1.
    assert m.matching_graph.distance(1, 2) == 2.
    with pytest.raises(ValueError):
        m = Matching(H, spacelike_weights=np.array([1.]))
    with pytest.raises(ValueError):
        m = Matching(H, spacelike_weights=np.array([1., -2.]))


def test_unweighted_stabiliser_graph_from_networkx():
    w = nx.Graph()
    w.add_edge(0, 1, fault_ids=0, weight=7.0)
    w.add_edge(0, 5, fault_ids=1, weight=14.0)
    w.add_edge(0, 2, fault_ids=2, weight=9.0)
    w.add_edge(1, 2, fault_ids=-1, weight=10.0)
    w.add_edge(1, 3, fault_ids=3, weight=15.0)
    w.add_edge(2, 5, fault_ids=4, weight=2.0)
    w.add_edge(2, 3, fault_ids=-1, weight=11.0)
    w.add_edge(3, 4, fault_ids=5, weight=6.0)
    w.add_edge(4, 5, fault_ids=6, weight=9.0)
    m = Matching(w)
    assert(m.num_fault_ids == 7)
    assert(m.num_detectors == 6)
    assert(m.matching_graph.shortest_path(3, 5) == [3, 2, 5])
    assert(m.matching_graph.distance(5, 0) == pytest.approx(11.0))
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
    g.add_edge(0, 1, fault_ids=0)
    g.add_edge(0, 2, fault_ids=1)
    g.add_edge(1, 2, fault_ids=2)
    m = Matching(g)
    assert(isinstance(m.matching_graph, MatchingGraph))
    assert(m.num_detectors == 3)
    assert(m.num_fault_ids == 3)
    assert(m.matching_graph.distance(0,2) == 1)
    assert(m.matching_graph.shortest_path(0,2) == [0,2])

    g = nx.Graph()
    g.add_edge(0, 1)
    g.add_edge(0, 2)
    g.add_edge(1, 2)
    m = Matching(g)
    assert(isinstance(m.matching_graph, MatchingGraph))
    assert(m.num_detectors == 3)
    assert(m.num_fault_ids == 0)
    assert(m.matching_graph.distance(0,2) == 1)
    assert(m.matching_graph.shortest_path(0,2) == [0,2])

    g = nx.Graph()
    g.add_edge(0, 1, weight=1.5)
    g.add_edge(0, 2, weight=1.7)
    g.add_edge(1, 2, weight=1.2)
    m = Matching(g)
    assert(isinstance(m.matching_graph, MatchingGraph))
    assert(m.num_detectors == 3)
    assert(m.num_fault_ids == 0)
    assert(m.matching_graph.distance(0,2) == pytest.approx(1.7))
    assert(m.matching_graph.shortest_path(0,2) == [0,2])


def test_repr():
    g = nx.Graph()
    g.add_edge(0, 1, fault_ids=0)
    g.add_edge(1, 2, fault_ids=1)
    g.add_edge(2, 3, fault_ids=2)
    g.nodes[0]['is_boundary'] = True
    g.nodes[3]['is_boundary'] = True
    g.add_edge(0, 3, weight=0.0)
    m = Matching(g)
    assert m.__repr__() == ("<pymatching.Matching object with "
                            "2 detectors, 2 boundary nodes, and 4 edges>")


def test_matching_edges_from_networkx():
    g = nx.Graph()
    g.add_edge(0, 1, fault_ids=0, weight=1.1, error_probability=0.1)
    g.add_edge(1, 2, fault_ids=1, weight=2.1, error_probability=0.2)
    g.add_edge(2, 3, fault_ids={2,3}, weight=0.9, error_probability=0.3)
    g.nodes[0]['is_boundary'] = True
    g.nodes[3]['is_boundary'] = True
    g.add_edge(0, 3, weight=0.0)
    m = Matching(g)
    es = list(m.edges())
    expected_edges = [
        (0,1,{'fault_ids': {0}, 'weight': 1.1, 'error_probability': 0.1}),
        (0,3,{'fault_ids': set(), 'weight': 0.0, 'error_probability': -1.0}),
        (1,2,{'fault_ids': {1}, 'weight': 2.1, 'error_probability': 0.2}),
        (2,3,{'fault_ids': {2,3}, 'weight': 0.9, 'error_probability': 0.3})
        
    ]
    assert es == expected_edges


def test_qubit_id_accepted_via_networkx():
    g = nx.Graph()
    g.add_edge(0, 1, qubit_id=0, weight=1.1, error_probability=0.1)
    g.add_edge(1, 2, qubit_id=1, weight=2.1, error_probability=0.2)
    g.add_edge(2, 3, qubit_id={2, 3}, weight=0.9, error_probability=0.3)
    g.nodes[0]['is_boundary'] = True
    g.nodes[3]['is_boundary'] = True
    g.add_edge(0, 3, weight=0.0)
    m = Matching(g)
    es = list(m.edges())
    expected_edges = [
        (0, 1, {'fault_ids': {0}, 'weight': 1.1, 'error_probability': 0.1}),
        (0, 3, {'fault_ids': set(), 'weight': 0.0, 'error_probability': -1.0}),
        (1, 2, {'fault_ids': {1}, 'weight': 2.1, 'error_probability': 0.2}),
        (2, 3, {'fault_ids': {2, 3}, 'weight': 0.9, 'error_probability': 0.3})

    ]
    assert es == expected_edges


def test_qubit_id_accepted_using_add_edge():
    m = Matching()
    m.add_edge(0, 1, qubit_id=0)
    m.add_edge(1, 2, qubit_id={1, 2})
    es = list(m.edges())
    expected_edges = [
        (0, 1, {'fault_ids': {0}, 'weight': 1.0, 'error_probability': -1.0}),
        (1, 2, {'fault_ids': {1, 2}, 'weight': 1.0, 'error_probability': -1.0})
    ]
    assert es == expected_edges


def test_add_edge_raises_value_error_if_qubit_id_and_fault_ids_both_supplied():
    with pytest.raises(ValueError):
        m = Matching()
        m.add_edge(0, 1, qubit_id=0, fault_ids=0)
        m.add_edge(1, 2, qubit_id=1, fault_ids=1)


def test_load_from_networkx_raises_value_error_if_qubit_id_and_fault_ids_both_supplied():
    with pytest.raises(ValueError):
        g = nx.Graph()
        g.add_edge(0, 1, qubit_id=0, fault_ids=0)
        g.add_edge(1, 2, qubit_id=1, fault_ids=1)
        m = Matching()
        m.load_from_networkx(g)


def test_matching_to_networkx():
    g = nx.Graph()
    g.add_edge(0, 1, fault_ids={0}, weight=1.1, error_probability=0.1)
    g.add_edge(1, 2, fault_ids={1}, weight=2.1, error_probability=0.2)
    g.add_edge(2, 3, fault_ids={2,3}, weight=0.9, error_probability=0.3)
    g.nodes[0]['is_boundary'] = True
    g.nodes[3]['is_boundary'] = True
    g.add_edge(0, 3, weight=0.0)
    m = Matching(g)

    g.edges[(0,3)]['fault_ids'] = set()
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
    g.add_edge(0, 1, fault_ids={0}, weight=1.1, error_probability=0.1)
    g.add_edge(1, 2, fault_ids={1}, weight=2.1, error_probability=0.2)
    g.add_edge(2, 3, fault_ids={2,3}, weight=0.9, error_probability=0.3)
    g.nodes[0]['is_boundary'] = True
    g.nodes[3]['is_boundary'] = True
    g.add_edge(0, 3, weight=0.0)
    m = Matching(g)
    plt.figure()
    m.draw()


def test_add_edge():
    m = Matching()
    m.add_edge(0, 1)
    m.add_edge(1, 2)
    assert m.num_nodes == 3
    assert m.num_edges == 2

    m = Matching()
    m.add_edge(0, 1, weight=0.123, error_probability=0.6)
    m.add_edge(1, 2, weight=0.6, error_probability=0.3, fault_ids=0)
    m.add_edge(2, 3, weight=0.01, error_probability=0.5, fault_ids={1, 2})
    expected = [(0, 1, {'fault_ids': set(), 'weight': 0.123, 'error_probability': 0.6}),
                (1, 2, {'fault_ids': {0}, 'weight': 0.6, 'error_probability': 0.3}),
                (2, 3, {'fault_ids': {1, 2}, 'weight': 0.01, 'error_probability': 0.5})]
    assert m.edges() == expected


def test_load_matching_from_dense_array():
    H = np.array([[1, 1, 0], [0, 1, 1]])
    m = Matching()
    m.load_from_check_matrix(H)


@pytest.mark.parametrize("t_weights,expected_edges",
                         [
                             (
                                     [0.5, 1.5],
                                     {((0, 1), 0.7), ((2, 3), 0.7), ((4, 5), 0.7),
                                      ((0, 6), 0.3), ((2, 6), 0.3), ((4, 6), 0.3),
                                      ((1, 6), 0.9), ((3, 6), 0.9), ((5, 6), 0.9),
                                      ((0, 2), 0.5), ((2, 4), 0.5), ((1, 3), 1.5), ((3, 5), 1.5)}
                             ),
                             (
                                     np.array([0.5, 1.5]),
                                     {((0, 1), 0.7), ((2, 3), 0.7), ((4, 5), 0.7),
                                      ((0, 6), 0.3), ((2, 6), 0.3), ((4, 6), 0.3),
                                      ((1, 6), 0.9), ((3, 6), 0.9), ((5, 6), 0.9),
                                      ((0, 2), 0.5), ((2, 4), 0.5), ((1, 3), 1.5), ((3, 5), 1.5)}
                             ),
                             (
                                 1.2,
                                 {((0, 1), 0.7), ((2, 3), 0.7), ((4, 5), 0.7),
                                  ((0, 6), 0.3), ((2, 6), 0.3), ((4, 6), 0.3),
                                  ((1, 6), 0.9), ((3, 6), 0.9), ((5, 6), 0.9),
                                  ((0, 2), 1.2), ((2, 4), 1.2), ((1, 3), 1.2), ((3, 5), 1.2)}
                             )
                         ]
                         )
def test_timelike_weights(t_weights, expected_edges):
    H = np.array([[1, 1, 0], [0, 1, 1]])
    m = Matching()
    m.load_from_check_matrix(H, spacelike_weights=np.array([0.3, 0.7, 0.9]),
                             timelike_weights=t_weights, repetitions=3)
    es = set((tuple(sorted([u, v])), d["weight"]) for u, v, d in m.edges())
    assert es == expected_edges


@pytest.mark.parametrize("t_weights", [[0.1, 0.01, 3], "A"])
def test_wrong_timelike_weights_raises_valueerror(t_weights):
    H = np.array([[1, 1, 0], [0, 1, 1]])
    with pytest.raises(ValueError):
        m = Matching()
        m.load_from_check_matrix(H, spacelike_weights=np.array([0.3, 0.7, 0.9]),
                                 timelike_weights=t_weights, repetitions=3)


@pytest.mark.parametrize("p_meas,expected_edges,repetitions",
                         [
                             (
                                [0.15, 0.25],
                                {((0, 1), 0.2), ((2, 3), 0.2), ((4, 5), 0.2),
                                 ((0, 6), 0.1), ((2, 6), 0.1), ((4, 6), 0.1),
                                 ((1, 6), 0.3), ((3, 6), 0.3), ((5, 6), 0.3),
                                 ((0, 2), 0.15), ((2, 4), 0.15), ((1, 3), 0.25), ((3, 5), 0.25)},
                                 3
                              ),
                             (
                                 np.array([0.15, 0.25]),
                                 {((0, 1), 0.2), ((2, 3), 0.2),
                                  ((0, 4), 0.1), ((2, 4), 0.1),
                                  ((1, 4), 0.3), ((3, 4), 0.3),
                                  ((0, 2), 0.15), ((1, 3), 0.25)},
                                 2
                             )
                         ]
                         )
def test_measurement_error_probabilities(p_meas, expected_edges, repetitions):
    m = Matching(
        [[1, 1, 0], [0, 1, 1]],
        error_probabilities=[0.1, 0.2, 0.3],
        measurement_error_probabilities=p_meas,
        repetitions=repetitions
    )
    es = set((tuple(sorted([u, v])), d["error_probability"]) for u, v, d in m.edges())
    assert es == expected_edges

    # Check measurement_error_probability also accepted
    m = Matching(
        [[1, 1, 0], [0, 1, 1]],
        error_probabilities=[0.1, 0.2, 0.3],
        measurement_error_probability=p_meas,
        repetitions=repetitions
    )
    es = set((tuple(sorted([u, v])), d["error_probability"]) for u, v, d in m.edges())
    assert es == expected_edges


@pytest.mark.parametrize("m_errors", [[0.1, 0.01, 3], "A"])
def test_wrong_measurement_error_probabilities_raises_valueerror(m_errors):
    H = np.array([[1, 1, 0], [0, 1, 1]])
    with pytest.raises(ValueError):
        m = Matching()
        m.load_from_check_matrix(H, spacelike_weights=np.array([0.3, 0.7, 0.9]),
                                 measurement_error_probabilities=m_errors, repetitions=3)


def test_measurement_error_probabilities_and_probability_raises_value_error():
    H = np.array([[1, 1, 0], [0, 1, 1]])
    with pytest.raises(ValueError):
        m = Matching()
        m.load_from_check_matrix(H, spacelike_weights=np.array([0.3, 0.7, 0.9]),
                                 measurement_error_probabilities=[0.1, 0.1], repetitions=3,
                                 measurement_error_probability=[0.1, 0.1])


def test_set_boundary_nodes():
    m = Matching([[1, 1, 0], [0, 1, 1]])
    m.set_boundary_nodes({1, 2, 4})
    assert m.boundary == {1, 2, 4}
