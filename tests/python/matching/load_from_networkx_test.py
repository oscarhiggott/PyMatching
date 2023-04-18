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

import numpy as np
import pytest
import networkx as nx

from pymatching import Matching
from pymatching._cpp_pymatching import MatchingGraph


def test_bad_fault_ids_raises_value_error():
    g = nx.Graph()
    g.add_edge(0, 1, fault_ids='test')
    with pytest.raises(TypeError):
        Matching(g)
    g = nx.Graph()
    g.add_edge(0, 1, fault_ids=[[1], [2]])
    with pytest.raises(TypeError):
        Matching(g)


def test_boundary_from_networkx():
    g = nx.Graph()
    g.add_edge(4, 0, fault_ids=0)
    g.add_edge(0, 1, fault_ids=1)
    g.add_edge(1, 2, fault_ids=2)
    g.add_edge(2, 3, fault_ids=3)
    g.add_edge(3, 4, fault_ids=4)
    g.nodes()[4]['is_boundary'] = True
    m = Matching.from_networkx(g)
    assert m.boundary == {4}
    assert np.array_equal(m.decode(np.array([1, 0, 0, 0])), np.array([1, 0, 0, 0, 0]))
    assert np.array_equal(m.decode(np.array([0, 1, 0, 0])), np.array([1, 1, 0, 0, 0]))
    assert np.array_equal(m.decode(np.array([0, 1, 1, 0])), np.array([0, 0, 1, 0, 0]))
    assert np.array_equal(m.decode(np.array([0, 0, 1, 0])), np.array([0, 0, 0, 1, 1]))


def test_boundaries_from_networkx():
    g = nx.Graph()
    g.add_edge(0, 1, fault_ids=0)
    g.add_edge(1, 2, fault_ids=1)
    g.add_edge(2, 3, fault_ids=2)
    g.add_edge(3, 4, fault_ids=3)
    g.add_edge(4, 5, fault_ids=4)
    g.add_edge(0, 5, fault_ids=-1, weight=0.0)
    g.nodes()[0]['is_boundary'] = True
    g.nodes()[5]['is_boundary'] = True
    m = Matching.from_networkx(g)
    assert m.boundary == {0, 5}
    assert np.array_equal(m.decode(np.array([0, 1, 0, 0, 0, 0])), np.array([1, 0, 0, 0, 0]))
    assert np.array_equal(m.decode(np.array([0, 0, 1, 0, 0])), np.array([1, 1, 0, 0, 0]))
    assert np.array_equal(m.decode(np.array([0, 0, 1, 1, 0])), np.array([0, 0, 1, 0, 0]))
    assert np.array_equal(m.decode(np.array([0, 0, 0, 1, 0])), np.array([0, 0, 0, 1, 1]))


def test_wrong_networkx_graph_type_raises_type_error():
    m = Matching()
    with pytest.raises(TypeError):
        m.load_from_networkx("test")


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
    assert (m.num_fault_ids == 7)
    assert (m.num_detectors == 6)
    assert (np.array_equal(
        m.decode(np.array([1, 0, 1, 0, 0, 0])),
        np.array([0, 0, 1, 0, 0, 0, 0]))
    )
    with pytest.raises(ValueError):
        m.decode(np.array([1, 1, 0]))
    with pytest.raises(ValueError):
        m.decode(np.array([1, 1, 1, 0, 0, 0]))
    assert (np.array_equal(
        m.decode(np.array([1, 0, 0, 0, 0, 1])),
        np.array([0, 0, 1, 0, 1, 0, 0]))
    )
    assert (np.array_equal(
        m.decode(np.array([0, 1, 0, 0, 0, 1])),
        np.array([0, 0, 0, 0, 1, 0, 0]))
    )


def test_mwpm_from_networkx():
    g = nx.Graph()
    g.add_edge(0, 1, fault_ids=0)
    g.add_edge(0, 2, fault_ids=1)
    g.add_edge(1, 2, fault_ids=2)
    m = Matching(g)
    assert (isinstance(m._matching_graph, MatchingGraph))
    assert (m.num_detectors == 3)
    assert (m.num_fault_ids == 3)

    g = nx.Graph()
    g.add_edge(0, 1)
    g.add_edge(0, 2)
    g.add_edge(1, 2)
    m = Matching(g)
    assert (isinstance(m._matching_graph, MatchingGraph))
    assert (m.num_detectors == 3)
    assert (m.num_fault_ids == 0)

    g = nx.Graph()
    g.add_edge(0, 1, weight=1.5)
    g.add_edge(0, 2, weight=1.7)
    g.add_edge(1, 2, weight=1.2)
    m = Matching(g)
    assert (isinstance(m._matching_graph, MatchingGraph))
    assert (m.num_detectors == 3)
    assert (m.num_fault_ids == 0)


def test_matching_edges_from_networkx():
    g = nx.Graph()
    g.add_edge(0, 1, fault_ids=0, weight=1.1, error_probability=0.1)
    g.add_edge(1, 2, fault_ids=1, weight=2.1, error_probability=0.2)
    g.add_edge(2, 3, fault_ids={2, 3}, weight=0.9, error_probability=0.3)
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


def test_load_from_networkx_raises_value_error_if_qubit_id_and_fault_ids_both_supplied():
    with pytest.raises(ValueError):
        g = nx.Graph()
        g.add_edge(0, 1, qubit_id=0, fault_ids=0)
        g.add_edge(1, 2, qubit_id=1, fault_ids=1)
        m = Matching()
        m.load_from_networkx(g)
