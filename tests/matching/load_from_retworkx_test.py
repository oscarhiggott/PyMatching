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
import retworkx as rx
import pytest

from pymatching import Matching
from pymatching._cpp_pymatching import MatchingGraph


def test_boundary_from_retworkx():
    g = rx.PyGraph()
    g.add_nodes_from([{} for _ in range(5)])
    g.add_edge(4, 0, dict(fault_ids=0))
    g.add_edge(0, 1, dict(fault_ids=1))
    g.add_edge(1, 2, dict(fault_ids=2))
    g.add_edge(2, 3, dict(fault_ids=3))
    g.add_edge(3, 4, dict(fault_ids=4))
    g[4]['is_boundary'] = True
    m = Matching(g)
    assert m.boundary == {4}
    assert np.array_equal(m.decode(np.array([1, 0, 0, 0])), np.array([1, 0, 0, 0, 0]))
    assert np.array_equal(m.decode(np.array([0, 1, 0, 0])), np.array([1, 1, 0, 0, 0]))
    assert np.array_equal(m.decode(np.array([0, 1, 1, 0])), np.array([0, 0, 1, 0, 0]))
    assert np.array_equal(m.decode(np.array([0, 0, 1, 0])), np.array([0, 0, 0, 1, 1]))


def test_boundaries_from_retworkx():
    g = rx.PyGraph()
    g.add_nodes_from([{} for _ in range(6)])
    g.add_edge(0, 1, dict(fault_ids=0))
    g.add_edge(1, 2, dict(fault_ids=1))
    g.add_edge(2, 3, dict(fault_ids=2))
    g.add_edge(3, 4, dict(fault_ids=3))
    g.add_edge(4, 5, dict(fault_ids=4))
    g.add_edge(0, 5, dict(fault_ids=-1, weight=0.0))
    g.nodes()[0]['is_boundary'] = True
    g.nodes()[5]['is_boundary'] = True
    m = Matching(g)
    assert m.boundary == {0, 5}
    assert np.array_equal(m.decode(np.array([0, 1, 0, 0, 0, 0])), np.array([1, 0, 0, 0, 0]))
    assert np.array_equal(m.decode(np.array([0, 0, 1, 0, 0])), np.array([1, 1, 0, 0, 0]))
    assert np.array_equal(m.decode(np.array([0, 0, 1, 1, 0])), np.array([0, 0, 1, 0, 0]))
    assert np.array_equal(m.decode(np.array([0, 0, 0, 1, 0])), np.array([0, 0, 0, 1, 1]))


def test_unweighted_stabiliser_graph_from_retworkx():
    w = rx.PyGraph()
    w.add_nodes_from([{} for _ in range(6)])
    w.add_edge(0, 1, dict(fault_ids=0, weight=7.0))
    w.add_edge(0, 5, dict(fault_ids=1, weight=14.0))
    w.add_edge(0, 2, dict(fault_ids=2, weight=9.0))
    w.add_edge(1, 2, dict(fault_ids=-1, weight=10.0))
    w.add_edge(1, 3, dict(fault_ids=3, weight=15.0))
    w.add_edge(2, 5, dict(fault_ids=4, weight=2.0))
    w.add_edge(2, 3, dict(fault_ids=-1, weight=11.0))
    w.add_edge(3, 4, dict(fault_ids=5, weight=6.0))
    w.add_edge(4, 5, dict(fault_ids=6, weight=9.0))
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


def test_mwpm_from_retworkx():
    g = rx.PyGraph()
    g.add_nodes_from([{} for _ in range(3)])
    g.add_edge(0, 1, dict(fault_ids=0))
    g.add_edge(0, 2, dict(fault_ids=1))
    g.add_edge(1, 2, dict(fault_ids=2))
    m = Matching(g)
    assert (isinstance(m._matching_graph, MatchingGraph))
    assert (m.num_detectors == 3)
    assert (m.num_fault_ids == 3)

    g = rx.PyGraph()
    g.add_nodes_from([{} for _ in range(3)])
    g.add_edge(0, 1, {})
    g.add_edge(0, 2, {})
    g.add_edge(1, 2, {})
    m = Matching(g)
    assert (isinstance(m._matching_graph, MatchingGraph))
    assert (m.num_detectors == 3)
    assert (m.num_fault_ids == 0)

    g = rx.PyGraph()
    g.add_nodes_from([{} for _ in range(3)])
    g.add_edge(0, 1, dict(weight=1.5))
    g.add_edge(0, 2, dict(weight=1.7))
    g.add_edge(1, 2, dict(weight=1.2))
    m = Matching(g)
    assert (isinstance(m._matching_graph, MatchingGraph))
    assert (m.num_detectors == 3)
    assert (m.num_fault_ids == 0)


def test_matching_edges_from_retworkx():
    g = rx.PyGraph()
    g.add_nodes_from([{} for _ in range(4)])
    g.add_edge(0, 1, dict(fault_ids=0, weight=1.1, error_probability=0.1))
    g.add_edge(1, 2, dict(fault_ids=1, weight=2.1, error_probability=0.2))
    g.add_edge(2, 3, dict(fault_ids={2, 3}, weight=0.9, error_probability=0.3))
    g[0]['is_boundary'] = True
    g[3]['is_boundary'] = True
    g.add_edge(0, 3, dict(weight=0.0))
    m = Matching(g)
    es = list(m.edges())
    expected_edges = [
        (0, 1, {'fault_ids': {0}, 'weight': 1.1, 'error_probability': 0.1}),
        (1, 2, {'fault_ids': {1}, 'weight': 2.1, 'error_probability': 0.2}),
        (2, 3, {'fault_ids': {2, 3}, 'weight': 0.9, 'error_probability': 0.3}),
        (0, 3, {'fault_ids': set(), 'weight': 0.0, 'error_probability': -1.0}),
    ]
    print(es)
    assert es == expected_edges


def test_qubit_id_accepted_via_retworkx():
    g = rx.PyGraph()
    g.add_nodes_from([{} for _ in range(4)])
    g.add_edge(0, 1, dict(qubit_id=0, weight=1.1, error_probability=0.1))
    g.add_edge(1, 2, dict(qubit_id=1, weight=2.1, error_probability=0.2))
    g.add_edge(2, 3, dict(qubit_id={2, 3}, weight=0.9, error_probability=0.3))
    g[0]['is_boundary'] = True
    g[3]['is_boundary'] = True
    g.add_edge(0, 3, dict(weight=0.0))
    m = Matching(g)
    es = list(m.edges())
    expected_edges = [
        (0, 1, {'fault_ids': {0}, 'weight': 1.1, 'error_probability': 0.1}),
        (1, 2, {'fault_ids': {1}, 'weight': 2.1, 'error_probability': 0.2}),
        (2, 3, {'fault_ids': {2, 3}, 'weight': 0.9, 'error_probability': 0.3}),
        (0, 3, {'fault_ids': set(), 'weight': 0.0, 'error_probability': -1.0})
    ]
    assert es == expected_edges


def test_load_from_retworkx_raises_value_error_if_qubit_id_and_fault_ids_both_supplied():
    with pytest.raises(ValueError):
        g = rx.PyGraph()
        g.add_nodes_from([{} for _ in range(3)])
        g.add_edge(0, 1, dict(qubit_id=0, fault_ids=0))
        g.add_edge(1, 2, dict(qubit_id=1, fault_ids=1))
        m = Matching()
        m.load_from_retworkx(g)


def test_load_from_retworkx_type_errors_raised():
    with pytest.raises(TypeError):
        m = Matching()
        m.load_from_retworkx("A")
    with pytest.raises(TypeError):
        g = rx.PyGraph()
        g.add_nodes_from([{} for _ in range(2)])
        g.add_edge(0, 1, dict(fault_ids={0, "a"}))
        m = Matching()
        m.load_from_retworkx(g)
    with pytest.raises(TypeError):
        g = rx.PyGraph()
        g.add_nodes_from([{} for _ in range(2)])
        g.add_edge(0, 1, dict(fault_ids=[[0], [2]]))
        m = Matching()
        m.load_from_retworkx(g)
