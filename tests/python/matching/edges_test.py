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

import pytest

from pymatching.matching import Matching


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


def test_add_edge_merge_strategy():
    m = Matching()
    m.add_edge(0, 10, fault_ids={0}, weight=1.2, error_probability=0.3)
    assert m.edges() == [(0, 10, {'fault_ids': {0}, 'weight': 1.2, 'error_probability': 0.3})]
    m.add_edge(0, 10, fault_ids={1}, weight=1.0, error_probability=0.35, merge_strategy="smallest-weight")
    assert m.edges() == [(0, 10, {'fault_ids': {1}, 'weight': 1.0, 'error_probability': 0.35})]
    with pytest.raises(ValueError):
        m.add_edge(0, 10, fault_ids={1}, weight=1.5, error_probability=0.6, merge_strategy="disallow")
    m.add_edge(0, 10, fault_ids={2}, weight=4.0, error_probability=0.2, merge_strategy="independent")
    es = m.edges()
    es[0][2]["weight"] = round(es[0][2]["weight"], 6)
    assert es == [(0, 10, {'fault_ids': {1}, 'weight': 0.958128, 'error_probability': 0.41})]
    m = Matching()
    m.add_edge(1, 10, fault_ids={0}, weight=2, error_probability=0.3)
    m.add_edge(1, 10, fault_ids={1}, weight=5, error_probability=0.1, merge_strategy="keep-original")
    assert m.edges() == [(1, 10, {'fault_ids': {0}, 'weight': 2, 'error_probability': 0.3})]
    m.add_edge(1, 10, fault_ids={2}, weight=5, error_probability=0.1, merge_strategy="replace")
    assert m.edges() == [(1, 10, {'fault_ids': {2}, 'weight': 5, 'error_probability': 0.1})]


def test_add_boundary_edge():
    m = Matching()
    m.add_boundary_edge(0, fault_ids={0}, weight=1.2, error_probability=0.3)
    assert m.edges() == [(0, None, {'fault_ids': {0}, 'weight': 1.2, 'error_probability': 0.3})]
    m.add_boundary_edge(0, fault_ids={1}, weight=1.0, error_probability=0.35, merge_strategy="smallest-weight")
    assert m.edges() == [(0, None, {'fault_ids': {1}, 'weight': 1.0, 'error_probability': 0.35})]
    with pytest.raises(ValueError):
        m.add_boundary_edge(0, fault_ids={1}, weight=1.5, error_probability=0.6, merge_strategy="disallow")
    m.add_boundary_edge(0, fault_ids={2}, weight=4.0, error_probability=0.2, merge_strategy="independent")
    es = m.edges()
    es[0][2]["weight"] = round(es[0][2]["weight"], 6)
    assert es == [(0, None, {'fault_ids': {1}, 'weight': 0.958128, 'error_probability': 0.41})]
    m = Matching()
    m.add_boundary_edge(1, fault_ids={0}, weight=2, error_probability=0.3)
    m.add_boundary_edge(1, fault_ids={1}, weight=5, error_probability=0.1, merge_strategy="keep-original")
    assert m.edges() == [(1, None, {'fault_ids': {0}, 'weight': 2, 'error_probability': 0.3})]
    m.add_boundary_edge(1, fault_ids={2}, weight=5, error_probability=0.1, merge_strategy="replace")
    assert m.edges() == [(1, None, {'fault_ids': {2}, 'weight': 5, 'error_probability': 0.1})]


def test_has_edge():
    m = Matching()
    m.add_edge(0, 1)
    m.add_edge(1, 2)
    m.add_boundary_edge(2)
    m.add_boundary_edge(5)
    assert m.has_edge(0, 1)
    assert not m.has_edge(7, 8)
    assert m.has_edge(1, 2)
    assert m.has_boundary_edge(2)
    assert m.has_boundary_edge(5)
    assert not m.has_boundary_edge(0)
    assert not m.has_boundary_edge(100)


def test_get_edge_data():
    m = Matching()
    m.add_edge(0, 1, {0, 1, 2}, 2.5, 0.1)
    assert m.get_edge_data(0, 1) == {"fault_ids": {0, 1, 2}, "weight": 2.5, "error_probability": 0.1}
    m.add_boundary_edge(5, {5, 6, 7}, 5.0, 0.34)
    assert m.get_boundary_edge_data(5) == {"fault_ids": {5, 6, 7}, "weight": 5.0, "error_probability": 0.34}


def test_large_edge_weight_not_added_to_graph():
    m = Matching()
    m.add_edge(0, 1)
    with pytest.warns(UserWarning):
        m.add_edge(1, 2, weight=9999999999)
    with pytest.warns(UserWarning):
        m.add_boundary_edge(3, weight=9999999999)
    with pytest.warns(UserWarning):
        m.add_edge(3, 4, weight=-9999999999)
    with pytest.warns(UserWarning):
        m.add_boundary_edge(5, weight=-9999999999)
    assert m.num_edges == 1
    assert m.edges() == [(0, 1, {"fault_ids": set(), "weight": 1.0, "error_probability": -1.0})]


def test_add_self_loop():
    m = Matching()
    m.add_edge(0, 1, weight=2)
    m.add_edge(2, 2, weight=5)
    m.add_edge(3, 4, weight=10)
    m.add_edge(4, 4, weight=11)
    assert m.edges() == [
        (0, 1, {"fault_ids": set(), "weight": 2.0, "error_probability": -1.0}),
        (2, 2, {"fault_ids": set(), "weight": 5.0, "error_probability": -1.0}),
        (3, 4, {"fault_ids": set(), "weight": 10.0, "error_probability": -1.0}),
        (4, 4, {"fault_ids": set(), "weight": 11.0, "error_probability": -1.0})
    ]
    with pytest.raises(ValueError):
        m.add_edge(4, 4, weight=12)
    m.add_edge(4, 4, weight=10, merge_strategy="smallest-weight")
    m.add_edge(4, 3, weight=14, merge_strategy="replace")
    assert m.edges() == [
        (0, 1, {"fault_ids": set(), "weight": 2.0, "error_probability": -1.0}),
        (2, 2, {"fault_ids": set(), "weight": 5.0, "error_probability": -1.0}),
        (3, 4, {"fault_ids": set(), "weight": 14.0, "error_probability": -1.0}),
        (4, 4, {"fault_ids": set(), "weight": 10.0, "error_probability": -1.0})
    ]
