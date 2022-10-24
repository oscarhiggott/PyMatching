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

import networkx as nx
import retworkx as rx

from pymatching import Matching


def test_matching_to_networkx():
    g = nx.Graph()
    g.add_edge(0, 1, fault_ids={0}, weight=1.1, error_probability=0.1)
    g.add_edge(1, 2, fault_ids={1}, weight=2.1, error_probability=0.2)
    g.add_edge(2, 3, fault_ids={2, 3}, weight=0.9, error_probability=0.3)
    g.nodes[0]['is_boundary'] = True
    g.nodes[3]['is_boundary'] = True
    g.add_edge(0, 3, weight=0.0)
    m = Matching(g)

    g.edges[(0, 3)]['fault_ids'] = set()
    g.edges[(0, 3)]['error_probability'] = -1.0
    g.nodes[1]['is_boundary'] = False
    g.nodes[2]['is_boundary'] = False

    g2 = m.to_networkx()

    assert g.nodes(data=True) == g2.nodes(data=True)
    gedges = [({s, t}, d) for (s, t, d) in g.edges(data=True)]
    g2edges = [({s, t}, d) for (s, t, d) in g2.edges(data=True)]
    assert sorted(gedges) == sorted(g2edges)

    m = Matching()
    m.add_boundary_edge(0, weight=2)
    m.add_edge(0, 1, weight=3)
    m.add_edge(1, 2, weight=4)
    g = m.to_networkx()
    es = list(g.edges(data=True))
    assert es == [(0, 3, {"weight": 2.0, "error_probability": -1, "fault_ids": set()}),
                  (0, 1, {"weight": 3.0, "error_probability": -1, "fault_ids": set()}),
                  (1, 2, {"weight": 4.0, "error_probability": -1, "fault_ids": set()})]
    assert sorted(list(g.nodes(data=True))) == [(0, {"is_boundary": False}), (1, {"is_boundary": False}),
                                                (2, {"is_boundary": False}), (3, {"is_boundary": True})]

    m = Matching()
    m.add_edge(0, 1)
    g = m.to_networkx()
    assert list(g.edges(data=True)) == [(0, 1, {"weight": 1.0, "error_probability": -1, "fault_ids": set()})]
    assert list(g.nodes(data=True)) == [(0, {"is_boundary": False}), (1, {"is_boundary": False})]


def test_matching_to_retworkx():
    g = rx.PyGraph()
    g.add_nodes_from([{} for _ in range(4)])
    g.add_edge(0, 1, dict(fault_ids={0}, weight=1.1, error_probability=0.1))
    g.add_edge(1, 2, dict(fault_ids={1}, weight=2.1, error_probability=0.2))
    g.add_edge(2, 3, dict(fault_ids={2, 3}, weight=0.9, error_probability=0.3))
    g[0]['is_boundary'] = True
    g[3]['is_boundary'] = True
    g.add_edge(0, 3, dict(weight=0.0))
    m = Matching(g)

    edge_0_3 = g.get_edge_data(0, 3)
    edge_0_3['fault_ids'] = set()
    edge_0_3['error_probability'] = -1.0
    g[1]['is_boundary'] = False
    g[2]['is_boundary'] = False

    g2 = m.to_retworkx()

    assert g.node_indices() == g2.node_indices()
    gedges = [({s, t}, d) for (s, t, d) in g.weighted_edge_list()]
    g2edges = [({s, t}, d) for (s, t, d) in g.weighted_edge_list()]
    assert sorted(gedges) == sorted(g2edges)

    m = Matching()
    m.add_boundary_edge(0, weight=2)
    m.add_edge(0, 1, weight=3)
    m.add_edge(1, 2, weight=4)
    g = m.to_retworkx()
    es = list(g.weighted_edge_list())
    assert es == [(0, 3, {"weight": 2.0, "error_probability": -1, "fault_ids": set()}),
                  (0, 1, {"weight": 3.0, "error_probability": -1, "fault_ids": set()}),
                  (1, 2, {"weight": 4.0, "error_probability": -1, "fault_ids": set()})]
    assert list(g.nodes()) == [{"is_boundary": False}, {"is_boundary": False},
                               {"is_boundary": False}, {"is_boundary": True}]

    m = Matching()
    m.add_edge(0, 1)
    g = m.to_retworkx()
    assert list(g.weighted_edge_list()) == [(0, 1, {"weight": 1.0, "error_probability": -1, "fault_ids": set()})]
    assert list(g.nodes()) == [{"is_boundary": False}, {"is_boundary": False}]


def test_negative_weight_edge_returned():
    m = Matching()
    m.add_edge(0, 1, weight=0.5, error_probability=0.3)
    m.add_edge(1, 2, weight=0.5, error_probability=0.3, fault_ids=0)
    m.add_edge(2, 3, weight=-0.5, error_probability=0.7, fault_ids={1, 2})
    expected = [(0, 1, {'fault_ids': set(), 'weight': 0.5, 'error_probability': 0.3}),
                (1, 2, {'fault_ids': {0}, 'weight': 0.5, 'error_probability': 0.3}),
                (2, 3, {'fault_ids': {1, 2}, 'weight': -0.5, 'error_probability': 0.7})]
    assert m.edges() == expected


def test_self_loop_to_networkx():
    m = Matching()
    m.add_edge(0, 0, weight=3)
    m.add_edge(1, 1, weight=2)
    m.add_edge(1, 2, weight=1)
    g = m.to_networkx()
    assert list(g.edges(data=True)) == [(0, 0, {'fault_ids': set(), 'weight': 3.0, 'error_probability': -1.0}),
                                        (1, 1, {'fault_ids': set(), 'weight': 2.0, 'error_probability': -1.0}),
                                        (1, 2, {'fault_ids': set(), 'weight': 1.0, 'error_probability': -1.0})]
