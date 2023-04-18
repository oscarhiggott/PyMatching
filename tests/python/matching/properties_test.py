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

from pymatching.matching import Matching


def test_set_boundary_nodes():
    m = Matching([[1, 1, 0], [0, 1, 1]])
    m.set_boundary_nodes({1, 2, 4})
    assert m.boundary == {1, 2, 4}


def test_set_min_num_fault_ids():
    m = Matching()
    assert m.num_fault_ids == 0
    m.load_from_check_matrix([[1, 1, 0], [0, 1, 1]])
    assert m.num_fault_ids == 3
    m.ensure_num_fault_ids(5)
    assert m.num_fault_ids == 5
    assert m.decode([1, 0]).shape[0] == 5
    m.load_from_check_matrix([[1, 1, 0], [0, 1, 1]])
    assert m.num_fault_ids == 3

    g = nx.Graph()
    g.add_edge(0, 1, fault_ids={3})
    m = Matching()
    m.load_from_networkx(g)
    assert m.num_fault_ids == 4
    assert m.decode([1, 1]).shape[0] == 4
    m.load_from_networkx(g, min_num_fault_ids=7)
    assert m.num_fault_ids == 7
    assert m.decode([1, 1]).shape[0] == 7
    m.load_from_networkx(g, min_num_fault_ids=2)
    assert m.num_fault_ids == 4
    assert m.decode([1, 1]).shape[0] == 4

    g = rx.PyGraph()
    g.add_nodes_from([{} for _ in range(2)])
    g.add_edge(0, 1, dict(fault_ids=3))
    m.load_from_retworkx(g)
    assert m.num_fault_ids == 4
    assert m.decode([1, 1]).shape[0] == 4
    m.load_from_retworkx(g, min_num_fault_ids=7)
    assert m.num_fault_ids == 7
    assert m.decode([1, 1]).shape[0] == 7
    m.load_from_retworkx(g, min_num_fault_ids=2)
    assert m.num_fault_ids == 4
    assert m.decode([1, 1]).shape[0] == 4
