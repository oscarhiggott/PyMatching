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

from pymatching import Matching


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
