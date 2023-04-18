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
import matplotlib.pyplot as plt

from pymatching import Matching


def test_draw_matching():
    g = nx.Graph()
    g.add_edge(0, 1, fault_ids={0}, weight=1.1, error_probability=0.1)
    g.add_edge(1, 2, fault_ids={1}, weight=2.1, error_probability=0.2)
    g.add_edge(2, 3, fault_ids={2, 3}, weight=0.9, error_probability=0.3)
    g.nodes[0]['is_boundary'] = True
    g.nodes[3]['is_boundary'] = True
    g.add_edge(0, 3, weight=0.0)
    m = Matching(g)
    plt.figure()
    m.draw()
