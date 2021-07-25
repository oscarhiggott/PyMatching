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

import sys

import networkx as nx
import numpy as np
import pytest

from pymatching import Matching
from pymatching._cpp_mwpm import local_matching


g = nx.Graph()
g.add_edge(0,1, qubit_id=0, weight=0.3)
g.add_edge(1,2, qubit_id=1, weight=0.1)
g.add_edge(2,3, qubit_id=2, weight=0.1)
g.add_edge(3,4, qubit_id=3, weight=0.1)
g.add_edge(4,5, qubit_id=4, weight=0.1)
m = Matching(g)


def test_dijkstra_nearest_neighbour_nodes():
    assert (set(m.matching_graph.get_nearest_neighbours(2, 3, [0,0,1,-1,2,-1]))
            == {(2, 0.0), (1, 0.1), (4, 0.2)})


def test_dijkstra_path():
    assert m.matching_graph.get_path(1, 4) == [1,2,3,4]
    assert m.matching_graph.get_path(4, 1) == [4,3,2,1]
    assert m.matching_graph.get_path(5, 3) == [5,4,3]


neighbour_match_fixtures = [
    ([1, 3], 3, [0, 1, 1, 0, 0]),
    ([0, 1, 2, 3, 4, 5], 3, [1, 0, 1, 0, 1]),
    ([0, 1, 2, 3, 4, 5], 11, [1, 0, 1, 0, 1])
]


@pytest.mark.parametrize("defects,num_neighbours,correction", neighbour_match_fixtures)
def test_neighbourhood_matching(defects,num_neighbours,correction):
    assert (np.array_equal(local_matching(m.matching_graph,
            np.array(defects), num_neighbours, False).correction, np.array(correction)))
