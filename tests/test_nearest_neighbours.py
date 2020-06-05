import networkx as nx
import numpy as np
import pytest

from pymatching import Matching
from pymatching._cpp_mwpm import decode_match_neighbourhood

g = nx.Graph()
g.add_edge(0,1, qubit_id=0, weight=0.1)
g.add_edge(1,2, qubit_id=1, weight=0.1)
g.add_edge(2,3, qubit_id=2, weight=0.1)
g.add_edge(3,4, qubit_id=3, weight=0.1)
g.add_edge(4,5, qubit_id=4, weight=0.1)
m = Matching(g)


def test_dijkstra_nearest_neighbour_nodes():
    assert (set(m.stabiliser_graph.get_nearest_neighbours(2, 3, [-1,0,1,-1,2,-1]))
             == {(2, 0.0), (1, 0.1), (4, 0.2)})


def test_dijkstra_path():
    assert m.stabiliser_graph.get_path(1, 4) == [1,2,3,4]
    assert m.stabiliser_graph.get_path(4, 1) == [4,3,2,1]
    assert m.stabiliser_graph.get_path(5, 3) == [5,4,3]

neighbour_match_fixtures = [
    ([1,3], 3, [0,1,1,0,0]),
    ([0,1,2,3,4,5], 2, [1,0,1,0,1]),
    ([0,1,2,3,4,5], 10, [1,0,1,0,1])
]
@pytest.mark.parametrize("defects,num_neighbours,correction", neighbour_match_fixtures)
def test_neighbourhood_matching(defects,num_neighbours,correction):
    assert (np.array_equal(decode_match_neighbourhood(m.stabiliser_graph, np.array(defects), 
            num_neighbours),
             np.array(correction)))

