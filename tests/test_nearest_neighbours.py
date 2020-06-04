import networkx as nx

from pymatching import Matching


def test_nearest_neighbour_nodes():
    g = nx.Graph()
    g.add_edge(0,1, qubit_id=0, weight=0.1)
    g.add_edge(1,2, qubit_id=1, weight=0.1)
    g.add_edge(2,3, qubit_id=2, weight=0.1)
    g.add_edge(3,4, qubit_id=3, weight=0.1)
    g.add_edge(4,5, qubit_id=4, weight=0.1)
    m = Matching(g)
    assert set(m.stabiliser_graph.nearest_bfs_neighbours(2, 3)) == {1,2,3}
    assert set(m.stabiliser_graph.nearest_bfs_neighbours(2, 5)) == {0,1,2,3,4}
    assert set(m.stabiliser_graph.nearest_bfs_neighbours(3, 1)) == {3}
    assert set(m.stabiliser_graph.nearest_bfs_neighbours(2, 10)) == {0,1,2,3,4,5}
    assert set(m.stabiliser_graph.nearest_bfs_neighbours(4, 0)) == set()
