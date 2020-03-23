import os

import pytest
import numpy as np
from scipy.sparse import csc_matrix, load_npz
import pytest

from mwpm._cpp_mwpm import (breadth_first_search, 
                            all_pairs_shortest_path, shortest_path,
                            decode, UnweightedStabiliserGraph,
                            WeightedStabiliserGraph)
from mwpm import MWPM

TEST_DIR = dir_path = os.path.dirname(os.path.realpath(__file__))


def test_cpp_syndrome_graph():
    fn = "css_toric_[[18,2,3]]_rank_deficient_Hx.npz"
    H = csc_matrix(load_npz(os.path.join(TEST_DIR, 'data', fn)))
    assert(np.count_nonzero(H.indptr[1:]-H.indptr[0:-1]-2) == 0)
    sg = UnweightedStabiliserGraph(H.indices)
    expected = [[1, 2, 3, 6], [0, 2, 4, 7], [1, 0, 5, 8], 
                [0, 4, 5, 6], [1, 3, 5, 7], [2, 4, 3, 8], 
                [3, 7, 8, 0], [4, 6, 8, 1], [5, 7, 6, 2]]
    assert(sg.adj_list==expected)


def test_breadth_first_search():
    g = [[1,10],[0,2],[1,3,5],[2,4],[3],[2],[7,10],
         [6,8],[7,9],[8,10],[9,0,6,11],[10]]
    res = breadth_first_search(g, 0)
    dist_exp = [0,1,2,3,4,3,2,3,3,2,1,2]
    parent_exp = [-1,0,1,2,3,2,10,6,9,10,0,10]
    assert(res.distance == dist_exp)
    assert(res.parent == parent_exp)


def test_all_pairs_shortest_path():
    g = [[1],[0,2],[1]]
    res = all_pairs_shortest_path(g)
    dists_exp = [[0,1,2],[1,0,1],[2,1,0]]
    parents_exp = [[-1,0,1],[1,-1,1],[1,2,-1]]
    assert(res.distances == dists_exp)
    assert(res.parents == parents_exp)


def test_shortest_path():
    g = [[1,10],[0,2],[1,3,5],[2,4],[3],[2],[7,10],
         [6,8],[7,9],[8,10],[9,0,6,11],[10]]
    res = breadth_first_search(g, 0)
    path = shortest_path(res.parent, 7)
    expected_path = [7,6,10,0]
    assert(path == expected_path)


def test_mwpm_decode_method():
    fn = "css_toric_[[18,2,3]]_rank_deficient_Hx.npz"
    H = load_npz(os.path.join(TEST_DIR, 'data', fn))
    m = MWPM(H)
    n = np.zeros(H.shape[1], dtype=int)
    n[5] = 1
    n[10] = 1
    z = H.dot(n) % 2
    c = m.decode(z)
    assert(np.array_equal(c,n))


noisy_fixtures = [
    (
        np.array([
            [0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        ]),
        np.array([
            [0,0,0,0,0,0,0,1,0],
            [0,0,0,0,0,0,1,0,0],
            [0,0,0,0,0,0,0,0,0]
        ]).T,
        np.array([0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0])    
    ),
    (
        np.array([
            [0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        ]),
        np.array([
            [0,0,0,0,0,1,0,0,0],
            [0,0,0,0,0,0,0,0,0],
            [0,0,0,0,0,0,0,0,0]
        ]).T,
        np.array([0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0])
    )
]

@pytest.mark.parametrize("n,z_err,c_expected", noisy_fixtures)
def test_mwpm_noisy_decode(n, z_err, c_expected):
    fn = "css_toric_[[18,2,3]]_rank_deficient_Hx.npz"
    H = load_npz(os.path.join(TEST_DIR, 'data', fn))
    m = MWPM(H)
    n_all = np.cumsum(n, 0) % 2
    z_noiseless = H.dot(n_all.T) % 2
    z_noisy = (z_noiseless + z_err) % 2
    z_noisy[:,1:] = (z_noisy[:,1:] - z_noisy[:,:-1]) % 2
    c = m.decode(z_noisy)
    assert(np.array_equal(c, c_expected))


distance_fixtures = [
    (2,11,1),
    (3,13,2),
    (2,8,1),
    (2,98, 11)
]


@pytest.mark.parametrize("node1,node2,expected", distance_fixtures)
def test_spacetime_distance(node1, node2, expected):
    fn = "css_toric_[[18,2,3]]_rank_deficient_Hx.npz"
    H = load_npz(os.path.join(TEST_DIR, 'data', fn))
    m = MWPM(H)
    d = m.stabiliser_graph.space_time_distance(node1, node2)
    assert(d == expected)


spacetime_path_fixtures = [
    (1,184,[1,4]),
    (2,62,[2,8])
]


@pytest.mark.parametrize("node1,node2,expected", spacetime_path_fixtures)
def test_spacetime_shortest_path(node1, node2, expected):
    fn = "css_toric_[[18,2,3]]_rank_deficient_Hx.npz"
    H = load_npz(os.path.join(TEST_DIR, 'data', fn))
    m = MWPM(H)
    path = m.stabiliser_graph.space_time_shortest_path(node1, node2)
    assert(path == expected)


def test_weighted_spacetime_shortest_path():
    w = WeightedStabiliserGraph(6)
    w.add_edge(0, 1, 0, 7.0)
    w.add_edge(0, 5, 1, 14.0)
    w.add_edge(0, 2, 2, 9.0)
    w.add_edge(1, 2, 3, 10.0)
    w.add_edge(1, 3, 4, 15.0)
    w.add_edge(2, 5, 5, 2.0)
    w.add_edge(2, 3, 6, 11.0)
    w.add_edge(3, 4, 7, 6.0)
    w.add_edge(4, 5, 8, 9.0)
    w.compute_all_pairs_shortest_paths()

    assert(w.qubit_id(3, 1) == 4)
    assert(w.distance(1, 2) == pytest.approx(10.0))
    assert(w.distance(5, 0) == pytest.approx(11.0))
    assert(w.shortest_path(3, 5) == [3, 2, 5])
    assert(w.shortest_path(4, 2) == [4, 5, 2])
    assert(w.get_num_qubits() == 9)
    assert(w.get_num_stabilisers() == 6)
