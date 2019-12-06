import numpy as np
from scipy.sparse import csc_matrix

from quantumcode import load_css_from_hdf5
from mwpm._cpp_mwpm import (_stabiliser_graph, breadth_first_search, 
                            all_pairs_shortest_path, shortest_path,
                            _decode)
from mwpm import MWPM


def test_cpp_syndrome_graph():
    fn = "css_toric_[[18,2,3]]_rank_deficient"
    css = load_css_from_hdf5(fn)
    H = csc_matrix(css.x_stabilisers)
    assert(np.count_nonzero(H.indptr[1:]-H.indptr[0:-1]-2) == 0)
    gd = _stabiliser_graph(H.indices, H.shape[0])
    expected = [[1, 2, 3, 6], [0, 2, 4, 7], [1, 0, 5, 8], 
                [0, 4, 5, 6], [1, 3, 5, 7], [2, 4, 3, 8], 
                [3, 7, 8, 0], [4, 6, 8, 1], [5, 7, 6, 2]]
    assert(gd.g==expected)


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


def test_cpp_decode():
    fn = "css_toric_[[18,2,3]]_rank_deficient"
    css = load_css_from_hdf5(fn)
    H = csc_matrix(css.x_stabilisers)
    gd = _stabiliser_graph(H.indices, H.shape[0])
    apsp = all_pairs_shortest_path(gd.g)
    n = np.zeros(H.shape[1], dtype=int)
    n[5] = 1
    n[10] = 1
    z = H.dot(n)
    defects = np.nonzero(z)[0]
    c = _decode(apsp, defects, gd.qubit, H.shape[1])
    assert(np.array_equal(c,n))


def test_mwpm_decode_method():
    fn = "css_toric_[[18,2,3]]_rank_deficient"
    H = load_css_from_hdf5(fn).x_stabilisers
    m = MWPM(H)
    n = np.zeros(H.shape[1], dtype=int)
    n[5] = 1
    n[10] = 1
    z = H.dot(n)
    c = m.decode(z)
    assert(np.array_equal(c,n))
