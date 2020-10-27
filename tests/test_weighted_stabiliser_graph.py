import pytest
from pymatching._cpp_mwpm import decode, WeightedStabiliserGraph


def test_weighted_spacetime_shortest_path():
    w = WeightedStabiliserGraph(6, [])
    w.add_edge(0, 1, {0}, 7.0)
    w.add_edge(0, 5, {1}, 14.0)
    w.add_edge(0, 2, {2}, 9.0)
    w.add_edge(1, 2, {3}, 10.0)
    w.add_edge(1, 3, {4}, 15.0)
    w.add_edge(2, 5, {5}, 2.0)
    w.add_edge(2, 3, {6}, 11.0)
    w.add_edge(3, 4, {7}, 6.0)
    w.add_edge(4, 5, {8}, 9.0)
    w.compute_all_pairs_shortest_paths()

    assert(w.qubit_ids(3, 1) == {4})
    assert(w.distance(1, 2) == pytest.approx(10.0))
    assert(w.distance(5, 0) == pytest.approx(11.0))
    assert(w.shortest_path(3, 5) == [3, 2, 5])
    assert(w.shortest_path(4, 2) == [4, 5, 2])
    assert(w.get_num_qubits() == 9)
    assert(w.get_num_stabilisers() == 6)


def test_weighted_num_qubits_and_stabilisers():
    w = WeightedStabiliserGraph(6, [])
    w.add_edge(0, 1, {0}, 7.0)
    w.add_edge(0, 5, {1}, 14.0)
    w.add_edge(0, 2, {2}, 9.0)
    w.add_edge(1, 2, set(), 10.0)
    w.add_edge(1, 3, {3}, 15.0)
    w.add_edge(2, 5, {4}, 2.0)
    w.add_edge(2, 3, set(), 11.0)
    w.add_edge(3, 4, {5}, 6.0)
    w.add_edge(4, 5, {6}, 9.0)
    w.compute_all_pairs_shortest_paths()
    assert(w.get_num_qubits() == 7)
    assert(w.get_num_stabilisers() == 6)


@pytest.mark.parametrize("num_loops", range(1, 10))
def test_num_connected_components(num_loops):
    loop_size = 3
    w = WeightedStabiliserGraph(num_loops*loop_size, [])
    q = 0
    for i in range(num_loops):
        for j in range(loop_size):
            n1 = j + i * loop_size
            n2 = ((j + 1) % loop_size) + i * loop_size
            w.add_edge(n1, n2, {q}, 1.0)
            q += 1
    assert w.get_num_connected_components() == num_loops

