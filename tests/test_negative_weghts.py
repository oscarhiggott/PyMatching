import pytest
import numpy as np

from pymatching import Matching


@pytest.mark.parametrize("nn", (None, 30))
def test_negative_weight_repetition_code(nn):
    m = Matching()
    m.add_edge(0, 1, 0, -1)
    m.add_edge(1, 2, 1, -1)
    m.add_edge(2, 3, 2, -1)
    m.add_edge(3, 4, 3, -1)
    m.add_edge(4, 5, 4, -1)
    m.add_edge(5, 0, 5, -1)
    c, w = m.decode([0, 1, 1, 0, 0, 0], return_weight=True, num_neighbours=nn)
    assert np.array_equal(c, np.array([1, 0, 1, 1, 1, 1]))
    assert w == -5


@pytest.mark.parametrize("nn", (None, 30))
def test_isolated_negative_weight(nn):
    m = Matching()
    m.add_edge(0, 1, 0, 1)
    m.add_edge(1, 2, 1, -10)
    m.add_edge(2, 3, 2, 1)
    m.add_edge(3, 0, 3, 1)
    c, w = m.decode([0, 1, 1, 0], return_weight=True, num_neighbours=nn)
    assert np.array_equal(c, np.array([0, 1, 0, 0]))
    assert w == -10


@pytest.mark.parametrize("nn", (None, 30))
def test_negative_and_positive_in_matching(nn):
    m = Matching()
    m.add_edge(0, 1, 0, 1)
    m.add_edge(1, 2, 1, -10)
    m.add_edge(2, 3, 2, 1)
    m.add_edge(3, 0, 3, 1)
    c, w = m.decode([0, 1, 0, 1], return_weight=True, num_neighbours=nn)
    assert np.array_equal(c, np.array([0, 1, 1, 0]))
    assert w == -9


def test_negative_weight_edge_returned():
    m = Matching()
    m.add_edge(0, 1, weight=0.5, error_probability=0.3)
    m.add_edge(1, 2, weight=0.5, error_probability=0.3, fault_ids=0)
    m.add_edge(2, 3, weight=-0.5, error_probability=0.7, fault_ids={1, 2})
    expected = [(0, 1, {'fault_ids': set(), 'weight': 0.5, 'error_probability': 0.3}),
                (1, 2, {'fault_ids': {0}, 'weight': 0.5, 'error_probability': 0.3}),
                (2, 3, {'fault_ids': {1, 2}, 'weight': -0.5, 'error_probability': 0.7})]
    assert m.edges() == expected
