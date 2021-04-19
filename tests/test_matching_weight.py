import numpy as np
from scipy.sparse import csr_matrix
import pytest

from pymatching import Matching


def repetition_code(n):
    row_ind, col_ind = zip(*((i, j) for i in range(n) for j in (i, (i+1)%n)))
    data = np.ones(2*n, dtype=np.uint8)
    return csr_matrix((data, (row_ind, col_ind)))


weight_fixtures = [
    (10, 20),
    (10, None),
    (15, 10),
    (15, None),
    (20, 1),
    (20, None)
]


@pytest.mark.parametrize("n,num_neighbours", weight_fixtures)
def test_matching_weight(n, num_neighbours):
    p = 0.4
    H = repetition_code(n)
    noise = np.random.rand(n) < p
    weights = np.random.rand(n)
    s = H@noise % 2
    m = Matching(H, spacelike_weights=weights)
    corr, weight = m.decode(s, num_neighbours=num_neighbours, return_weight=True)
    expected_weight = np.sum(weights[corr==1])
    assert expected_weight == pytest.approx(weight)
