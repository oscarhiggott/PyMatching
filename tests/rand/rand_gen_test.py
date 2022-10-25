from scipy.sparse import csr_matrix
import numpy as np
from pymatching import Matching, rand_float


def test_add_noise_without_error_probabilities_returns_none():
    m = Matching(csr_matrix(np.array([[1, 1, 0], [0, 1, 1]])))
    assert m.add_noise() is None
    m = Matching(csr_matrix(np.array([[1, 1, 0], [0, 1, 1]])),
                 error_probabilities=np.array([0.5, 0.7, -0.1]))
    assert m.add_noise() is None


def test_rand_float():
    N = 1000
    s = sum(rand_float(0., 1.) for i in range(N))
    assert 430 < s < 570
