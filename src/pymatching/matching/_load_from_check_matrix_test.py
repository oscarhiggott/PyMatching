# Copyright 2022 Oscar Higgott

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#      http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import numpy as np
from scipy.sparse import csr_matrix, csc_matrix
import pytest

from pymatching import Matching


def test_boundary_from_check_matrix():
    H = csr_matrix(np.array([[1, 1, 0, 0, 0], [0, 1, 1, 0, 0],
                             [0, 0, 1, 1, 0], [0, 0, 0, 1, 1]]))
    m = Matching(H=H)   # Checks `H` still accepted as keyword argument despite being deprecated
    assert m.boundary == {4}
    assert np.array_equal(m.decode(np.array([1, 0, 0, 0])), np.array([1, 0, 0, 0, 0]))
    assert np.array_equal(m.decode(np.array([0, 1, 0, 0])), np.array([1, 1, 0, 0, 0]))
    assert np.array_equal(m.decode(np.array([0, 1, 1, 0])), np.array([0, 0, 1, 0, 0]))
    assert np.array_equal(m.decode(np.array([0, 0, 1, 0])), np.array([0, 0, 0, 1, 1]))


def test_nonzero_matrix_elements_not_one_raises_value_error():
    H = csr_matrix(np.array([[0, 1.01, 1.01], [1.01, 1.01, 0]]))
    with pytest.raises(ValueError):
        Matching(H)


def test_too_many_checks_per_qubit_raises_value_error():
    H = csr_matrix(np.array([[1, 1, 0, 0], [1, 0, 1, 0], [1, 0, 0, 1]]))
    with pytest.raises(ValueError):
        Matching(H)


def test_wrong_check_matrix_type_raises_type_error():
    with pytest.raises(TypeError):
        Matching("test")
    m = Matching()
    with pytest.raises(TypeError):
        m.load_from_check_matrix("test")


def test_error_probability_from_array():
    H = csr_matrix(np.array([[1, 1, 0, 0, 0], [0, 1, 1, 0, 0],
                             [0, 0, 1, 1, 0], [0, 0, 0, 1, 1]]))
    m = Matching(H, error_probabilities=np.array([0., 0., 0., 0., 1.]))
    assert np.array_equal(m.add_noise()[0], np.array([0, 0, 0, 0, 1]))
    assert np.array_equal(m.add_noise()[1], np.array([0, 0, 0, 1, 0]))
    m = Matching(H, error_probabilities=np.array([0., 0., 0., 0., 0.]))
    assert np.array_equal(m.add_noise()[0], np.array([0, 0, 0, 0, 0]))
    assert np.array_equal(m.add_noise()[1], np.array([0, 0, 0, 0, 0]))
    m = Matching(H, error_probabilities=0.0)
    assert np.array_equal(m.add_noise()[0], np.array([0, 0, 0, 0, 0]))
    assert np.array_equal(m.add_noise()[1], np.array([0, 0, 0, 0, 0]))
    m = Matching(H, error_probabilities=1.0)
    assert np.array_equal(m.add_noise()[0], np.array([1, 1, 1, 1, 1]))
    assert np.array_equal(m.add_noise()[1], np.array([0, 0, 0, 0, 0]))


def test_weighted_mwpm_from_array():
    H = csc_matrix([[1, 0], [1, 1], [0, 1]])
    m = Matching(H, spacelike_weights=np.array([1., 2.]))
    with pytest.raises(ValueError):
        m = Matching(H, spacelike_weights=np.array([1.]))


def test_load_matching_from_dense_array():
    H = np.array([[1, 1, 0], [0, 1, 1]])
    m = Matching()
    m.load_from_check_matrix(H)


@pytest.mark.parametrize("t_weights,expected_edges",
                         [
                             (
                                     [0.5, 1.5],
                                     {((0, 1), 0.7), ((2, 3), 0.7), ((4, 5), 0.7),
                                      ((0, 6), 0.3), ((2, 6), 0.3), ((4, 6), 0.3),
                                      ((1, 6), 0.9), ((3, 6), 0.9), ((5, 6), 0.9),
                                      ((0, 2), 0.5), ((2, 4), 0.5), ((1, 3), 1.5), ((3, 5), 1.5)}
                             ),
                             (
                                     np.array([0.5, 1.5]),
                                     {((0, 1), 0.7), ((2, 3), 0.7), ((4, 5), 0.7),
                                      ((0, 6), 0.3), ((2, 6), 0.3), ((4, 6), 0.3),
                                      ((1, 6), 0.9), ((3, 6), 0.9), ((5, 6), 0.9),
                                      ((0, 2), 0.5), ((2, 4), 0.5), ((1, 3), 1.5), ((3, 5), 1.5)}
                             ),
                             (
                                     1.2,
                                     {((0, 1), 0.7), ((2, 3), 0.7), ((4, 5), 0.7),
                                      ((0, 6), 0.3), ((2, 6), 0.3), ((4, 6), 0.3),
                                      ((1, 6), 0.9), ((3, 6), 0.9), ((5, 6), 0.9),
                                      ((0, 2), 1.2), ((2, 4), 1.2), ((1, 3), 1.2), ((3, 5), 1.2)}
                             )
                         ]
                         )
def test_timelike_weights(t_weights, expected_edges):
    H = np.array([[1, 1, 0], [0, 1, 1]])
    m = Matching()
    m.load_from_check_matrix(H, spacelike_weights=np.array([0.3, 0.7, 0.9]),
                             timelike_weights=t_weights, repetitions=3)
    es = set((tuple(sorted([u, v])), d["weight"]) for u, v, d in m.edges())
    assert es == expected_edges


@pytest.mark.parametrize("t_weights", [[0.1, 0.01, 3], "A"])
def test_wrong_timelike_weights_raises_valueerror(t_weights):
    H = np.array([[1, 1, 0], [0, 1, 1]])
    with pytest.raises(ValueError):
        m = Matching()
        m.load_from_check_matrix(H, spacelike_weights=np.array([0.3, 0.7, 0.9]),
                                 timelike_weights=t_weights, repetitions=3)


@pytest.mark.parametrize("p_meas,expected_edges,repetitions",
                         [
                             (
                                     [0.15, 0.25],
                                     {((0, 1), 0.2), ((2, 3), 0.2), ((4, 5), 0.2),
                                      ((0, 6), 0.1), ((2, 6), 0.1), ((4, 6), 0.1),
                                      ((1, 6), 0.3), ((3, 6), 0.3), ((5, 6), 0.3),
                                      ((0, 2), 0.15), ((2, 4), 0.15), ((1, 3), 0.25), ((3, 5), 0.25)},
                                     3
                             ),
                             (
                                     np.array([0.15, 0.25]),
                                     {((0, 1), 0.2), ((2, 3), 0.2),
                                      ((0, 4), 0.1), ((2, 4), 0.1),
                                      ((1, 4), 0.3), ((3, 4), 0.3),
                                      ((0, 2), 0.15), ((1, 3), 0.25)},
                                     2
                             )
                         ]
                         )
def test_measurement_error_probabilities(p_meas, expected_edges, repetitions):
    m = Matching(
        [[1, 1, 0], [0, 1, 1]],
        error_probabilities=[0.1, 0.2, 0.3],
        measurement_error_probabilities=p_meas,
        repetitions=repetitions
    )
    es = set((tuple(sorted([u, v])), d["error_probability"]) for u, v, d in m.edges())
    assert es == expected_edges

    # Check measurement_error_probability also accepted
    m = Matching(
        [[1, 1, 0], [0, 1, 1]],
        error_probabilities=[0.1, 0.2, 0.3],
        measurement_error_probability=p_meas,
        repetitions=repetitions
    )
    es = set((tuple(sorted([u, v])), d["error_probability"]) for u, v, d in m.edges())
    assert es == expected_edges


@pytest.mark.parametrize("m_errors", [[0.1, 0.01, 3], "A"])
def test_wrong_measurement_error_probabilities_raises_valueerror(m_errors):
    H = np.array([[1, 1, 0], [0, 1, 1]])
    with pytest.raises(ValueError):
        m = Matching()
        m.load_from_check_matrix(H, spacelike_weights=np.array([0.3, 0.7, 0.9]),
                                 measurement_error_probabilities=m_errors, repetitions=3)


def test_measurement_error_probabilities_and_probability_raises_value_error():
    H = np.array([[1, 1, 0], [0, 1, 1]])
    with pytest.raises(ValueError):
        m = Matching()
        m.load_from_check_matrix(H, spacelike_weights=np.array([0.3, 0.7, 0.9]),
                                 measurement_error_probabilities=[0.1, 0.1], repetitions=3,
                                 measurement_error_probability=[0.1, 0.1])
