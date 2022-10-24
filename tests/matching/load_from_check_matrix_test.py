# Copyright 2022 PyMatching Contributors

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
from pymatching._cpp_pymatching import sparse_column_check_matrix_to_matching_graph


def test_boundary_from_check_matrix():
    H = csr_matrix(np.array([[1, 1, 0, 0, 0], [0, 1, 1, 0, 0],
                             [0, 0, 1, 1, 0], [0, 0, 0, 1, 1]]))
    m = Matching(H=H)  # Checks `H` still accepted as keyword argument despite being deprecated
    assert m.boundary == {4}
    assert np.array_equal(m.decode(np.array([1, 0, 0, 0])), np.array([1, 0, 0, 0, 0]))
    assert np.array_equal(m.decode(np.array([0, 1, 0, 0])), np.array([1, 1, 0, 0, 0]))
    assert np.array_equal(m.decode(np.array([0, 1, 1, 0])), np.array([0, 0, 1, 0, 0]))
    assert np.array_equal(m.decode(np.array([0, 0, 1, 0])), np.array([0, 0, 0, 1, 1]))


def test_nonzero_matrix_elements_not_one_raises_value_error():
    H = csr_matrix(np.array([[0, 2.01, 2.01], [1.01, 1.01, 0]]))
    with pytest.raises(ValueError):
        Matching(H)


def test_too_many_checks_per_column_raises_value_error():
    H = csr_matrix(np.array([[1, 1, 0, 0],
                             [1, 0, 1, 0],
                             [1, 0, 0, 1]]))
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
    Matching(H, spacelike_weights=np.array([1., 2.]))
    with pytest.raises(ValueError):
        Matching(H, spacelike_weights=np.array([1.]))


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
                             ),
                             (
                                     0.17,
                                     {((0, 1), 0.2), ((2, 3), 0.2),
                                      ((0, 4), 0.1), ((2, 4), 0.1),
                                      ((1, 4), 0.3), ((3, 4), 0.3),
                                      ((0, 2), 0.17), ((1, 3), 0.17)},
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
        m.load_from_check_matrix(H, weights=np.array([0.3, 0.7, 0.9]),
                                 measurement_error_probabilities=m_errors, repetitions=3)


def test_measurement_error_probabilities_and_probability_raises_value_error():
    H = np.array([[1, 1, 0], [0, 1, 1]])
    with pytest.raises(ValueError):
        m = Matching()
        m.load_from_check_matrix(H, spacelike_weights=np.array([0.3, 0.7, 0.9]),
                                 measurement_error_probabilities=[0.1, 0.1], repetitions=3,
                                 measurement_error_probability=[0.1, 0.1])


def test_cpp_csc_matrix_to_matching_graph():
    H = csc_matrix(np.array([[1, 1, 0, 0, 0, 0],
                             [0, 1, 1, 0, 0, 1],
                             [0, 0, 1, 1, 0, 1],
                             [0, 0, 0, 1, 1, 0]]))
    weights = np.array([1.0, 2.0, 3.0, 2.0, 1.5, 5.0])
    error_probabilities = np.array([0.4, 0.3, 0.4, 0.2, 0.1, 0.01])
    g = sparse_column_check_matrix_to_matching_graph(H, weights, error_probabilities, merge_strategy="smallest-weight",
                                                     use_virtual_boundary_node=False)
    assert g.get_edges() == [
        (0, 4, {"fault_ids": {0}, "weight": 1.0, "error_probability": 0.4}),
        (0, 1, {"fault_ids": {1}, "weight": 2.0, "error_probability": 0.3}),
        (1, 2, {"fault_ids": {2}, "weight": 3.0, "error_probability": 0.4}),
        (2, 3, {"fault_ids": {3}, "weight": 2.0, "error_probability": 0.2}),
        (3, 4, {"fault_ids": {4}, "weight": 1.5, "error_probability": 0.1}),
    ]
    assert g.get_boundary() == {4}
    assert g.get_num_nodes() == 5
    with pytest.raises(ValueError):
        # Check that scipy.sparse.csr_matrix is not accepted
        g = sparse_column_check_matrix_to_matching_graph(csr_matrix(H.T), weights, error_probabilities,
                                                         merge_strategy="smallest-weight",
                                                         use_virtual_boundary_node=False)
    g = sparse_column_check_matrix_to_matching_graph(H, weights, error_probabilities, merge_strategy="replace",
                                                     use_virtual_boundary_node=True)
    assert g.get_edges() == [
        (0, None, {"fault_ids": {0}, "weight": 1.0, "error_probability": 0.4}),
        (0, 1, {"fault_ids": {1}, "weight": 2.0, "error_probability": 0.3}),
        (1, 2, {"fault_ids": {5}, "weight": 5.0, "error_probability": 0.01}),
        (2, 3, {"fault_ids": {3}, "weight": 2.0, "error_probability": 0.2}),
        (3, None, {"fault_ids": {4}, "weight": 1.5, "error_probability": 0.1}),
    ]
    assert g.get_boundary() == set()
    assert g.get_num_nodes() == 4
    p_meas = np.array([0.1, 0.2, 0.15, 0.25])
    t_weights = np.log((1 - p_meas) / p_meas)
    g = sparse_column_check_matrix_to_matching_graph(H, weights, error_probabilities, merge_strategy="smallest-weight",
                                                     use_virtual_boundary_node=False, num_repetitions=3,
                                                     timelike_weights=t_weights,
                                                     measurement_error_probabilities=t_weights
                                                     )

    with pytest.raises(ValueError):
        p_meas = np.array([[0.1, 0.2], [0.15, 0.25]])
        t_weights = np.log((1 - p_meas) / p_meas)
        g = sparse_column_check_matrix_to_matching_graph(H, weights, error_probabilities,
                                                         merge_strategy="smallest-weight",
                                                         use_virtual_boundary_node=False, num_repetitions=3,
                                                         timelike_weights=t_weights,
                                                         measurement_error_probabilities=t_weights
                                                         )


def test_load_from_check_matrix_with_faults():
    H = csc_matrix(np.array([[1, 1, 0, 0, 0],
                             [0, 1, 1, 0, 0],
                             [0, 0, 1, 1, 0],
                             [0, 0, 0, 1, 1]]))
    F = csc_matrix(np.array([[0, 0, 1, 1, 1],
                             [1, 1, 1, 0, 0]]))
    m = Matching()
    m.load_from_check_matrix(H=H, faults_matrix=F)
    assert m.num_fault_ids == 2
    assert m.edges() == [(0, 4, {'error_probability': -1.0, 'fault_ids': {1}, 'weight': 1.0}),
                         (0, 1, {'error_probability': -1.0, 'fault_ids': {1}, 'weight': 1.0}),
                         (1, 2, {'error_probability': -1.0, 'fault_ids': {0, 1}, 'weight': 1.0}),
                         (2, 3, {'error_probability': -1.0, 'fault_ids': {0}, 'weight': 1.0}),
                         (3, 4, {'error_probability': -1.0, 'fault_ids': {0}, 'weight': 1.0})]
    m = Matching.from_check_matrix([[1, 1, 0], [0, 1, 1]], faults_matrix=[[1, 1, 1]], use_virtual_boundary_node=True)
    assert m.num_fault_ids == 1
    assert m.edges() == [(0, None, {'error_probability': -1.0, 'fault_ids': {0}, 'weight': 1.0}),
                         (0, 1, {'error_probability': -1.0, 'fault_ids': {0}, 'weight': 1.0}),
                         (1, None, {'error_probability': -1.0, 'fault_ids': {0}, 'weight': 1.0})]
    m = Matching()
    m.load_from_check_matrix([[1, 1, 0], [0, 1, 1]], faults_matrix=csc_matrix((0, 3), dtype=np.uint8),
                             use_virtual_boundary_node=True)
    assert m.num_fault_ids == 0
    assert m.edges() == [(0, None, {'error_probability': -1.0, 'fault_ids': set(), 'weight': 1.0}),
                         (0, 1, {'error_probability': -1.0, 'fault_ids': set(), 'weight': 1.0}),
                         (1, None, {'error_probability': -1.0, 'fault_ids': set(), 'weight': 1.0})]
    H = csc_matrix(np.array([[1, 1, 0],
                             [0, 1, 1],
                             [0, 0, 1]]))
    F = csc_matrix(np.array([[0, 0, 1],
                             [1, 1, 1],
                             [1, 0, 1],
                             [0, 0, 1],
                             [1, 1, 0]]))
    m = Matching.from_check_matrix(check_matrix=H, faults_matrix=F, use_virtual_boundary_node=True)
    assert m.num_fault_ids == 5
    assert m.edges() == [(0, None, {'error_probability': -1.0, 'fault_ids': {1, 2, 4}, 'weight': 1.0}),
                         (0, 1, {'error_probability': -1.0, 'fault_ids': {1, 4}, 'weight': 1.0}),
                         (1, 2, {'error_probability': -1.0, 'fault_ids': {0, 1, 2, 3}, 'weight': 1.0})]


def test_no_check_matrix_raises_value_error():
    with pytest.raises(ValueError):
        m = Matching()
        m.load_from_check_matrix()


def test_from_empty_check_matrix():
    m = Matching.from_check_matrix([[0, 0, 0], [0, 0, 0]], faults_matrix=[[0, 0, 0], [1, 0, 1], [1, 1, 0]])
    assert m.num_fault_ids == 3
    assert m.num_edges == 0
    assert m.num_detectors == 2


def test_type_error_from_check_matrix():
    with pytest.raises(TypeError):
        Matching.from_check_matrix([[0, 1, 1]], faults_matrix=["A", 1])


def test_spacelike_weights_and_weights_raises_value_error():
    with pytest.raises(ValueError):
        Matching.from_check_matrix([[0, 1, 1]], weights=[1, 1, 1], spacelike_weights=[2, 2, 2])
