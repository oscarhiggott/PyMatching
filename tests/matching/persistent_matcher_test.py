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
import pytest
from scipy.sparse import csc_matrix, hstack, kron, eye

from pymatching import Matching, PersistentMatcher


def repetition_code(n):
    return csc_matrix(np.eye(n - 1, n, dtype=np.uint8) + np.eye(n - 1, n, 1, dtype=np.uint8))


def toric_code_x_checks(L):
    h = csc_matrix(np.eye(L, dtype=np.uint8) + np.roll(np.eye(L, dtype=np.uint8), 1, axis=1))
    return csc_matrix(hstack([kron(h, eye(L)), kron(eye(L), h.T)]), dtype=np.uint8)


def surface_code_x_checks(d):
    h = repetition_code(d)
    m, n = h.shape
    return csc_matrix(hstack([kron(h, eye(n)), kron(eye(m), h.T)]), dtype=np.uint8)


def with_parallel_edges(h, rng):
    """Duplicate some columns, so that several faults share an edge (merged by smallest weight)."""
    cols = rng.choice(h.shape[1], size=h.shape[1] // 3, replace=False)
    return csc_matrix(hstack([h, h[:, cols]]), dtype=np.uint8)


def random_faults_matrix(rng, num_observables, num_cols):
    return csc_matrix((rng.random((num_observables, num_cols)) < 0.2).astype(np.uint8))


def assert_matches_fresh_matching(h, faults_matrix, rng, num_weight_sets=8, shots=15, negative=False):
    pm = PersistentMatcher(h, faults_matrix=faults_matrix)
    for _ in range(num_weight_sets):
        w = rng.uniform(0.1, 10.0, h.shape[1])
        if negative:
            w[rng.random(h.shape[1]) < 0.1] *= -1
        pm.set_weights(w)
        m = Matching.from_check_matrix(h, weights=w, faults_matrix=faults_matrix)
        for _ in range(shots):
            errors = (rng.random(h.shape[1]) < 0.1).astype(np.uint8)
            syndrome = (h @ errors) % 2
            assert np.array_equal(pm.decode(syndrome), m.decode(syndrome))


@pytest.mark.parametrize("negative", [False, True])
@pytest.mark.parametrize(
    "make_h",
    [
        lambda rng: repetition_code(9),
        lambda rng: surface_code_x_checks(5),
        lambda rng: toric_code_x_checks(4),
        lambda rng: with_parallel_edges(surface_code_x_checks(4), rng),
    ],
)
def test_persistent_matcher_returns_the_correction_of_a_fresh_matching(make_h, negative):
    rng = np.random.default_rng(0)
    assert_matches_fresh_matching(make_h(rng), None, rng, negative=negative)


@pytest.mark.parametrize("num_observables", [1, 3, 64, 70])
def test_persistent_matcher_with_faults_matrix(num_observables):
    rng = np.random.default_rng(1)
    h = with_parallel_edges(surface_code_x_checks(5), rng)
    f = random_faults_matrix(rng, num_observables, h.shape[1])
    assert_matches_fresh_matching(h, f, rng, negative=True)


def test_persistent_matcher_integral_weights():
    # integral weights are not rescaled before discretisation, as in Matching
    rng = np.random.default_rng(2)
    h = surface_code_x_checks(5)
    pm = PersistentMatcher(h)
    for _ in range(5):
        w = rng.integers(1, 20, h.shape[1]).astype(float)
        pm.set_weights(w)
        m = Matching.from_check_matrix(h, weights=w)
        for _ in range(10):
            syndrome = (h @ (rng.random(h.shape[1]) < 0.1).astype(np.uint8)) % 2
            assert np.array_equal(pm.decode(syndrome), m.decode(syndrome))


def test_persistent_matcher_scalar_and_default_weights():
    h = repetition_code(7)
    syndrome = np.array([0, 1, 0, 0, 1, 0])
    expected = Matching.from_check_matrix(h).decode(syndrome)
    assert np.array_equal(PersistentMatcher(h).decode(syndrome), expected)
    pm = PersistentMatcher(h, weights=2.5)
    assert np.array_equal(pm.decode(syndrome), expected)
    assert pm.num_detectors == 6
    assert pm.num_fault_ids == 7


def test_persistent_matcher_input_validation():
    h = repetition_code(5)
    pm = PersistentMatcher(h)
    with pytest.raises(ValueError):
        pm.set_weights(np.ones(4))
    with pytest.raises(ValueError):
        pm.decode(np.zeros(5))
    with pytest.raises(ValueError):
        PersistentMatcher(h, faults_matrix=np.ones((1, 4)))
    with pytest.raises(ValueError):
        PersistentMatcher(2 * h)
    with pytest.raises(ValueError):
        PersistentMatcher(np.ones((3, 2), dtype=np.uint8))  # a column with three non-zero elements


@pytest.mark.parametrize("bad", [np.nan, np.inf, -np.inf])
def test_persistent_matcher_rejects_non_finite_weights(bad):
    h = repetition_code(5)
    w = np.ones(h.shape[1])
    w[2] = bad
    with pytest.raises(ValueError, match="finite"):
        PersistentMatcher(h, weights=w)
    with pytest.raises(ValueError, match="finite"):
        PersistentMatcher(h, weights=bad)
    pm = PersistentMatcher(h)
    with pytest.raises(ValueError, match="finite"):
        pm.set_weights(w)
    with pytest.raises(ValueError, match="finite"):
        pm.set_weights(bad)


@pytest.mark.parametrize("num_observables", [None, 1, 3])
def test_persistent_matcher_num_fault_ids_counts_observables(num_observables):
    h = repetition_code(5)
    f = None if num_observables is None else np.ones((num_observables, h.shape[1]), dtype=np.uint8)
    pm = PersistentMatcher(h, faults_matrix=f)
    m = Matching.from_check_matrix(h, faults_matrix=f)
    assert pm.num_fault_ids == m.num_fault_ids
    assert pm.decode(np.zeros(h.shape[0], dtype=np.uint8)).shape == (pm.num_fault_ids,)
