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
import stim

import pymatching

# Edges (0, 1) and (2, 3) each have marginal probability 0.2, and P((2, 3) | (0, 1)) = 0.08 / 0.2 = 0.4.
DEM = stim.DetectorErrorModel("""
    error(0.14285714285714285) D0 D1
    error(0.14285714285714285) D2 D3
    error(0.08) D0 D1 ^ D2 D3
""")
SYNDROME = np.array([1, 1, 1, 1], dtype=np.uint8)


def test_default_reweighting_strength_is_one():
    m = pymatching.Matching.from_detector_error_model(DEM, enable_correlations=True)
    shots = np.array([[0, 0, 0, 0], [1, 1, 0, 0], [0, 0, 1, 1], [1, 1, 1, 1]], dtype=np.uint8)
    default = m.decode_batch(shots, return_weights=True, enable_correlations=True)
    explicit = m.decode_batch(shots, return_weights=True, enable_correlations=True, reweighting_strength=1.0)
    np.testing.assert_array_equal(default[0], explicit[0])
    np.testing.assert_array_equal(default[1], explicit[1])
    np.testing.assert_array_equal(
        m.decode_to_edges_array(SYNDROME, enable_correlations=True),
        m.decode_to_edges_array(SYNDROME, enable_correlations=True, reweighting_strength=1.0),
    )


@pytest.mark.parametrize("strength", [0.6, 1.0, 1.5])
def test_reweighting_strength_scales_conditional_probability(strength):
    m = pymatching.Matching.from_detector_error_model(DEM, enable_correlations=True)
    _, weight = m.decode(SYNDROME, return_weight=True, enable_correlations=True, reweighting_strength=strength)
    p = min(0.5, strength * 0.4)
    assert weight == pytest.approx(2 * np.log((1 - p) / p), abs=1e-3)


def test_changing_reweighting_strength_between_calls():
    m = pymatching.Matching.from_detector_error_model(DEM, enable_correlations=True)
    weights = [
        m.decode(SYNDROME, return_weight=True, enable_correlations=True, reweighting_strength=s)[1]
        for s in (1.0, 0.6, 1.0)
    ]
    assert weights[0] == weights[2] != weights[1]


def test_reweighting_strength_must_be_positive():
    m = pymatching.Matching.from_detector_error_model(DEM, enable_correlations=True)
    with pytest.raises(ValueError):
        m.decode(SYNDROME, enable_correlations=True, reweighting_strength=0)
