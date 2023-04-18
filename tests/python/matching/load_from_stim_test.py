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

from pathlib import Path

import pytest

from pymatching.matching import Matching


def test_load_from_stim_objects():
    stim = pytest.importorskip("stim")
    c = stim.Circuit.generated("surface_code:rotated_memory_x", distance=5, rounds=5,
                               after_clifford_depolarization=0.01,
                               before_measure_flip_probability=0.01,
                               after_reset_flip_probability=0.01,
                               before_round_data_depolarization=0.01)
    dem = c.detector_error_model(decompose_errors=True)
    m = Matching.from_detector_error_model(dem)
    assert m.num_detectors == dem.num_detectors
    assert m.num_fault_ids == dem.num_observables
    assert m.num_edges == 502
    m2 = Matching(dem)
    assert m2.num_detectors == dem.num_detectors
    assert m2.num_fault_ids == dem.num_observables
    assert m2.num_edges == 502
    m3 = Matching.from_stim_circuit(c)
    assert m3.num_detectors == dem.num_detectors
    assert m3.num_fault_ids == dem.num_observables
    assert m3.num_edges == 502


def test_load_from_stim_files(data_dir: Path):
    circuit_path = data_dir / "negative_weight_circuit.stim"
    m = Matching.from_stim_circuit_file(str(circuit_path))
    assert m.num_detectors == 2
    assert m.num_edges == 2
    assert m.num_fault_ids == 1
    dem_path = data_dir / "negative_weight_circuit.dem"
    m2 = Matching.from_detector_error_model_file(str(dem_path))
    assert m2.edges() == m.edges()
    with pytest.raises(ValueError):
        Matching.from_stim_circuit_file("fake_filename.stim")
    with pytest.raises(ValueError):
        Matching.from_detector_error_model_file("fake_filename.dem")
    with pytest.raises(ValueError):
        Matching.from_stim_circuit_file(str(dem_path))
    with pytest.raises(IndexError):
        Matching.from_detector_error_model_file(str(circuit_path))


def test_load_from_stim_wrong_type_raises_type_error():
    stim = pytest.importorskip("stim")
    c = stim.Circuit.generated("surface_code:rotated_memory_x", distance=3, rounds=1,
                               after_clifford_depolarization=0.01)
    with pytest.raises(TypeError):
        Matching.from_detector_error_model(c)
    with pytest.raises(TypeError):
        Matching.from_stim_circuit(c.detector_error_model(decompose_errors=True))


def test_load_from_dem_without_stim_raises_type_error():
    try:
        import stim  # noqa
    except ImportError:
        with pytest.raises(TypeError):
            Matching.from_detector_error_model("test")
        with pytest.raises(TypeError):
            Matching.from_stim_circuit("test")
