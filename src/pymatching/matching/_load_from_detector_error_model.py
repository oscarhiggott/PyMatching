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

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    import pymatching

import stim

from pymatching._cpp_pymatching import detector_error_model_to_matching_graph


def load_from_detector_error_model(self: 'pymatching.Matching', model: stim.DetectorErrorModel) -> None:
    """
    Load from a `stim.DetectorErrorModel`.

    A `stim.DetectorErrorModel` (DEM) describes a circuit-level noise model in a quantum error correction protocol,
    and is defined in the
    Stim documentation: https://github.com/quantumlib/Stim/blob/main/doc/file_format_dem_detector_error_model.md.
    When loading from a DEM, there is a one-to-one correspondence with a detector in the DEM and a
    node in the `pymatching.Matching` graph, and each graphlike error in the DEM becomes an edge (or merged into
    a parallel edge) in the `pymatching.Matching` graph.
    A error instruction in the DEM is graphlike if it causes either one or two detection events, and can be
    either its own DEM instruction, or within a suggested decomposition of a larger DEM instruction.
    Error instruction in the DEM that cause more than two detection events and do not have a suggested
    decomposition into edges are ignored.
    There set of `fault_ids` assigned to a `pymatching.Matching` graph edge is the set of
    `logical_observable` indices associated with the corresponding graphlike fault mechanism in the DEM.
    Parallel edges are merged, with weights chosen on the assumption that the error mechanisms associated with the
    parallel edges are independent.
    If parallel edges have different `logical_observable` indices, this implies the code has distance 2, and only
     the `logical_observable` indices associated with the first added parallel edge are kept for the merged edge.
    If you are loading a `pymatching.Matching` graph from a DEM, you may be interested in
    using the sinter Python package for monte carlo sampling: https://pypi.org/project/sinter/.
    Parameters
    ----------
    model

    Returns
    -------

    """
    self._matching_graph = detector_error_model_to_matching_graph(str(model))
