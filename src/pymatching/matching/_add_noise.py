# Copyright 2020 Oscar Higgott

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#      http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from typing import Union, Tuple, TYPE_CHECKING
if TYPE_CHECKING:
    import pymatching

import numpy as np


def add_noise(self: 'pymatching.Matching') -> Union[Tuple[np.ndarray, np.ndarray], None]:
    """Add noise by flipping edges in the matching graph with
    a probability given by the error_probility edge attribute.
    The ``error_probability`` must be set for all edges for this
    method to run, otherwise it returns `None`.
    All boundary nodes are always given a 0 syndrome.
    Returns
    -------
    numpy.ndarray of dtype int
        Noise vector (binary numpy int array of length self.num_fault_ids)
    numpy.ndarray of dtype int
        Syndrome vector (binary numpy int array of length
        self.num_detectors if there is no boundary, or self.num_detectors+len(self.boundary)
        if there are boundary nodes)
    """
    if not self._matching_graph.all_edges_have_error_probabilities():
        return None
    return self._matching_graph.add_noise()
