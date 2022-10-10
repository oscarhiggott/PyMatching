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

from typing import Union, List, TYPE_CHECKING
if TYPE_CHECKING:
    import pymatching

import numpy as np
import scipy
from scipy.sparse import csc_matrix

from pymatching._cpp_pymatching import MatchingGraph, sparse_column_check_matrix_to_matching_graph


def load_from_check_matrix(self: 'pymatching.Matching',
                           H: Union[scipy.sparse.spmatrix, np.ndarray, List[List[int]]],
                           weights: Union[float, np.ndarray, List[float]] = None,
                           error_probabilities: Union[float, np.ndarray, List[float]] = None,
                           repetitions: int = None,
                           timelike_weights: Union[float, np.ndarray, List[float]] = None,
                           measurement_error_probabilities: Union[float, np.ndarray, List[float]] = None,
                           *,
                           merge_strategy: str = "smallest-weight",
                           use_virtual_boundary_node: bool = False,
                           **kwargs
                           ) -> None:
    """
    Load a matching graph from a check matrix
    Parameters
    ----------
    H : `scipy.spmatrix` or `numpy.ndarray` or List[List[int]]
        The quantum code to be decoded with minimum-weight perfect
        matching, given as a binary check matrix (scipy sparse
        matrix or numpy.ndarray)
    weights : float or numpy.ndarray, optional
        If `H` is given as a scipy or numpy array, `weights` gives the weights
        of edges in the matching graph corresponding to columns of `H`.
        If `weights` is a numpy.ndarray, it should be a 1D array with length
        equal to `H.shape[1]`. If weights is a float, it is used as the weight for all
        edges corresponding to columns of `H`. By default None, in which case
        all weights are set to 1.0
        This argument was renamed from `spacelike_weights` in PyMatching v2.0, but
        `spacelike_weights` is still accepted in place of `weights` for backward compatibility.
    error_probabilities : float or numpy.ndarray, optional
        The probabilities with which an error occurs on each edge associated with a
        column of H. If a
        single float is given, the same error probability is used for each
        column. If a numpy.ndarray of floats is given, it must have a
        length equal to the number of columns in H. This parameter is only
        needed for the Matching.add_noise method, and not for decoding.
        By default None
    repetitions : int, optional
        The number of times the stabiliser measurements are repeated, if
        the measurements are noisy. By default None
    timelike_weights : float or numpy.ndarray, optional
        If `repetitions>1`, `timelike_weights` gives the weight of
        timelike edges. If a float is given, all timelike edges weights are set to
        the same value. If a numpy array of size `(H.shape[0],)` is given, the
        edge weight for each vertical timelike edge associated with the `i`th check (row)
        of `H` is set to `timelike_weights[i]`. By default None, in which case all
        timelike weights are set to 1.0
    measurement_error_probabilities : float or numpy.ndarray, optional
        If `repetitions>1`, gives the probability of a measurement
        error to be used for the add_noise method. If a float is given, all measurement
        errors are set to the same value. If a numpy array of size `(H.shape[0],)` is given,
        the error probability for each vertical timelike edge associated with the `i`th check
        (row) of `H` is set to `measurement_error_probabilities[i]`. This argument can also be
        given using the keyword argument `measurement_error_probability` to maintain backward
        compatibility with previous versions of Pymatching. By default None
    merge_strategy: str, optional
        Which strategy to use when adding an edge (`node1`, `node2`) that is already in the graph. The available
        options are "disallow", "independent", "smallest-weight", "keep-original" and "replace". "disallow" raises a
        `ValueError` if the edge (`node1`, `node2`) is already present. The "independent" strategy assumes that
        the existing edge (`node1`, `node2`) and the edge being added represent independent error mechanisms, and
        they are merged into a new edge with updated weights and error_probabilities accordingly (it is assumed
        that each weight represents the log-likelihood ratio log((1-p)/p) where p is the `error_probability` and
        where the natural logarithm is used. The fault_ids associated with the existing edge are kept only, since
        the code has distance 2 if parallel edges have different fault_ids anyway). The "smallest-weight" strategy
        keeps only the new edge if it has a smaller weight than the existing edge, otherwise the graph is left
        unchanged. The "keep-original" strategy keeps only the existing edge, and ignores the edge being added.
        The "replace" strategy always keeps the edge being added, replacing the existing edge.
        By default, "smallest-weight"
    use_virtual_boundary_node: bool, optional
        This option determines how columns are handled if they contain only a single 1 (representing a boundary edge).
        Consider a column contains a single 1 at row index i. If `use_virtual_boundary_node=False`, then this column
        will be handled by adding an edge `(i, H.shape[0])`, and marking the node `H.shape[0]` as a boundary node with
        `Matching.set_boundary(H.shape[0])`. The resulting graph will contain `H.shape[0]+1` nodes, the largest of
        which is the boundary node. If `use_virtual_boundary_node=True` then instead the boundary is a virtual node, and
        this column is handled with `Matching.add_boundary_edge(i, ...)`. The resulting graph will contain `H.shape[0]`
        nodes, and there is no boundary node. Both options are handled identically by the decoder, although
        `use_virtual_boundary_node=True` is recommended since it is simpler (with a one-to-one correspondence between
         nodes and rows of H), and is also slightly more efficient. By default, False (for backward compatibility)
    Examples
    --------
    >>> import pymatching
    >>> m = pymatching.Matching([[1, 1, 0, 0], [0, 1, 1, 0], [0, 0, 1, 1]])
    >>> m
    <pymatching.Matching object with 3 detectors, 1 boundary node, and 4 edges>

    Matching objects can also be initialised from a sparse scipy matrix:
    >>> import pymatching
    >>> from scipy.sparse import csc_matrix
    >>> H = csc_matrix([[1, 1, 0], [0, 1, 1]])
    >>> m = pymatching.Matching(H)
    >>> m
    <pymatching.Matching object with 2 detectors, 1 boundary node, and 3 edges>
    """
    if not isinstance(H, csc_matrix):
        try:
            H = csc_matrix(H)
        except TypeError:
            raise TypeError("H must be convertible to a `scipy.csc_matrix`")
    num_edges = H.shape[1]

    slw = kwargs.get("spacelike_weights")
    if weights is None and slw is not None:
        weights = slw
    elif weights is not None and slw is not None:
        raise ValueError("Both `weights` and `spacelike_weights` were provided as arguments, but these "
                         "two arguments are equivalent. Please provide only `weights` as an argument, as "
                         "the `spacelike_weights` argument has been deprecated.")

    weights = 1.0 if weights is None else weights
    if isinstance(weights, (int, float, np.integer, np.floating)):
        weights = np.ones(num_edges, dtype=float) * weights
    weights = np.asarray(weights)

    if error_probabilities is None:
        error_probabilities = np.ones(num_edges) * -1
    elif isinstance(error_probabilities, (int, float)):
        error_probabilities = np.ones(num_edges) * error_probabilities

    H.eliminate_zeros()

    repetitions = 1 if repetitions is None else repetitions

    if repetitions > 1:
        timelike_weights = 1.0 if timelike_weights is None else timelike_weights
        if isinstance(timelike_weights, (int, float, np.integer, np.floating)):
            timelike_weights = np.ones(H.shape[0], dtype=float) * timelike_weights
        elif isinstance(timelike_weights, (np.ndarray, list)):
            timelike_weights = np.array(timelike_weights, dtype=float)
        else:
            raise ValueError("timelike_weights should be a float or a 1d numpy array")

        mep = kwargs.get("measurement_error_probability")
        if measurement_error_probabilities is not None and mep is not None:
            raise ValueError("Both `measurement_error_probabilities` and `measurement_error_probability` "
                             "were provided as arguments. Please "
                             "provide `measurement_error_probabilities` instead of `measurement_error_probability` "
                             "as an argument, as use of `measurement_error_probability` has been deprecated.")
        if measurement_error_probabilities is None and mep is not None:
            measurement_error_probabilities = mep

        p_meas = measurement_error_probabilities if measurement_error_probabilities is not None else -1
        if isinstance(p_meas, (int, float, np.integer, np.floating)):
            p_meas = np.ones(H.shape[0], dtype=float)
        elif isinstance(p_meas, (np.ndarray, list)):
            p_meas = np.array(p_meas, dtype=float)
        else:
            raise ValueError("measurement_error_probabilities should be a float or 1d numpy array")
    else:
        timelike_weights = None
        p_meas = None

    self._matching_graph = sparse_column_check_matrix_to_matching_graph(H, weights, error_probabilities, merge_strategy,
                                                                        use_virtual_boundary_node, repetitions,
                                                                        timelike_weights, p_meas)
