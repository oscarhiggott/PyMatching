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

from typing import List, Optional, Union

import numpy as np
from scipy.sparse import csc_matrix, identity, spmatrix

from pymatching._cpp_pymatching import PersistentMatcher as _CppPersistentMatcher


def _to_binary_csc(matrix: Union[csc_matrix, spmatrix, np.ndarray, List[List[int]]], name: str) -> csc_matrix:
    m = csc_matrix(matrix)
    m.eliminate_zeros()
    if m.nnz and not np.all(m.data == 1):
        raise ValueError(f"{name} must be a binary matrix")
    m = m.astype(np.uint8)
    m.sort_indices()
    return m


class PersistentMatcher:
    """A matching graph built once from a check matrix, whose edge weights are updated in place.

    Decoding many syndromes of the same check matrix with different weights normally means building a new
    :py:class:`pymatching.Matching` for every weight vector, and building the graph can cost far more than the
    decode itself. A `PersistentMatcher` builds the graph once; :py:meth:`set_weights` then rewrites the edge
    weights in place, and :py:meth:`decode` gives the same result as::

        pymatching.Matching.from_check_matrix(check_matrix, weights=weights,
                                              faults_matrix=faults_matrix).decode(syndrome)

    that is, with the default ``merge_strategy="smallest-weight"``, no virtual boundary node and a single
    repetition. Weights may be negative, as for `Matching`.

    Parameters
    ----------
    check_matrix : scipy.sparse.csc_matrix or numpy.ndarray
        A binary check matrix with at most two non-zero elements per column.
    weights : float or numpy.ndarray, optional
        The initial weight of each column (fault). A float sets every column to that weight. By default 1.
    faults_matrix : scipy.sparse.csc_matrix or numpy.ndarray, optional
        A binary matrix with one column per column of `check_matrix`: ``faults_matrix[i, j] == 1`` if fault ``j``
        flips logical observable ``i``. :py:meth:`decode` then returns the predicted observables. By default the
        identity, so :py:meth:`decode` returns the correction, as `Matching.decode` does.
    """

    def __init__(
        self,
        check_matrix: Union[csc_matrix, spmatrix, np.ndarray, List[List[int]]],
        weights: Optional[Union[float, np.ndarray, List[float]]] = None,
        *,
        faults_matrix: Optional[Union[csc_matrix, spmatrix, np.ndarray, List[List[int]]]] = None,
    ) -> None:
        h = _to_binary_csc(check_matrix, "check_matrix")
        num_rows, num_cols = h.shape
        if faults_matrix is None:
            f = identity(num_cols, dtype=np.uint8, format="csc")
        else:
            f = _to_binary_csc(faults_matrix, "faults_matrix")
            if f.shape[1] != num_cols:
                raise ValueError(
                    f"faults_matrix has {f.shape[1]} columns but check_matrix has {num_cols}; they must be equal"
                )
        self._num_rows = num_rows
        self._num_cols = num_cols
        self._num_observables = f.shape[0]
        self._matcher = _CppPersistentMatcher(
            h.indptr.astype(np.int64),
            h.indices.astype(np.int64),
            num_rows,
            num_cols,
            f.indptr.astype(np.int64),
            f.indices.astype(np.int64),
            f.shape[0],
            self._weight_array(1.0 if weights is None else weights),
        )

    def _weight_array(self, weights: Union[float, np.ndarray, List[float]]) -> np.ndarray:
        if np.isscalar(weights):
            w = np.full(self._num_cols, float(weights))
        else:
            w = np.ascontiguousarray(weights, dtype=np.float64)
            if w.shape != (self._num_cols,):
                raise ValueError(f"weights must have shape ({self._num_cols},), not {w.shape}")
        if not np.all(np.isfinite(w)):
            raise ValueError("weights must be finite")
        return w

    def set_weights(self, weights: Union[float, np.ndarray, List[float]]) -> None:
        """Set the weight of every column (fault).

        Parameters
        ----------
        weights : float or numpy.ndarray
            One weight per column of the check matrix, or a float for all of them.
        """
        self._matcher.set_weights(self._weight_array(weights))

    def decode(self, syndrome: Union[np.ndarray, List[int]]) -> np.ndarray:
        """Decode a syndrome with the current weights.

        Parameters
        ----------
        syndrome : numpy.ndarray
            A binary vector with one element per row of the check matrix.

        Returns
        -------
        numpy.ndarray
            The predicted observables (one per row of `faults_matrix`); the correction if no `faults_matrix` was
            given.
        """
        z = np.asarray(syndrome)
        if z.shape != (self._num_rows,):
            raise ValueError(f"syndrome must have shape ({self._num_rows},), not {z.shape}")
        return self._matcher.decode(np.flatnonzero(z).astype(np.uint64))

    @property
    def num_detectors(self) -> int:
        """The number of detectors (rows of the check matrix)."""
        return self._num_rows

    @property
    def num_fault_ids(self) -> int:
        """The number of fault IDs (observables): the rows of `faults_matrix`, or the columns of the check matrix
        if no `faults_matrix` was given."""
        return self._num_observables
