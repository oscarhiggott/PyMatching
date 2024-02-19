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

from typing import Union, List, TYPE_CHECKING, Tuple, Set, Dict, Optional
import warnings

import numpy as np
import pymatching
import networkx as nx
from scipy.sparse import csc_matrix, spmatrix
import matplotlib.cbook

if TYPE_CHECKING:
    import stim  # pragma: no cover
    import rustworkx as rx  # pragma: no cover

import pymatching._cpp_pymatching as _cpp_pm


class Matching:
    """
    A class for constructing matching graphs and decoding using the minimum-weight perfect matching decoder.
    The matching graph can be constructed using the `Matching.add_edge` and `Matching.add_boundary_edge`
    methods. Alternatively, it can be loaded from a parity check matrix (a `scipy.sparse` matrix or `numpy.ndarray`
    with one or two non-zero elements in each column), a NetworkX or rustworkx graph, or from
    a `stim.DetectorErrorModel`.
    """

    def __init__(self,
                 graph: Union[csc_matrix, np.ndarray, "rx.PyGraph", nx.Graph, List[
                     List[int]], 'stim.DetectorErrorModel', spmatrix] = None,
                 weights: Union[float, np.ndarray, List[float]] = None,
                 error_probabilities: Union[float, np.ndarray, List[float]] = None,
                 repetitions: int = None,
                 timelike_weights: Union[float, np.ndarray, List[float]] = None,
                 measurement_error_probabilities: Union[float, np.ndarray, List[float]] = None,
                 **kwargs
                 ):
        r"""Constructor for the Matching class

        Parameters
        ----------
        graph : `scipy.spmatrix` or `numpy.ndarray` or `networkx.Graph` or `stim.DetectorErrorModel`, optional
            The matching graph to be decoded with minimum-weight perfect
            matching, given either as a binary parity check matrix (scipy sparse
            matrix or numpy.ndarray), a NetworkX or rustworkx graph, or a Stim DetectorErrorModel.
            Each edge in the NetworkX or rustworkx graph can have optional
            attributes ``fault_ids``, ``weight`` and ``error_probability``.
            ``fault_ids`` should be an int or a set of ints.
            Each fault id corresponds to a self-inverse fault that is flipped when the
            corresponding edge is flipped. These self-inverse faults could correspond to
            physical Pauli errors (physical frame changes)
            or to the logical observables that are flipped by the fault
            (a logical frame change, equivalent to an obersvable ID in an error instruction in a Stim
            detector error model). The `fault_ids` attribute was previously named `qubit_id` in an
            earlier version of PyMatching, and `qubit_id` is still accepted instead of `fault_ids` in order
            to maintain backward compatibility.
            Each ``weight`` attribute should be a non-negative float. If
            every edge is assigned an error_probability between zero and one,
            then the ``add_noise`` method can be used to simulate noise and
            flip edges independently in the graph. By default, None
        weights : float or numpy.ndarray, optional
            If `graph` is given as a scipy or numpy array, `weights` gives the weights
            of edges in the matching graph corresponding to columns of `graph`.
            If weights is a numpy.ndarray, it should be a 1D array with length
            equal to `graph.shape[1]`. If weights is a float, it is used as the weight for all
            edges corresponding to columns of `graph`. By default None, in which case
            all weights are set to 1.0
            This argument was renamed from `spacelike_weights` in PyMatching v2.0, but
            `spacelike_weights` is still accepted in place of `weights` for backward compatibility.
        error_probabilities : float or numpy.ndarray, optional
            The probabilities with which an error occurs on each edge corresponding
            to a column of the check matrix. If a
            single float is given, the same error probability is used for each
            edge. If a numpy.ndarray of floats is given, it must have a
            length equal to the number of columns in the check matrix. This parameter is only
            needed for the Matching.add_noise method, and not for decoding.
            By default None
        repetitions : int, optional
            The number of times the stabiliser measurements are repeated, if
            the measurements are noisy. This option is only used if `check_matrix` is
            provided as a check matrix, not a NetworkX graph. By default None
        timelike_weights : float, optional
            If `check_matrix` is given as a scipy or numpy array and `repetitions>1`,
            `timelike_weights` gives the weight of timelike edges.
            If a float is given, all timelike edges weights are set to
            the same value. If a numpy array of size `(check_matrix.shape[0],)` is given, the
            edge weight for each vertical timelike edge associated with the `i`th check (row)
            of `check_matrix` is set to `timelike_weights[i]`. By default None, in which case all
            timelike weights are set to 1.0
        measurement_error_probabilities : float, optional
            If `check_matrix` is given as a scipy or numpy array and `repetitions>1`,
            gives the probability of a measurement error to be used for
            the add_noise method. If a float is given, all measurement
            errors are set to the same value. If a numpy array of size `(check_matrix.shape[0],)` is given,
            the error probability for each vertical timelike edge associated with the `i`th check
            (row) of `check_matrix` is set to `measurement_error_probabilities[i]`. By default None
        **kwargs
            The remaining keyword arguments are passed to `Matching.load_from_check_matrix` if `graph` is a
            check matrix.

        Examples
        --------
        >>> import pymatching
        >>> import math
        >>> m = pymatching.Matching()
        >>> m.add_edge(0, 1, fault_ids={0}, weight=0.1)
        >>> m.add_edge(1, 2, fault_ids={1}, weight=0.15)
        >>> m.add_edge(2, 3, fault_ids={2, 3}, weight=0.2)
        >>> m.add_edge(0, 3, fault_ids={4}, weight=0.1)
        >>> m.set_boundary_nodes({3})
        >>> m
        <pymatching.Matching object with 3 detectors, 1 boundary node, and 4 edges>

        Matching objects can also be created from a check matrix (provided as a scipy.sparse matrix,
        dense numpy array, or list of lists):
        >>> import pymatching
        >>> m = pymatching.Matching([[1, 1, 0, 0], [0, 1, 1, 0], [0, 0, 1, 1]])
        >>> m
        <pymatching.Matching object with 3 detectors, 1 boundary node, and 4 edges>
        """
        self._matching_graph = _cpp_pm.MatchingGraph()
        if graph is None:
            graph = kwargs.get("H")
            if graph is None:
                return
            del kwargs["H"]
        # Networkx graph
        if isinstance(graph, nx.Graph):
            self.load_from_networkx(graph)
            return
        # Rustworkx PyGraph
        try:
            import rustworkx as rx
            if isinstance(graph, rx.PyGraph):
                self.load_from_rustworkx(graph)
                return
        except ImportError:  # pragma no cover
            pass
        # stim.DetectorErrorModel
        try:
            import stim
            if isinstance(graph, stim.DetectorErrorModel):
                self._load_from_detector_error_model(graph)
                return
        except ImportError:  # pragma no cover
            pass
        # scipy.csc_matrix
        try:
            graph = csc_matrix(graph)
        except TypeError:
            raise TypeError("The type of the input graph is not recognised. `graph` must be "
                            "a scipy.sparse or numpy matrix, networkx or rustworkx graph, or "
                            "stim.DetectorErrorModel.")
        self.load_from_check_matrix(graph, weights, error_probabilities,
                                    repetitions, timelike_weights, measurement_error_probabilities,
                                    **kwargs)

    def add_noise(self) -> Union[Tuple[np.ndarray, np.ndarray], None]:
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

    def _syndrome_array_to_detection_events(self, z: Union[np.ndarray, List[int]]) -> np.ndarray:
        try:
            z = np.array(z, dtype=np.uint8)
        except ValueError:
            raise ValueError("Syndrome must be of type numpy.ndarray or "
                             "convertible to numpy.ndarray, not {}".format(z))
        if len(z.shape) == 1 and (self.num_detectors <= z.shape[0]
                                  <= self.num_detectors + len(self.boundary)):
            detection_events = z.nonzero()[0]
        elif len(z.shape) == 2 and z.shape[0] * z.shape[1] == self.num_detectors:
            times, checks = z.T.nonzero()
            detection_events = times * z.shape[0] + checks
        else:
            raise ValueError("The shape ({}) of the syndrome vector z is not valid.".format(z.shape))
        return detection_events

    def decode(self,
               z: Union[np.ndarray, List[bool], List[int]],
               _legacy_num_neighbours: int = None,
               _legacy_return_weight: bool = None,
               *,
               return_weight: bool = False,
               **kwargs
               ) -> Union[np.ndarray, Tuple[np.ndarray, int]]:
        r"""
        Decode the syndrome `z` using minimum-weight perfect matching

        Parameters
        ----------
        z : numpy.ndarray
            A binary syndrome vector to decode. The number of elements in
            `z` should equal the number of nodes in the matching graph. If
            `z` is a 1D array, then `z[i]` is the syndrome at node `i` of
            the matching graph. If `z` is 2D then `z[i,j]` is the difference
            (modulo 2) between the (noisy) measurement of stabiliser `i` in time
            step `j+1` and time step `j` (for the case where the matching graph is
            constructed from a check matrix with `repetitions>1`).
        return_weight : bool, optional
            If `return_weight==True`, the sum of the weights of the edges in the
            minimum weight perfect matching is also returned. By default False

        Returns
        -------
        correction : numpy.ndarray or list[int]
            A 1D numpy array of ints giving the minimum-weight correction operator as a
            binary vector. The number of elements in `correction` is one greater than
            the largest fault ID. The ith element of `correction` is 1 if the
            minimum-weight perfect matching (MWPM) found by PyMatching contains an odd
            number of edges that have `i` as one of the `fault_ids`, and is 0 otherwise.
            If each edge in the matching graph is assigned a unique integer in its
            `fault_ids` attribute, then the locations of nonzero entries in `correction`
            correspond to the edges in the MWPM. However, `fault_ids` can instead be used,
            for example, to store IDs of the physical or logical frame changes that occur
            when an edge flips (see the documentation for ``Matching.add_edge`` for more information).
        weight : float
            Present only if `return_weight==True`.
            The sum of the weights of the edges in the minimum-weight perfect
            matching.

        Raises
        ------
        ValueError
            If there is no error consistent with the provided syndrome. Occurs if the syndrome has odd parity in the
            support of a connected component without a boundary.

        Examples
        --------
        >>> import pymatching
        >>> import numpy as np
        >>> check_matrix = np.array([[1, 1, 0, 0, 0],
        ...               [0, 1, 1, 0, 0],
        ...               [0, 0, 1, 1, 0],
        ...               [0, 0, 0, 1, 1]])
        >>> m = pymatching.Matching(check_matrix)
        >>> z = np.array([0, 1, 0, 0])
        >>> m.decode(z)
        array([1, 1, 0, 0, 0], dtype=uint8)

        Each bit in the correction provided by ``Matching.decode`` corresponds to a
        fault_ids. The index of a bit in a correction corresponds to its fault_ids.
        For example, here an error on edge (0, 1) flips fault_ids 2 and 3, as
        inferred by the minimum-weight correction:

        >>> import pymatching
        >>> m = pymatching.Matching()
        >>> m.add_edge(0, 1, fault_ids={2, 3})
        >>> m.add_edge(1, 2, fault_ids=1)
        >>> m.add_edge(2, 0, fault_ids=0)
        >>> m.decode([1, 1, 0])
        array([0, 0, 1, 1], dtype=uint8)

        To decode with a phenomenological noise model (qubits and measurements both suffering
        bit-flip errors), you can provide a check matrix and number of syndrome repetitions to
        construct a matching graph with a time dimension (where nodes in consecutive time steps
        are connected by an edge), and then decode with a 2D syndrome
        (dimension 0 is space, dimension 1 is time):

        >>> import pymatching
        >>> import numpy as np
        >>> np.random.seed(0)
        >>> check_matrix = np.array([[1, 1, 0, 0],
        ...               [0, 1, 1, 0],
        ...               [0, 0, 1, 1]])
        >>> m = pymatching.Matching(check_matrix, repetitions=5)
        >>> data_qubit_noise = (np.random.rand(4, 5) < 0.1).astype(np.uint8)
        >>> print(data_qubit_noise)
        [[0 0 0 0 0]
         [0 0 0 0 0]
         [0 0 0 0 1]
         [1 1 0 0 0]]
        >>> cumulative_noise = (np.cumsum(data_qubit_noise, 1) % 2).astype(np.uint8)
        >>> syndrome = check_matrix@cumulative_noise % 2
        >>> print(syndrome)
        [[0 0 0 0 0]
         [0 0 0 0 1]
         [1 0 0 0 1]]
        >>> syndrome[:,:-1] ^= (np.random.rand(3, 4) < 0.1).astype(np.uint8)
        >>> # Take the parity of consecutive timesteps to construct a difference syndrome:
        >>> syndrome[:,1:] = syndrome[:,:-1] ^ syndrome[:,1:]
        >>> m.decode(syndrome)
        array([0, 0, 1, 0], dtype=uint8)
        """

        if _legacy_num_neighbours is not None:
            warnings.warn("The ``num_neighbours`` argument no longer has any effect in PyMatching v2.0.0 or later,"
                          " since it introduced an approximation that is no longer relevant or necessary. Providing "
                          "``num_neighbours`` as the second positional argument to ``Matching.decode`` will raise an "
                          "exception in a future version of PyMatching", DeprecationWarning, stacklevel=2)
        if _legacy_return_weight is not None:
            warnings.warn("The ``return_weights`` argument was provided as a positional argument, but in a future "
                          "version of PyMatching, it will be required to provide ``return_weights`` as a keyword "
                          "argument.", DeprecationWarning, stacklevel=2)
            return_weight = _legacy_return_weight
        detection_events = self._syndrome_array_to_detection_events(z)
        correction, weight = self._matching_graph.decode(detection_events)
        if return_weight:
            return correction, weight
        else:
            return correction

    def decode_batch(
            self,
            shots: np.ndarray,
            *,
            return_weights: bool = False,
            bit_packed_shots: bool = False,
            bit_packed_predictions: bool = False) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
        """
        Decode from a 2D `shots` array containing a batch of syndrome measurements. A faster
        alternative to using `pymatching.Matching.decode` and iterating over the shots in Python.

        Parameters
        ----------
        shots : np.ndarray
            A 2D numpy array of shots to decode, of `dtype=np.uint8`.

            If `bit_packed_shots==False`, then
            `shots` should have shape `shots.shape=(num_shots, syndrome_length)`, where `num_shots` is the
            number of shots (samples), and `syndrome_length` is the length of the binary syndrome vector to be
            decoded for each shot. If `len(self.boundary)==0` (e.g. if there is no boundary, or only a virtual
            boundary node, the default when loading from stim) then `syndrome_length=self.num_detectors`.
            However, `syndrome_length` is permitted to be as high as `self.num_nodes` in case the graph contains
            detectors nodes with an index larger than `self.num_detectors-1` (when `len(self.boundary)>0`).

            If `bit_packed_shots==True` then `shots` should have shape
            `shots.shape=(num_shots, math.ceil(syndrome_length / 8))`. Bit packing should be done using little endian
            order on the last axis (like ``np.packbits(data, bitorder='little', axis=1)``), so that the bit for
            detection event `m` in shot `s` can be found at ``(dets[s, m // 8] >> (m % 8)) & 1``.
        return_weights : bool
            If True, then also return a numpy array containing the weights of the solutions for all the shots.
            By default, False.
        bit_packed_shots : bool
            Set to `True` to provide `shots` as a bit-packed array, such that the bit for
            detection event `m` in shot `s` can be found at ``(dets[s, m // 8] >> (m % 8)) & 1``.
        bit_packed_predictions : bool
            Set to `True` if the returned predictions should be bit-packed, with the bit for fault id `m` in
            shot `s` in ``(obs[s, m // 8] >> (m % 8)) & 1``

        Returns
        -------
        predictions: np.ndarray
            The batch of predictions output by the decoder, a binary numpy array of `dtype=np.uint8` and with shape
            `predictions.shape=(num_shots, self.num_fault_ids)`. `predictions[i, j]=1` iff the decoder predicts that
            fault id `j` was flipped in the shot `i`.
        weights: np.ndarray
            The weights of the MWPM solutions, a numpy array of `dtype=float`. `weights[i]` is the weight of the
            MWPM solution in shot `i`.

        Examples
        --------
        >>> import pymatching
        >>> import stim
        >>> circuit = stim.Circuit.generated("surface_code:rotated_memory_x",
        ...                                  distance=5,
        ...                                  rounds=5,
        ...                                  after_clifford_depolarization=0.005)
        >>> model = circuit.detector_error_model(decompose_errors=True)
        >>> matching = pymatching.Matching.from_detector_error_model(model)
        >>> sampler = circuit.compile_detector_sampler()
        >>> syndrome, actual_observables = sampler.sample(shots=10000, separate_observables=True)
        >>> syndrome.shape
        (10000, 120)
        >>> actual_observables.shape
        (10000, 1)
        >>> predicted_observables = matching.decode_batch(syndrome)
        >>> predicted_observables.shape
        (10000, 1)
        >>> num_errors = np.sum(np.any(predicted_observables != actual_observables, axis=1))

        We can also decode bit-packed shots, and return bit-packed predictions:
        >>> import pymatching
        >>> import stim
        >>> circuit = stim.Circuit.generated("surface_code:rotated_memory_x",
        ...                                  distance=5,
        ...                                  rounds=5,
        ...                                  after_clifford_depolarization=0.005)
        >>> model = circuit.detector_error_model(decompose_errors=True)
        >>> matching = pymatching.Matching.from_detector_error_model(model)
        >>> sampler = circuit.compile_detector_sampler()
        >>> syndrome, actual_observables = sampler.sample(shots=10000, separate_observables=True, bit_packed=True)
        >>> syndrome.shape
        (10000, 15)
        >>> actual_observables.shape
        (10000, 1)
        >>> predicted_observables = matching.decode_batch(syndrome, bit_packed_shots=True, bit_packed_predictions=True)
        >>> predicted_observables.shape
        (10000, 1)
        >>> num_errors = np.sum(np.any(predicted_observables != actual_observables, axis=1))
        """
        predictions, weights = self._matching_graph.decode_batch(
            shots,
            bit_packed_predictions=bit_packed_predictions,
            bit_packed_shots=bit_packed_shots
        )
        if return_weights:
            return predictions, weights
        else:
            return predictions

    def decode_to_edges_array(self,
                              syndrome: Union[np.ndarray, List[bool], List[int]]
                              ) -> np.ndarray:
        """
        Decode the syndrome `syndrome` using minimum-weight perfect matching, returning the edges in the
        solution, given as pairs of detector node indices in a numpy array.

        Parameters
        ----------
        syndrome : numpy.ndarray
            A binary syndrome vector to decode. The number of elements in
            `syndrome` should equal the number of nodes in the matching graph. If
            `syndrome` is a 1D array, then `syndrome[i]` is the syndrome at node `i` of
            the matching graph. If `syndrome` is 2D then `syndrome[i,j]` is the difference
            (modulo 2) between the (noisy) measurement of stabiliser `i` in time
            step `j+1` and time step `j` (for the case where the matching graph is
            constructed from a check matrix with `repetitions>1`).

        Returns
        -------
        numpy.ndarray
            A 2D array `edges` giving the edges in the matching solution as pairs of detector nodes (or as a detector
            node and the boundary, for a boundary edge). If there are `num_predicted_edges` edges then the shape of
            `edges` is `edges.shape=(num_predicted_edges, 2)`, and edge `i` is between detector node `edges[i, 0]`
            and detector node `edges[i, 1]`. For a boundary edge `i` between a detector node `k` and the boundary
            (either a boundary node or the virtual boundary node), then `pairs[i,0]` is `k`, and `pairs[i,1]=-1`
            denotes the boundary (the boundary is always denoted by -1 and is always in the second column).

        Examples
        --------
        >>> import pymatching
        >>> m = pymatching.Matching()
        >>> m.add_boundary_edge(0)
        >>> m.add_edge(0, 1)
        >>> m.add_edge(1, 2)
        >>> m.add_edge(2, 3)
        >>> m.add_edge(3, 4)
        >>> m.add_edge(4, 5)
        >>> m.add_edge(5, 6)
        >>> edges = m.decode_to_edges_array([0, 1, 0, 0, 1, 0, 1])
        >>> print(edges)
        [[ 0  1]
         [ 0 -1]
         [ 5  4]
         [ 5  6]]
        """
        detection_events = self._syndrome_array_to_detection_events(syndrome)
        return self._matching_graph.decode_to_edges_array(detection_events)

    def decode_to_matched_dets_array(self,
                                     syndrome: Union[np.ndarray, List[bool], List[int]]
                                     ) -> np.ndarray:
        """
        Decode the syndrome `syndrome` using minimum-weight perfect matching, returning the pairs of
        matched detection events (or detection events matched to the boundary) as a 2D numpy array.
        Each pair of matched detection events returned by this method corresponds to a shortest path
        between the detection events in the solution to the problem: if you instead want the set of
        all edges in the solution (pairs of detector nodes), use `Matching.decode_to_edges` instead.
        Note that, unlike `Matching.decode`, `Matching.decode_batch` and `Matching.decode_to_edges_array`,
        this method currently only supports non-negative edge weights.

        Parameters
        ----------
        syndrome : numpy.ndarray
            A binary syndrome vector to decode. The number of elements in
            `syndrome` should equal the number of nodes in the matching graph. If
            `syndrome` is a 1D array, then `syndrome[i]` is the syndrome at node `i` of
            the matching graph. If `syndrome` is 2D then `syndrome[i,j]` is the difference
            (modulo 2) between the (noisy) measurement of stabiliser `i` in time
            step `j+1` and time step `j` (for the case where the matching graph is
            constructed from a check matrix with `repetitions>1`).

        Returns
        -------
        numpy.ndarray
            An 2D array `pairs` giving the endpoints of the paths between detection events in the solution of the
            matching. If there are `num_paths` paths then the shape of `pairs` is `pairs.shape=(num_paths, 2)`, and
            path `i` starts at detection event `pairs[i,0]` and ends at detection event `pairs[i,1]`. For a path `i`
            connecting a detection event to the boundary (either a boundary node or the virtual boundary node), then
            `pairs[i,0]` is the index of the detection event, and `pairs[i,1]=-1` denotes the boundary.

        Examples
        --------
        >>> import pymatching
        >>> m = pymatching.Matching()
        >>> m.add_boundary_edge(0)
        >>> m.add_edge(0, 1)
        >>> m.add_edge(1, 2)
        >>> m.add_edge(2, 3)
        >>> m.add_edge(3, 4)
        >>> m.add_edge(4, 5)
        >>> m.add_edge(5, 6)
        >>> matched_dets = m.decode_to_matched_dets_array([0, 1, 0, 0, 1, 0, 1])
        >>> print(matched_dets)
        [[ 1 -1]
         [ 4  6]]
        """
        detection_events = self._syndrome_array_to_detection_events(syndrome)
        return self._matching_graph.decode_to_matched_detection_events_array(detection_events)

    def decode_to_matched_dets_dict(self,
                                    syndrome: Union[np.ndarray, List[bool], List[int]]
                                    ) -> Union[np.ndarray, Tuple[np.ndarray, int]]:
        """
        Decode the syndrome `syndrome` using minimum-weight perfect matching, returning a dictionary
        giving the detection event that each detection event was matched to (or None if it was matched
        to the boundary). Note that (unlike `Matching.decode`), this method currently only supports non-negative
        edge weights.

        Parameters
        ----------
        syndrome : numpy.ndarray
            A binary syndrome vector to decode. The number of elements in
            `syndrome` should equal the number of nodes in the matching graph. If
            `syndrome` is a 1D array, then `syndrome[i]` is the syndrome at node `i` of
            the matching graph. If `syndrome` is 2D then `syndrome[i,j]` is the difference
            (modulo 2) between the (noisy) measurement of stabiliser `i` in time
            step `j+1` and time step `j` (for the case where the matching graph is
            constructed from a check matrix with `repetitions>1`).

        Returns
        -------
        dict
            A dictionary `mate` giving the detection event that each detection event is matched to (or `None` if
            it is matched to the boundary). If detection event `i` is matched to detection event `j`, then
            `mate[i]=j`. If detection event `i` is matched to the boundary (either a boundary node or the virtual boundary
            node), then `mate[i]=None`.

        Examples
        --------
        >>> import pymatching
        >>> m = pymatching.Matching()
        >>> m.add_boundary_edge(0)
        >>> m.add_edge(0, 1)
        >>> m.add_edge(1, 2)
        >>> m.add_edge(2, 3)
        >>> m.add_edge(3, 4)
        >>> d = m.decode_to_matched_dets_dict([1, 0, 0, 1, 1])
        >>> d[3]
        4
        >>> d
        {0: None, 3: 4, 4: 3}
        """
        detection_events = self._syndrome_array_to_detection_events(syndrome)
        return self._matching_graph.decode_to_matched_detection_events_dict(detection_events)

    def draw(self) -> None:
        """Draw the matching graph using matplotlib
        Draws the matching graph as a matplotlib graph. Detector nodes are
        filled grey and boundary nodes are filled white. The line thickness of each
        edge is determined from its weight (with min and max thicknesses of 0.2 pts
        and 2 pts respectively).
        Each node is labelled with its id/index, and each edge is labelled with its `fault_ids`.
        Note that you may need to call `plt.figure()` before and `plt.show()` after calling
        this function.
        """
        # Ignore matplotlib deprecation warnings from networkx.draw_networkx
        warnings.filterwarnings("ignore", category=matplotlib.MatplotlibDeprecationWarning)
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        G = self.to_networkx()
        pos = nx.spectral_layout(G, weight=None)
        c = "#bfbfbf"
        ncolors = ['w' if n[1]['is_boundary'] else c for n in G.nodes(data=True)]
        nx.draw_networkx_nodes(G, pos=pos, node_color=ncolors, edgecolors=c)
        nx.draw_networkx_labels(G, pos=pos)
        weights = np.array([e[2]['weight'] for e in G.edges(data=True)])
        normalised_weights = 0.2 + 2 * weights / np.max(weights)
        nx.draw_networkx_edges(G, pos=pos, width=normalised_weights)

        def qid_to_str(qid):
            if len(qid) == 0:
                return ""
            elif len(qid) == 1:
                return str(qid.pop())
            else:
                return str(qid)

        edge_labels = {(s, t): qid_to_str(d['fault_ids']) for (s, t, d) in G.edges(data=True)}
        nx.draw_networkx_edge_labels(G, pos=pos, edge_labels=edge_labels)

    def __repr__(self) -> str:
        m = self.num_detectors
        b = len(self.boundary)
        e = self._matching_graph.get_num_edges()
        return "<pymatching.Matching object with " \
               "{} detector{}, " \
               "{} boundary node{}, " \
               "and {} edge{}>".format(
                m, 's' if m != 1 else '', b, 's' if b != 1 else '',
                e, 's' if e != 1 else '')

    def add_edge(
            self,
            node1: int,
            node2: int,
            fault_ids: Union[int, Set[int]] = None,
            weight: float = 1.0,
            error_probability: float = None,
            *,
            merge_strategy: str = "disallow",
            **kwargs
    ) -> None:
        """
        Add an edge to the matching graph

        Parameters
        ----------
        node1: int
            The index of node1 in the new edge (node1, node2)
        node2: int
            The index of node2 in the new edge (node1, node2)
        fault_ids: set[int] or int, optional
            The indices of any self-inverse faults which are flipped when the edge is flipped, and which should be tracked.
            This could correspond to the IDs of physical Pauli errors that occur when this
            edge flips (physical frame changes). Alternatively,
            this attribute can be used to store the IDs of any logical observables that are
            flipped when an error occurs on an edge (logical frame changes). In earlier versions of PyMatching, this
            attribute was instead named `qubit_id` (since for CSS codes and physical frame changes, there can be
            a one-to-one correspondence between each fault ID and physical qubit ID). For backward
            compatibility, `qubit_id` can still be used instead of `fault_ids` as a keyword argument.
            By default None
        weight: float, optional
            The weight of the edge. The weight can be positive or negative, but its absolute value cannot exceed
            the maximum absolute edge weight of 2**24-1=16,777,215. If the absolute value of the weight exceeds this
            value, the edge will not be added to the graph and a warning will be raised. By default 1.0
        error_probability: float, optional
            The probability that the edge is flipped. This is used by the `add_noise()` method
            to sample from the distribution defined by the matching graph (in which each edge
            is flipped independently with the corresponding `error_probability`). By default None
        merge_strategy: str, optional
            Which strategy to use if the edge (`node1`, `node2`) is already in the graph. The available options
            are "disallow", "independent", "smallest-weight", "keep-original" and "replace". "disallow" raises a
            `ValueError` if the edge (`node1`, `node2`) is already present. The "independent" strategy assumes that
            the existing edge (`node1`, `node2`) and the edge being added represent independent error mechanisms, and
            they are merged into a new edge with updated weights and error_probabilities accordingly (it is assumed
            that each weight represents the log-likelihood ratio log((1-p)/p) where p is the `error_probability` and
            where the natural logarithm is used. The fault_ids associated with the existing edge are kept only, since
            where the natural logarithm is used. The fault_ids associated with the existing edge are kept only, since
            the code has distance 2 if parallel edges have different fault_ids anyway). The "smallest-weight" strategy
            keeps only the new edge if it has a smaller weight than the existing edge, otherwise the graph is left
            unchanged. The "keep-original" strategy keeps only the existing edge, and ignores the edge being added.
            The "replace" strategy always keeps the edge being added, replacing the existing edge.
            By default, "disallow"

        Examples
        --------
        >>> import pymatching
        >>> m = pymatching.Matching()
        >>> m.add_edge(0, 1)
        >>> m.add_edge(1, 2)
        >>> print(m.num_edges)
        2
        >>> print(m.num_nodes)
        3
        >>> import math
        >>> m = pymatching.Matching()
        >>> m.add_edge(0, 1, fault_ids=2, weight=math.log((1-0.05)/0.05), error_probability=0.05)
        >>> m.add_edge(1, 2, fault_ids=0, weight=math.log((1-0.1)/0.1), error_probability=0.1)
        >>> m.add_edge(2, 0, fault_ids={1, 2}, weight=math.log((1-0.2)/0.2), error_probability=0.2)
        >>> m
        <pymatching.Matching object with 3 detectors, 0 boundary nodes, and 3 edges>
        >>> m = pymatching.Matching()
        >>> m.add_edge(0, 1, fault_ids=0, weight=2)
        >>> m.add_edge(0, 1, fault_ids=1, weight=1, merge_strategy="smallest-weight")
        >>> m.add_edge(0, 1, fault_ids=2, weight=3, merge_strategy="smallest-weight")
        >>> m.edges()
        [(0, 1, {'fault_ids': {1}, 'weight': 1.0, 'error_probability': -1.0})]
        """
        if fault_ids is not None and "qubit_id" in kwargs:
            raise ValueError("Both `fault_ids` and `qubit_id` were provided as arguments. Please "
                             "provide `fault_ids` instead of `qubit_id` as an argument, as use of `qubit_id` has "
                             "been deprecated.")
        if fault_ids is None and "qubit_id" in kwargs:
            fault_ids = kwargs["qubit_id"]
        if isinstance(fault_ids, (int, np.integer)):
            fault_ids = set() if fault_ids == -1 else {int(fault_ids)}
        fault_ids = set() if fault_ids is None else fault_ids
        error_probability = error_probability if error_probability is not None else -1
        self._matching_graph.add_edge(node1, node2, fault_ids, weight,
                                      error_probability, merge_strategy)

    def add_boundary_edge(
            self,
            node: int,
            fault_ids: Union[int, Set[int]] = None,
            weight: float = 1.0,
            error_probability: float = None,
            *,
            merge_strategy: str = "disallow",
            **kwargs
    ) -> None:
        """
        Add an edge connecting `node` to the boundary

        Parameters
        ----------
        node: int
            The index of the node to be connected to the boundary with a boundary edge
        fault_ids: set[int] or int, optional
            The IDs of any self-inverse faults which are flipped when the edge is flipped, and which should be tracked.
            This could correspond to the IDs of physical Pauli errors that occur when this
            edge flips (physical frame changes). Alternatively,
            this attribute can be used to store the IDs of any logical observables that are
            flipped when an error occurs on an edge (logical frame changes). By default None
        weight: float, optional
            The weight of the edge. The weight can be positive or negative, but its absolute value cannot exceed
            the maximum absolute edge weight of 2**24-1=16,777,215. If the absolute value of the weight exceeds this
            value, the edge will not be added to the graph and a warning will be raised. By default 1.0
        error_probability: float, optional
            The probability that the edge is flipped. This is used by the `add_noise()` method
            to sample from the distribution defined by the matching graph (in which each edge
            is flipped independently with the corresponding `error_probability`). By default None
        merge_strategy: str, optional
            Which strategy to use if the edge (`node1`, `node2`) is already in the graph. The available options
            are "disallow", "independent", "smallest-weight", "keep-original" and "replace". "disallow" raises a
            `ValueError` if the edge (`node1`, `node2`) is already present. The "independent" strategy assumes that
            the existing edge (`node1`, `node2`) and the edge being added represent independent error mechanisms, and
            they are merged into a new edge with updated weights and error_probabilities accordingly (it is assumed
            that each weight represents the log-likelihood ratio log((1-p)/p) where p is the `error_probability` and
            where the natural logarithm is used. The fault_ids associated with the existing edge are kept only, since
            where the natural logarithm is used. The fault_ids associated with the existing edge are kept only, since
            the code has distance 2 if parallel edges have different fault_ids anyway). The "smallest-weight" strategy
            keeps only the new edge if it has a smaller weight than the existing edge, otherwise the graph is left
            unchanged. The "keep-original" strategy keeps only the existing edge, and ignores the edge being added.
            The "replace" strategy always keeps the edge being added, replacing the existing edge.
            By default, "disallow"

        Examples
        --------
        >>> import pymatching
        >>> m = pymatching.Matching()
        >>> m.add_boundary_edge(0)
        >>> m.add_edge(0, 1)
        >>> print(m.num_edges)
        2
        >>> print(m.num_nodes)
        2
        >>> import math
        >>> m = pymatching.Matching()
        >>> m.add_boundary_edge(0, fault_ids={0}, weight=math.log((1-0.05)/0.05), error_probability=0.05)
        >>> m.add_edge(0, 1, fault_ids={1}, weight=math.log((1-0.1)/0.1), error_probability=0.1)
        >>> m.add_boundary_edge(1, fault_ids={2}, weight=math.log((1-0.2)/0.2), error_probability=0.2)
        >>> m
        <pymatching.Matching object with 2 detectors, 0 boundary nodes, and 3 edges>
        >>> m = pymatching.Matching()
        >>> m.add_boundary_edge(0, fault_ids=0, weight=2)
        >>> m.add_boundary_edge(0, fault_ids=1, weight=1, merge_strategy="smallest-weight")
        >>> m.add_boundary_edge(0, fault_ids=2, weight=3, merge_strategy="smallest-weight")
        >>> m.edges()
        [(0, None, {'fault_ids': {1}, 'weight': 1.0, 'error_probability': -1.0})]
        >>> m.boundary  # Using Matching.add_boundary_edge, no boundary nodes are added (the boundary is a virtual node)
        set()
        """
        if isinstance(fault_ids, (int, np.integer)):
            fault_ids = set() if fault_ids == -1 else {int(fault_ids)}
        fault_ids = set() if fault_ids is None else fault_ids
        error_probability = error_probability if error_probability is not None else -1
        self._matching_graph.add_boundary_edge(node, fault_ids, weight,
                                               error_probability, merge_strategy)

    def has_edge(self, node1: int, node2: int) -> bool:
        """
        Returns True if edge `(node1, node2)` is in the graph.

        Parameters
        ----------
        node1: int
            The index of the first node
        node2: int
            The index of the second node

        Returns
        -------
        bool
            True if the edge `(node1, node2)` is in the graph, otherwise False.
        """
        return self._matching_graph.has_edge(node1, node2)

    def has_boundary_edge(self, node: int) -> bool:
        """
        Returns True if the boundary edge `(node,)` is in the graph. Note: this method does
        not check if `node` is connected to a boundary node in `Matching.boundary`; it only
        checks if `node` is connected to the virtual boundary node (i.e. whether there is a boundary
        edge `(node,)` present).

        Parameters
        ----------
        node: int
            The index of the node

        Returns
        -------
        bool
            True if the boundary edge `(node,)` is present, otherwise False.

        """
        return self._matching_graph.has_boundary_edge(node)

    def get_edge_data(self, node1: int, node2: int) -> Dict[str, Union[Set[int], float]]:
        """
        Returns the edge data associated with the edge `(node1, node2)`.

        Parameters
        ----------
        node1: int
            The index of the first node
        node2: int
            The index of the second node

        Returns
        -------
        dict
            A dictionary with keys `fault_ids`, `weight` and `error_probability`, and values giving the respective
            edge attributes
        """
        return self._matching_graph.get_edge_data(node1, node2)

    def get_boundary_edge_data(self, node: int) -> Dict[str, Union[Set[int], float]]:
        """
        Returns the edge data associated with the boundary edge `(node,)`.

        Parameters
        ----------
        node: int
            The index of the node

        Returns
        -------
        dict
            A dictionary with keys `fault_ids`, `weight` and `error_probability`, and values giving the respective
            boundary edge attributes
        """
        return self._matching_graph.get_boundary_edge_data(node)

    def edges(self) -> List[Tuple[int, Optional[int], Dict]]:
        """Edges of the matching graph
        Returns a list of edges of the matching graph. Each edge is a
        tuple ``(source, target, attr)`` where `source` and `target` are ints corresponding to the
        indices of the source and target nodes, and `attr` is a dictionary containing the
        attributes of the edge.
        The dictionary `attr` has keys `fault_ids` (a set of ints), `weight` (the weight of the edge,
        set to 1.0 if not specified), and `error_probability`
        (the error probability of the edge, set to -1 if not specified).

        Returns
        -------
        List of (int, int, dict) tuples
            A list of edges of the matching graph
        """
        return self._matching_graph.get_edges()

    @staticmethod
    def from_check_matrix(
            check_matrix: Union[csc_matrix, spmatrix, np.ndarray, List[List[int]]],
            weights: Union[float, np.ndarray, List[float]] = None,
            error_probabilities: Union[float, np.ndarray, List[float]] = None,
            repetitions: int = None,
            timelike_weights: Union[float, np.ndarray, List[float]] = None,
            measurement_error_probabilities: Union[float, np.ndarray, List[float]] = None,
            *,
            faults_matrix: Union[csc_matrix, spmatrix, np.ndarray, List[List[int]]] = None,
            merge_strategy: str = "smallest-weight",
            use_virtual_boundary_node: bool = False,
            **kwargs
    ) -> 'pymatching.Matching':
        r"""
        Load a matching graph from a check matrix

        Parameters
        ----------
        check_matrix : `scipy.csc_matrix` or `numpy.ndarray` or List[List[int]]
            The quantum code to be decoded with minimum-weight perfect
            matching, given as a binary check matrix (scipy sparse
            matrix or numpy.ndarray)
        weights : float or numpy.ndarray, optional
            If `check_matrix` is given as a scipy or numpy array, `weights` gives the weights
            of edges in the matching graph corresponding to columns of `check_matrix`.
            If `weights` is a numpy.ndarray, it should be a 1D array with length
            equal to `check_matrix.shape[1]`. If weights is a float, it is used as the weight for all
            edges corresponding to columns of `check_matrix`. By default None, in which case
            all weights are set to 1.0
            This argument was renamed from `spacelike_weights` in PyMatching v2.0, but
            `spacelike_weights` is still accepted in place of `weights` for backward compatibility.
        error_probabilities : float or numpy.ndarray, optional
            The probabilities with which an error occurs on each edge associated with a
            column of check_matrix. If a
            single float is given, the same error probability is used for each
            column. If a numpy.ndarray of floats is given, it must have a
            length equal to the number of columns in check_matrix. This parameter is only
            needed for the Matching.add_noise method, and not for decoding.
            By default None
        repetitions : int, optional
            The number of times the stabiliser measurements are repeated, if
            the measurements are noisy. By default None
        timelike_weights : float or numpy.ndarray, optional
            If `repetitions>1`, `timelike_weights` gives the weight of
            timelike edges. If a float is given, all timelike edges weights are set to
            the same value. If a numpy array of size `(check_matrix.shape[0],)` is given, the
            edge weight for each vertical timelike edge associated with the `i`th check (row)
            of `check_matrix` is set to `timelike_weights[i]`. By default None, in which case all
            timelike weights are set to 1.0
        measurement_error_probabilities : float or numpy.ndarray, optional
            If `repetitions>1`, gives the probability of a measurement
            error to be used for the add_noise method. If a float is given, all measurement
            errors are set to the same value. If a numpy array of size `(check_matrix.shape[0],)` is given,
            the error probability for each vertical timelike edge associated with the `i`th check
            (row) of `check_matrix` is set to `measurement_error_probabilities[i]`. This argument can also be
            given using the keyword argument `measurement_error_probability` to maintain backward
            compatibility with previous versions of Pymatching. By default None
        faults_matrix: `scipy.csc_matrix` or `numpy.ndarray` or List[List[int]], optional
            A binary array of faults, which can be used to set the `fault_ids` for each edge in the
            constructed matching graph. The `fault_ids` attribute of the edge corresponding to column
            `j` of `check_matrix` includes fault id `i` if and only if `faults[i,j]==1`. Therefore, the number
            of columns in `faults` must match the number of columns in `check_matrix`. By default, `faults` is just
            set to the identity matrix, in which case the edge corresponding to column `j` of `check_matrix`
            has `fault_ids={j}`. As an example, if `check_matrix` corresponds to the X check matrix of
            a CSS stabiliser code, then you could set `faults` to the X logical operators: in this case
            the output of `Matching.decode` will be a binary array `correction` where `correction[i]==1`
            if the decoder predicts that the logical operator corresponding to row `i` of `faults` was flipped,
            given the observed syndrome.
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
            will be handled by adding an edge `(i, check_matrix.shape[0])`, and marking the node `check_matrix.shape[0]` as a boundary node with
            `Matching.set_boundary(check_matrix.shape[0])`. The resulting graph will contain `check_matrix.shape[0]+1` nodes, the largest of
            which is the boundary node. If `use_virtual_boundary_node=True` then instead the boundary is a virtual node, and
            this column is handled with `Matching.add_boundary_edge(i, ...)`. The resulting graph will contain `check_matrix.shape[0]`
            nodes, and there is no boundary node. Both options are handled identically by the decoder, although
            `use_virtual_boundary_node=True` is recommended since it is simpler (with a one-to-one correspondence between
            nodes and rows of check_matrix), and is also slightly more efficient. By default, False (for backward compatibility)

        Examples
        --------

        >>> import pymatching
        >>> m = pymatching.Matching.from_check_matrix([[1, 1, 0, 0], [0, 1, 1, 0], [0, 0, 1, 1]])
        >>> m
        <pymatching.Matching object with 3 detectors, 1 boundary node, and 4 edges>

        Matching objects can also be initialised from a sparse scipy matrix:

        >>> import pymatching
        >>> from scipy.sparse import csc_matrix
        >>> check_matrix = csc_matrix([[1, 1, 0], [0, 1, 1]])
        >>> m = pymatching.Matching.from_check_matrix(check_matrix)
        >>> m
        <pymatching.Matching object with 2 detectors, 1 boundary node, and 3 edges>

        """
        m = pymatching.Matching()
        m.load_from_check_matrix(
            check_matrix=check_matrix,
            weights=weights,
            error_probabilities=error_probabilities,
            repetitions=repetitions,
            timelike_weights=timelike_weights,
            measurement_error_probabilities=measurement_error_probabilities,
            faults_matrix=faults_matrix,
            merge_strategy=merge_strategy,
            use_virtual_boundary_node=use_virtual_boundary_node,
            **kwargs
        )
        return m

    def load_from_check_matrix(self,
                               check_matrix: Union[csc_matrix, spmatrix, np.ndarray, List[List[int]]] = None,
                               weights: Union[float, np.ndarray, List[float]] = None,
                               error_probabilities: Union[float, np.ndarray, List[float]] = None,
                               repetitions: int = None,
                               timelike_weights: Union[float, np.ndarray, List[float]] = None,
                               measurement_error_probabilities: Union[float, np.ndarray, List[float]] = None,
                               *,
                               faults_matrix: Union[csc_matrix, spmatrix, np.ndarray, List[List[int]]] = None,
                               merge_strategy: str = "smallest-weight",
                               use_virtual_boundary_node: bool = False,
                               **kwargs
                               ) -> None:
        """
        Load a matching graph from a check matrix

        Parameters
        ----------
        check_matrix : `scipy.csc_matrix` or `numpy.ndarray` or List[List[int]]
            The quantum code to be decoded with minimum-weight perfect
            matching, given as a binary check matrix (scipy sparse
            matrix or numpy.ndarray)
        weights : float or numpy.ndarray, optional
            If `check_matrix` is given as a scipy or numpy array, `weights` gives the weights
            of edges in the matching graph corresponding to columns of `check_matrix`.
            If `weights` is a numpy.ndarray, it should be a 1D array with length
            equal to `check_matrix.shape[1]`. If weights is a float, it is used as the weight for all
            edges corresponding to columns of `check_matrix`. By default None, in which case
            all weights are set to 1.0
            This argument was renamed from `spacelike_weights` in PyMatching v2.0, but
            `spacelike_weights` is still accepted in place of `weights` for backward compatibility.
        error_probabilities : float or numpy.ndarray, optional
            The probabilities with which an error occurs on each edge associated with a
            column of check_matrix. If a
            single float is given, the same error probability is used for each
            column. If a numpy.ndarray of floats is given, it must have a
            length equal to the number of columns in check_matrix. This parameter is only
            needed for the Matching.add_noise method, and not for decoding.
            By default None
        repetitions : int, optional
            The number of times the stabiliser measurements are repeated, if
            the measurements are noisy. By default None
        timelike_weights : float or numpy.ndarray, optional
            If `repetitions>1`, `timelike_weights` gives the weight of
            timelike edges. If a float is given, all timelike edges weights are set to
            the same value. If a numpy array of size `(check_matrix.shape[0],)` is given, the
            edge weight for each vertical timelike edge associated with the `i`th check (row)
            of `check_matrix` is set to `timelike_weights[i]`. By default None, in which case all
            timelike weights are set to 1.0
        measurement_error_probabilities : float or numpy.ndarray, optional
            If `repetitions>1`, gives the probability of a measurement
            error to be used for the add_noise method. If a float is given, all measurement
            errors are set to the same value. If a numpy array of size `(check_matrix.shape[0],)` is given,
            the error probability for each vertical timelike edge associated with the `i`th check
            (row) of `check_matrix` is set to `measurement_error_probabilities[i]`. This argument can also be
            given using the keyword argument `measurement_error_probability` to maintain backward
            compatibility with previous versions of Pymatching. By default None
        faults_matrix: `scipy.csc_matrix` or `numpy.ndarray` or List[List[int]], optional
            A binary array of faults, which can be used to set the `fault_ids` for each edge in the
            constructed matching graph. The `fault_ids` attribute of the edge corresponding to column
            `j` of `check_matrix` includes fault id `i` if and only if `faults[i,j]==1`. Therefore, the number
            of columns in `faults` must match the number of columns in `check_matrix`. By default, `faults` is just
            set to the identity matrix, in which case the edge corresponding to column `j` of `check_matrix`
            has `fault_ids={j}`. As an example, if `check_matrix` corresponds to the X check matrix of
            a CSS stabiliser code, then you could set `faults` to the X logical operators: in this case
            the output of `Matching.decode` will be a binary array `correction` where `correction[i]==1`
            if the decoder predicts that the logical operator corresponding to row `i` of `faults` was flipped,
            given the observed syndrome.
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
            will be handled by adding an edge `(i, check_matrix.shape[0])`, and marking the node `check_matrix.shape[0]` as a boundary node with
            `Matching.set_boundary(check_matrix.shape[0])`. The resulting graph will contain `check_matrix.shape[0]+1` nodes, the largest of
            which is the boundary node. If `use_virtual_boundary_node=True` then instead the boundary is a virtual node, and
            this column is handled with `Matching.add_boundary_edge(i, ...)`. The resulting graph will contain `check_matrix.shape[0]`
            nodes, and there is no boundary node. Both options are handled identically by the decoder, although
            `use_virtual_boundary_node=True` is recommended since it is simpler (with a one-to-one correspondence between
            nodes and rows of check_matrix), and is also slightly more efficient. By default, False (for backward compatibility)

        Examples
        --------
        >>> import pymatching
        >>> m = pymatching.Matching()
        >>> m.load_from_check_matrix([[1, 1, 0, 0], [0, 1, 1, 0], [0, 0, 1, 1]])
        >>> m
        <pymatching.Matching object with 3 detectors, 1 boundary node, and 4 edges>

        Matching objects can also be initialised from a sparse scipy matrix:
        >>> import pymatching
        >>> from scipy.sparse import csc_matrix
        >>> check_matrix = csc_matrix([[1, 1, 0], [0, 1, 1]])
        >>> m = pymatching.Matching()
        >>> m.load_from_check_matrix(check_matrix)
        >>> m
        <pymatching.Matching object with 2 detectors, 1 boundary node, and 3 edges>
        """
        if check_matrix is None:
            check_matrix = kwargs.get("H", None)
            if check_matrix is None:
                raise ValueError("No check_matrix provided")

        if not isinstance(check_matrix, csc_matrix):
            try:
                check_matrix = csc_matrix(check_matrix)
            except TypeError:
                raise TypeError("`check_matrix` must be convertible to a `scipy.sparse.csc_matrix`")

        if faults_matrix is not None:
            try:
                faults_matrix = csc_matrix(faults_matrix)
            except TypeError:
                raise TypeError("`faults` must be convertible to `scipy.sparse.csc_matrix`")

        num_edges = check_matrix.shape[1]

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

        check_matrix.eliminate_zeros()

        repetitions = 1 if repetitions is None else repetitions

        if repetitions > 1:
            timelike_weights = 1.0 if timelike_weights is None else timelike_weights
            if isinstance(timelike_weights, (int, float, np.integer, np.floating)):
                timelike_weights = np.ones(check_matrix.shape[0], dtype=float) * timelike_weights
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
                p_meas = np.ones(check_matrix.shape[0], dtype=float) * p_meas
            elif isinstance(p_meas, (np.ndarray, list)):
                p_meas = np.array(p_meas, dtype=float)
            else:
                raise ValueError("measurement_error_probabilities should be a float or 1d numpy array")
        else:
            timelike_weights = None
            p_meas = None
        self._matching_graph = _cpp_pm.sparse_column_check_matrix_to_matching_graph(check_matrix, weights,
                                                                                    error_probabilities,
                                                                                    merge_strategy,
                                                                                    use_virtual_boundary_node,
                                                                                    repetitions,
                                                                                    timelike_weights, p_meas,
                                                                                    faults_matrix)

    @staticmethod
    def from_detector_error_model(model: 'stim.DetectorErrorModel') -> 'pymatching.Matching':
        """
        Constructs a `pymatching.Matching` object by loading from a `stim.DetectorErrorModel`.

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
        the `logical_observable` indices associated with the first added parallel edge are kept for the merged edge.
        If you are loading a `pymatching.Matching` graph from a DEM, you may be interested in
        using the sinter Python package for monte carlo sampling: https://pypi.org/project/sinter/.

        Parameters
        ----------
        model : stim.DetectorErrorModel
            A stim DetectorErrorModel, with all error mechanisms either graphlike, or decomposed into graphlike
            error mechanisms

        Returns
        -------
        pymatching.Matching
            A `pymatching.Matching` object representing the graphlike error mechanisms in `model`

        Examples
        --------
        >>> import stim
        >>> import pymatching
        >>> circuit = stim.Circuit.generated("surface_code:rotated_memory_x",
        ...                                  distance=5,
        ...                                  rounds=5,
        ...                                  after_clifford_depolarization=0.005)
        >>> model = circuit.detector_error_model(decompose_errors=True)
        >>> matching = pymatching.Matching.from_detector_error_model(model)
        >>> matching
        <pymatching.Matching object with 120 detectors, 0 boundary nodes, and 502 edges>
        """
        m = Matching()
        m._load_from_detector_error_model(model)
        return m

    @staticmethod
    def from_detector_error_model_file(dem_path: str) -> 'pymatching.Matching':
        """
        Construct a `pymatching.Matching` by loading from a stim DetectorErrorModel file path.

        Parameters
        ----------
        dem_path : str
            The path of the detector error model file

        Returns
        -------
        pymatching.Matching
            A `pymatching.Matching` object representing the graphlike error mechanisms in the stim DetectorErrorModel
            in the file `dem_path`
        """
        m = Matching()
        m._matching_graph = _cpp_pm.detector_error_model_file_to_matching_graph(dem_path)
        return m

    @staticmethod
    def from_stim_circuit(circuit: 'stim.Circuit') -> 'pymatching.Matching':
        """
        Constructs a `pymatching.Matching` object by loading from a `stim.Circuit`

        Parameters
        ----------
        circuit : stim.Circuit
            A stim circuit containing error mechanisms that are all either graphlike, or decomposable into
            graphlike error mechanisms

        Returns
        -------
        pymatching.Matching
            A `pymatching.Matching` object representing the graphlike error mechanisms in `circuit`, with any hyperedge
            error mechanisms decomposed into graphlike error mechanisms. Parallel edges are merged using
            `merge_strategy="independent"`.


        Examples
        --------
        >>> import stim
        >>> import pymatching
        >>> circuit = stim.Circuit.generated("surface_code:rotated_memory_x",
        ...                                  distance=5,
        ...                                  rounds=5,
        ...                                  after_clifford_depolarization=0.005)
        >>> matching = pymatching.Matching.from_stim_circuit(circuit)
        >>> matching
        <pymatching.Matching object with 120 detectors, 0 boundary nodes, and 502 edges>
        """
        try:
            import stim
        except ImportError:  # pragma no cover
            raise TypeError(
                f"`circuit` must be a `stim.Circuit. Instead, got: {type(circuit)}.`"
                "The 'stim' package also isn't installed and is required for this method. \n"
                "To install stim using pip, run `pip install stim`."
            )
        if not isinstance(circuit, stim.Circuit):
            raise TypeError(f"`circuit` must be a `stim.Circuit`. Instead, got {type(circuit)}")
        m = Matching()
        m._matching_graph = _cpp_pm.detector_error_model_to_matching_graph(
            str(circuit.detector_error_model(decompose_errors=True))
        )
        return m

    @staticmethod
    def from_stim_circuit_file(stim_circuit_path: str) -> 'pymatching.Matching':
        """
        Construct a `pymatching.Matching` by loading from a stim circuit file path.

        Parameters
        ----------
        stim_circuit_path : str
            The path of the stim circuit file

        Returns
        -------
        pymatching.Matching
            A `pymatching.Matching` object representing the graphlike error mechanisms in the stim circuit
            in the file `stim_circuit_path`, with any hyperedge error mechanisms decomposed into graphlike error
            mechanisms. Parallel edges are merged using `merge_strategy="independent"`.
        """
        m = Matching()
        m._matching_graph = _cpp_pm.stim_circuit_file_to_matching_graph(stim_circuit_path)
        return m

    def _load_from_detector_error_model(self, model: 'stim.DetectorErrorModel') -> None:
        try:
            import stim
        except ImportError:  # pragma no cover
            raise TypeError(
                f"`model` must be a `stim.DetectorErrorModel. Instead, got: {type(model)}.`"
                "The 'stim' package also isn't installed and is required for this method. \n"
                "To install stim using pip, run `pip install stim`."
            )
        if not isinstance(model, stim.DetectorErrorModel):
            raise TypeError(f"'model' must be `stim.DetectorErrorModel`. Instead, got: {type(model)}")
        self._matching_graph = _cpp_pm.detector_error_model_to_matching_graph(str(model))

    @staticmethod
    def from_networkx(graph: nx.Graph, *, min_num_fault_ids: int = None) -> 'pymatching.Matching':
        r"""
        Returns a new `pymatching.Matching` object from a NetworkX graph

        Parameters
        ----------
        graph : networkx.Graph
            Each edge in the NetworkX graph can have optional
            attributes ``fault_ids``, ``weight`` and ``error_probability``.
            ``fault_ids`` should be an int or a set of ints.
            Each fault id corresponds to a self-inverse fault that is flipped when the
            corresponding edge is flipped. These self-inverse faults could correspond to
            physical Pauli errors (physical frame changes)
            or to the logical observables that are flipped by the fault
            (a logical frame change, equivalent to an obersvable ID in an error instruction in a Stim
            detector error model).
            The `fault_ids` attribute determines how the solution is output via `pymatching.Matching.decode`:
            the binary `correction` array has length `pymatching.Matching.num_fault_ids`, and `correction[i]`
            is 1 if and only if an odd number of edges in the MWPM solution have `i` in their `fault_ids` attribute.
            The `fault_ids` attribute was previously named `qubit_id` in an
            earlier version of PyMatching, and `qubit_id` is still accepted instead of `fault_ids` in order
            to maintain backward compatibility.
            Each ``weight`` attribute should be a non-negative float. If
            every edge is assigned an error_probability between zero and one,
            then the ``add_noise`` method can be used to simulate noise and
            flip edges independently in the graph.
        min_num_fault_ids: int
            Sets the minimum number of fault ids in the matching graph. Let `max_id` be the maximum fault id assigned to
            any of the edges in the graph. Then setting this argument will ensure that
            `Matching.num_fault_ids=max(min_num_fault_ids, max_id)`. Note that `Matching.num_fault_ids` sets the length
            of the correction array output by `Matching.decode`.

        Examples
        --------
        >>> import pymatching
        >>> import networkx as nx
        >>> import math
        >>> g = nx.Graph()
        >>> g.add_edge(0, 1, fault_ids=0, weight=math.log((1-0.1)/0.1), error_probability=0.1)
        >>> g.add_edge(1, 2, fault_ids=1, weight=math.log((1-0.15)/0.15), error_probability=0.15)
        >>> g.nodes[0]['is_boundary'] = True
        >>> g.nodes[2]['is_boundary'] = True
        >>> m = pymatching.Matching.from_networkx(g)
        >>> m
        <pymatching.Matching object with 1 detector, 2 boundary nodes, and 2 edges>
        """
        m = Matching()
        m.load_from_networkx(
            graph=graph, min_num_fault_ids=min_num_fault_ids
        )
        return m

    def load_from_networkx(self, graph: nx.Graph, *, min_num_fault_ids: int = None) -> None:
        r"""
        Load a matching graph from a NetworkX graph into a `pymatching.Matching` object

        Parameters
        ----------
        graph : networkx.Graph
            Each edge in the NetworkX graph can have optional
            attributes ``fault_ids``, ``weight`` and ``error_probability``.
            ``fault_ids`` should be an int or a set of ints.
            Each fault id corresponds to a self-inverse fault that is flipped when the
            corresponding edge is flipped. These self-inverse faults could correspond to
            physical Pauli errors (physical frame changes)
            or to the logical observables that are flipped by the fault
            (a logical frame change, equivalent to an obersvable ID in an error instruction in a Stim
            detector error model).
            The `fault_ids` attribute determines how the solution is output via `pymatching.Matching.decode`:
            the binary `correction` array has length `pymatching.Matching.num_fault_ids`, and `correction[i]`
            is 1 if and only if an odd number of edges in the MWPM solution have `i` in their `fault_ids` attribute.
            The `fault_ids` attribute was previously named `qubit_id` in an
            earlier version of PyMatching, and `qubit_id` is still accepted instead of `fault_ids` in order
            to maintain backward compatibility.
            Each ``weight`` attribute should be a non-negative float. If
            every edge is assigned an error_probability between zero and one,
            then the ``add_noise`` method can be used to simulate noise and
            flip edges independently in the graph.
        min_num_fault_ids: int
            Sets the minimum number of fault ids in the matching graph. Let `max_id` be the maximum fault id assigned to
            any of the edges in the graph. Then setting this argument will ensure that
            `Matching.num_fault_ids=max(min_num_fault_ids, max_id)`. Note that `Matching.num_fault_ids` sets the length
            of the correction array output by `Matching.decode`.

        Examples
        --------
        >>> import pymatching
        >>> import networkx as nx
        >>> import math
        >>> g = nx.Graph()
        >>> g.add_edge(0, 1, fault_ids=0, weight=math.log((1-0.1)/0.1), error_probability=0.1)
        >>> g.add_edge(1, 2, fault_ids=1, weight=math.log((1-0.15)/0.15), error_probability=0.15)
        >>> g.nodes[0]['is_boundary'] = True
        >>> g.nodes[2]['is_boundary'] = True
        >>> m = pymatching.Matching(g)
        >>> m
        <pymatching.Matching object with 1 detector, 2 boundary nodes, and 2 edges>
        """

        if not isinstance(graph, nx.Graph):
            raise TypeError("G must be a NetworkX graph")
        boundary = {i for i, attr in graph.nodes(data=True)
                    if attr.get("is_boundary", False)}
        num_nodes = graph.number_of_nodes()
        all_fault_ids = set()
        num_fault_ids = 0 if min_num_fault_ids is None else min_num_fault_ids
        g = _cpp_pm.MatchingGraph(num_nodes, num_fault_ids)
        g.set_boundary(boundary)
        for (u, v, attr) in graph.edges(data=True):
            u, v = int(u), int(v)
            if "fault_ids" in attr and "qubit_id" in attr:
                raise ValueError("Both `fault_ids` and `qubit_id` were provided as edge attributes, however use "
                                 "of `qubit_id` has been deprecated in favour of `fault_ids`. Please only supply "
                                 "`fault_ids` as an edge attribute.")
            if "fault_ids" not in attr and "qubit_id" in attr:
                fault_ids = attr["qubit_id"]  # Still accept qubit_id as well for now
            else:
                fault_ids = attr.get("fault_ids", set())
            if isinstance(fault_ids, (int, np.integer)):
                fault_ids = {int(fault_ids)} if fault_ids != -1 else set()
            else:
                try:
                    fault_ids = set(fault_ids)
                    if not all(isinstance(q, (int, np.integer)) for q in fault_ids):
                        raise TypeError("fault_ids must be a set of ints, not {}".format(fault_ids))
                except TypeError:
                    raise TypeError(
                        "fault_ids property must be an int or a set of int"
                        " (or convertible to a set), not {}".format(fault_ids))
            all_fault_ids = all_fault_ids | fault_ids
            weight = attr.get("weight", 1)  # Default weight is 1 if not provided
            e_prob = attr.get("error_probability", -1)
            # Note: NetworkX graphs do not support parallel edges (merge strategy is redundant)
            g.add_edge(u, v, fault_ids, weight, e_prob, merge_strategy="smallest-weight")
        self._matching_graph = g

    def load_from_retworkx(self, graph: "rx.PyGraph", *, min_num_fault_ids: int = None) -> None:
        r"""
        Load a matching graph from a retworkX graph. This method is deprecated since the retworkx package has been
        renamed to rustworkx. Please use ``pymatching.Matching.load_from_rustworkx`` instead.
        """
        warnings.warn("`pymatching.Matching.load_from_retworkx` is now deprecated since the `retworkx` library has been "
                      "renamed to `rustworkx`. Please use `pymatching.Matching.load_from_rustworkx` instead.", DeprecationWarning, stacklevel=2)
        self.load_from_rustworkx(graph=graph, min_num_fault_ids=min_num_fault_ids)

    def load_from_rustworkx(self, graph: "rx.PyGraph", *, min_num_fault_ids: int = None) -> None:
        r"""
        Load a matching graph from a rustworkX graph

        Parameters
        ----------
        graph : rustworkx.PyGraph
            Each edge in the rustworkx graph can have dictionary payload with keys
            ``fault_ids``, ``weight`` and ``error_probability``. ``fault_ids`` should be
            an int or a set of ints. Each fault id corresponds to a self-inverse fault
            that is flipped when the corresponding edge is flipped. These self-inverse
            faults could correspond to physical Pauli errors (physical frame changes)
            or to the logical observables that are flipped by the fault
            (a logical frame change, equivalent to an obersvable ID in an error instruction in a Stim
            detector error model). The `fault_ids` attribute was previously named `qubit_id` in an
            earlier version of PyMatching, and `qubit_id` is still accepted instead of `fault_ids` in order
            to maintain backward compatibility.
            Each ``weight`` attribute should be a non-negative float. If
            every edge is assigned an error_probability between zero and one,
            then the ``add_noise`` method can be used to simulate noise and
            flip edges independently in the graph.
        min_num_fault_ids: int
            Sets the minimum number of fault ids in the matching graph. Let `max_id` be the maximum fault id assigned to
            any of the edges in the graph. Then setting this argument will ensure that
            `Matching.num_fault_ids=max(min_num_fault_ids, max_id)`. Note that `Matching.num_fault_ids` sets the length
            of the correction array output by `Matching.decode`.

        Examples
        --------
        >>> import pymatching
        >>> import rustworkx as rx
        >>> import math
        >>> g = rx.PyGraph()
        >>> matching = g.add_nodes_from([{} for _ in range(3)])
        >>> edge_a =g.add_edge(0, 1, dict(fault_ids=0, weight=math.log((1-0.1)/0.1), error_probability=0.1))
        >>> edge_b = g.add_edge(1, 2, dict(fault_ids=1, weight=math.log((1-0.15)/0.15), error_probability=0.15))
        >>> g[0]['is_boundary'] = True
        >>> g[2]['is_boundary'] = True
        >>> m = pymatching.Matching(g)
        >>> m
        <pymatching.Matching object with 1 detector, 2 boundary nodes, and 2 edges>
        """
        try:
            import rustworkx as rx
        except ImportError:  # pragma no cover
            raise ImportError("rustworkx must be installed to use Matching.load_from_rustworkx")
        if not isinstance(graph, rx.PyGraph):
            raise TypeError("G must be a rustworkx graph")
        boundary = {i for i in graph.node_indices() if graph[i].get("is_boundary", False)}
        num_nodes = len(graph)
        num_fault_ids = 0 if min_num_fault_ids is None else min_num_fault_ids
        g = _cpp_pm.MatchingGraph(num_nodes, num_fault_ids)
        g.set_boundary(boundary)
        for (u, v, attr) in graph.weighted_edge_list():
            u, v = int(u), int(v)
            if "fault_ids" in attr and "qubit_id" in attr:
                raise ValueError("Both `fault_ids` and `qubit_id` were provided as edge attributes, however use "
                                 "of `qubit_id` has been deprecated in favour of `fault_ids`. Please only supply "
                                 "`fault_ids` as an edge attribute.")
            if "fault_ids" not in attr and "qubit_id" in attr:
                fault_ids = attr["qubit_id"]  # Still accept qubit_id as well for now
            else:
                fault_ids = attr.get("fault_ids", set())
            if isinstance(fault_ids, (int, np.integer)):
                fault_ids = {int(fault_ids)} if fault_ids != -1 else set()
            else:
                try:
                    fault_ids = set(fault_ids)
                    if not all(isinstance(q, (int, np.integer)) for q in fault_ids):
                        raise TypeError("fault_ids must be a set of ints, not {}".format(fault_ids))
                except TypeError:
                    raise TypeError(
                        "fault_ids property must be an int or a set of int"
                        " (or convertible to a set), not {}".format(fault_ids))
            weight = attr.get("weight", 1)  # Default weight is 1 if not provided
            e_prob = attr.get("error_probability", -1)
            # Note: rustworkx graphs do not support parallel edges (merge strategy is redundant)
            g.add_edge(u, v, fault_ids, weight, e_prob, merge_strategy="smallest-weight")
        self._matching_graph = g

    def to_networkx(self) -> nx.Graph:
        """Convert to NetworkX graph
        Returns a NetworkX graph corresponding to the matching graph. Each edge
        has attributes `fault_ids`, `weight` and `error_probability` and each node has
        the attribute `is_boundary`.

        Returns
        -------
        NetworkX.Graph
            NetworkX Graph corresponding to the matching graph
        """
        graph = nx.Graph()
        num_nodes = self.num_nodes
        has_virtual_boundary = False
        for u, v, data in self.edges():
            if v is None:
                graph.add_edge(u, num_nodes, **data)
                has_virtual_boundary = True
            else:
                graph.add_edge(u, v, **data)
        boundary = self.boundary
        for i in graph.nodes:
            is_boundary = i in boundary
            graph.nodes[i]['is_boundary'] = is_boundary
        if has_virtual_boundary:
            graph.nodes[num_nodes]['is_boundary'] = True
        return graph

    def to_retworkx(self) -> "rx.PyGraph":
        """Deprecated, use ``pymatching.Matching.to_rustworkx`` instead (since the `retworkx` package has been renamed to `rustworkx`).
        This method just calls ``pymatching.Matching.to_rustworkx`` and returns a ``rustworkx.PyGraph``, which is now just the preferred name for
         ``retworkx.PyGraph``. Note that in the future, only the `rustworkx` package name will be supported,
         see: https://pypi.org/project/retworkx/.
        """
        warnings.warn("`pymatching.Matching.to_retworkx` is now deprecated since the `retworkx` library has been "
                      "renamed to `rustworkx`. Please use `pymatching.Matching.to_rustworkx` instead.", DeprecationWarning, stacklevel=2)
        return self.to_rustworkx()

    def to_rustworkx(self) -> "rx.PyGraph":
        """Convert to rustworkx graph
        Returns a rustworkx graph object corresponding to the matching graph. Each edge
        payload is a ``dict`` with keys `fault_ids`, `weight` and `error_probability` and
        each node has a ``dict`` payload with the key ``is_boundary`` and the value is
        a boolean.

        Returns
        -------
        rustworkx.PyGraph
            rustworkx graph corresponding to the matching graph
        """
        try:
            import rustworkx as rx
        except ImportError:  # pragma no cover
            raise ImportError("rustworkx must be installed to use Matching.to_rustworkx.")

        graph = rx.PyGraph(multigraph=False)
        num_nodes = self.num_nodes
        has_virtual_boundary = False
        edges = []
        for u, v, data in self.edges():
            if v is None:
                edges.append((u, num_nodes, data))
                has_virtual_boundary = True
            else:
                edges.append((u, v, data))
        graph.add_nodes_from([{} for _ in range(num_nodes + has_virtual_boundary)])
        graph.extend_from_weighted_edge_list(edges)
        boundary = self.boundary
        for i in graph.node_indices():
            is_boundary = i in boundary
            graph[i]['is_boundary'] = is_boundary
        if has_virtual_boundary:
            graph[num_nodes]["is_boundary"] = True
        return graph

    def set_boundary_nodes(self, nodes: Set[int]) -> None:
        """
        Set boundary nodes in the matching graph. This defines the
        nodes in `nodes` to be boundary nodes.

        Parameters
        ----------
        nodes: set[int]
            The IDs of the nodes to be set as boundary nodes

        Examples
        --------
        >>> import pymatching
        >>> m = pymatching.Matching()
        >>> m.add_edge(0, 1)
        >>> m.add_edge(1, 2)
        >>> m.set_boundary_nodes({0, 2})
        >>> m.boundary
        {0, 2}
        >>> m
        <pymatching.Matching object with 1 detector, 2 boundary nodes, and 2 edges>
        """
        self._matching_graph.set_boundary(nodes)

    def ensure_num_fault_ids(self, min_num_fault_ids: int) -> None:
        """
        Set the minimum number of fault ids in the matching graph.

        Let `max_id` be the maximum fault id assigned to any of the edges in a `pymatching.Matching` graph `m`.
        Then setting `m.ensure_num_fault_ids(n)` will ensure that `Matching.num_fault_ids=max(n, max_id)`.
        Note that `Matching.num_fault_ids` sets the length of the correction array output by `Matching.decode`.

        Parameters
        ----------
        min_num_fault_ids: int
            The required minimum number of fault ids in the matching graph

        """
        self._matching_graph.set_min_num_observables(min_num_fault_ids)

    @property
    def num_fault_ids(self) -> int:
        """
        The number of fault IDs defined in the matching graph

        Returns
        -------
        int
            Number of fault IDs
        """
        return self._matching_graph.get_num_observables()

    @property
    def boundary(self) -> Set[int]:
        """Return the indices of the boundary nodes.
        Note that this property is a copy of the set of boundary nodes.
        In-place modification of the set Matching.boundary will not
        change the boundary nodes of the matching graph - boundary nodes should
        instead be set or updated using the `Matching.set_boundary_nodes` method.

        Returns
        -------
        set of int
            The indices of the boundary nodes
        """
        return self._matching_graph.get_boundary()

    @property
    def num_nodes(self) -> int:
        """
        The number of nodes in the matching graph

        Returns
        -------
        int
            The number of nodes
        """
        return self._matching_graph.get_num_nodes()

    @property
    def num_edges(self) -> int:
        """
        The number of edges in the matching graph

        Returns
        -------
        int
            The number of edges
        """
        return self._matching_graph.get_num_edges()

    @property
    def num_detectors(self) -> int:
        """
        The number of detectors in the matching graph. A
        detector is a node that can have a non-trivial syndrome
        (i.e. it is a node that is not a boundary node).

        Returns
        -------
        int
            The number of detectors
        """
        return self._matching_graph.get_num_detectors()
