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

from typing import Union, List

import numpy as np
import networkx as nx
import retworkx as rx
import scipy
import stim
from scipy.sparse import csc_matrix

from pymatching._cpp_pymatching import MatchingGraph, detector_error_model_to_matching_graph


def load_from_networkx(self, graph: nx.Graph) -> None:
    r"""
    Load a matching graph from a NetworkX graph
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
        detector error model). The `fault_ids` attribute was previously named `qubit_id` in an
        earlier version of PyMatching, and `qubit_id` is still accepted instead of `fault_ids` in order
        to maintain backward compatibility.
        Each ``weight`` attribute should be a non-negative float. If
        every edge is assigned an error_probability between zero and one,
        then the ``add_noise`` method can be used to simulate noise and
        flip edges independently in the graph.
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
    g = MatchingGraph(num_nodes)
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
                    raise ValueError("fault_ids must be a set of ints, not {}".format(fault_ids))
            except:
                raise ValueError(
                    "fault_ids property must be an int or a set of int" \
                    " (or convertible to a set), not {}".format(fault_ids))
        all_fault_ids = all_fault_ids | fault_ids
        weight = attr.get("weight", 1)  # Default weight is 1 if not provided
        e_prob = attr.get("error_probability", -1)
        # Note: NetworkX graphs do not support parallel edges (merge strategy is redundant)
        g.add_edge(u, v, fault_ids, weight, e_prob, merge_strategy="smallest-weight")
    self._matching_graph = g


def load_from_retworkx(self, graph: rx.PyGraph) -> None:
    r"""
    Load a matching graph from a retworkX graph
    Parameters
    ----------
    graph : retworkx.PyGraph
        Each edge in the retworkx graph can have dictionary payload with keys
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
    Examples
    --------
    >>> import pymatching
    >>> import retworkx as rx
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
    if not isinstance(graph, rx.PyGraph):
        raise TypeError("G must be a retworkx graph")
    boundary = {i for i in graph.node_indices() if graph[i].get("is_boundary", False)}
    num_nodes = len(graph)
    g = MatchingGraph(num_nodes)
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
                    raise ValueError("fault_ids must be a set of ints, not {}".format(fault_ids))
            except:
                raise ValueError(
                    "fault_ids property must be an int or a set of int" \
                    " (or convertible to a set), not {}".format(fault_ids))
        weight = attr.get("weight", 1)  # Default weight is 1 if not provided
        e_prob = attr.get("error_probability", -1)
        # Note: retworkx graphs do not support parallel edges (merge strategy is redundant)
        g.add_edge(u, v, fault_ids, weight, e_prob, merge_strategy="smallest-weight")
    self._matching_graph = g


def load_from_check_matrix(self,
                           H: Union[scipy.sparse.spmatrix, np.ndarray, List[List[int]]],
                           weights: Union[float, np.ndarray, List[float]] = None,
                           error_probabilities: Union[float, np.ndarray, List[float]] = None,
                           repetitions: int = None,
                           timelike_weights: Union[float, np.ndarray, List[float]] = None,
                           measurement_error_probabilities: Union[float, np.ndarray, List[float]] = None,
                           *,
                           merge_strategy: str = "smallest-weight",
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
        options are "disallow", "independent", "smallest-weight", "first-only" and "last-only". "disallow" raises a
        `ValueError` if the edge (`node1`, `node2`) is already present. The "independent" strategy assumes that
        the existing edge (`node1`, `node2`) and the edge being added represent independent error mechanisms, and
        they are merged into a new edge with updated weights and error_probabilities accordingly (it is assumed
        that each weight represents the log-likelihood ratio log((1-p)/p) where p is the `error_probability` and
        where the natural logarithm is used. The fault_ids associated with the existing edge are kept only, since
        the code has distance 2 if parallel edges have different fault_ids anyway). The "smallest-weight" strategy
        keeps only the new edge if it has a smaller weight than the existing edge, otherwise the graph is left
        unchanged. The "first-only" strategy keeps only the existing edge, and ignores the edge being added.
        The "last-only" strategy always keeps the edge being added, replacing the existing edge.
        By default, "smallest-weight"
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
    try:
        H = csc_matrix(H)
    except TypeError:
        raise TypeError("H must be convertible to a scipy.csc_matrix")
    unique_elements = np.unique(H.data)
    if len(unique_elements) > 1 or unique_elements[0] != 1:
        raise ValueError("Nonzero elements in the parity check matrix" \
                         " must be 1, not {}.".format(unique_elements))
    H = H.astype(np.uint8)
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

    column_weights = np.asarray(H.sum(axis=0))[0]
    unique_column_weights = np.unique(column_weights)
    if np.setdiff1d(unique_column_weights, np.array([1, 2])).size > 0:
        raise ValueError("Each column of H must have weight "
                         "1 or 2, not {}".format(unique_column_weights))
    H.eliminate_zeros()
    H.sort_indices()
    num_fault_ids = H.shape[1]

    if weights.shape[0] != num_fault_ids:
        raise ValueError("Weights array must have num_fault_ids elements")

    timelike_weights = 1.0 if timelike_weights is None else timelike_weights
    if isinstance(timelike_weights, (int, float, np.integer, np.floating)):
        timelike_weights = np.ones(H.shape[0], dtype=float) * timelike_weights
    elif isinstance(timelike_weights, (np.ndarray, list)):
        timelike_weights = np.array(timelike_weights, dtype=float)
        if timelike_weights.shape != (H.shape[0],):
            raise ValueError("timelike_weights should have the same number of elements as there are rows in H")
    else:
        raise ValueError("timelike_weights should be a float or a 1d numpy array")

    repetitions = 1 if repetitions is None else repetitions

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
        if p_meas.shape != (H.shape[0],):
            raise ValueError("measurement_error_probabilities should have dimensions {}"
                             " not {}".format((H.shape[0],), p_meas.shape))
    else:
        raise ValueError("measurement_error_probabilities should be a float or 1d numpy array")

    boundary = {H.shape[0] * repetitions} if 1 in unique_column_weights else set()
    self._matching_graph = MatchingGraph(H.shape[0] * repetitions + len(boundary))
    self._matching_graph.set_boundary(boundary=boundary)
    for t in range(repetitions):
        for i in range(len(H.indptr) - 1):
            s, e = H.indptr[i:i + 2]
            v1 = H.indices[s] + H.shape[0] * t
            v2 = H.indices[e - 1] + H.shape[0] * t if e - s == 2 else next(iter(boundary))
            self._matching_graph.add_edge(v1, v2, {i}, weights[i],
                                          error_probabilities[i], merge_strategy=merge_strategy)
    for t in range(repetitions - 1):
        for i in range(H.shape[0]):
            self._matching_graph.add_edge(i + t * H.shape[0], i + (t + 1) * H.shape[0],
                                          set(), timelike_weights[i], p_meas[i], merge_strategy=merge_strategy)


def load_from_detector_error_model(self, model: stim.DetectorErrorModel) -> None:
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
