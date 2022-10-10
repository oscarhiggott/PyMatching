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

from pymatching._cpp_pymatching import MatchingGraph as _MatchingGraph


class Matching:
    """A class for constructing matching graphs and decoding using the minimum-weight perfect matching decoder.
    The matching graph can be constructed using the `Matching.add_edge` and `Matching.add_boundary_edge`
    methods. Alternatively, it can be loaded from a parity check matrix (a `scipy.sparse` matrix or `numpy.ndarray`
    with one or two non-zero elements in each column), a NetworkX or retworkx graph, or from
    a `stim.DetectorErrorModel`.
    """
    from pymatching.matching._edges import (add_edge, add_boundary_edge, edges, has_edge, has_boundary_edge,
                                            get_edge_data, get_boundary_edge_data)
    from pymatching.matching._add_noise import add_noise
    from pymatching.matching._decode import decode
    from pymatching.matching._draw import draw
    from pymatching.matching._load_from_check_matrix import load_from_check_matrix
    from pymatching.matching._load_from_detector_error_model import load_from_detector_error_model
    from pymatching.matching._load_from_networkx import load_from_networkx
    from pymatching.matching._load_from_retworkx import load_from_retworkx
    from pymatching.matching._output_graph import (to_networkx, to_retworkx)
    from pymatching.matching._properties import (set_boundary_nodes, num_fault_ids, boundary, num_nodes,
                              num_edges, num_detectors)

    def __init__(self,
                 graph: Union[scipy.sparse.spmatrix, np.ndarray, rx.PyGraph, nx.Graph, List[
                     List[int]], stim.DetectorErrorModel] = None,
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
            matrix or numpy.ndarray), a NetworkX or retworkx graph, or a Stim DetectorErrorModel.
            Each edge in the NetworkX or retworkx graph can have optional
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
            the measurements are noisy. This option is only used if `H` is
            provided as a check matrix, not a NetworkX graph. By default None
        timelike_weights : float, optional
            If `H` is given as a scipy or numpy array and `repetitions>1`,
            `timelike_weights` gives the weight of timelike edges.
            If a float is given, all timelike edges weights are set to
            the same value. If a numpy array of size `(H.shape[0],)` is given, the
            edge weight for each vertical timelike edge associated with the `i`th check (row)
            of `H` is set to `timelike_weights[i]`. By default None, in which case all
            timelike weights are set to 1.0
        measurement_error_probabilities : float, optional
            If `H` is given as a scipy or numpy array and `repetitions>1`,
            gives the probability of a measurement error to be used for
            the add_noise method. If a float is given, all measurement
            errors are set to the same value. If a numpy array of size `(H.shape[0],)` is given,
            the error probability for each vertical timelike edge associated with the `i`th check
            (row) of `H` is set to `measurement_error_probabilities[i]`. By default None
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
        self._matching_graph = _MatchingGraph()
        if graph is None:
            graph = kwargs.get("H")
            if graph is None:
                return
            del kwargs["H"]
        if isinstance(graph, nx.Graph):
            self.load_from_networkx(graph)
        elif isinstance(graph, rx.PyGraph):
            self.load_from_retworkx(graph)
        elif isinstance(graph, stim.DetectorErrorModel):
            self.load_from_detector_error_model(graph)
        else:
            try:
                graph = csc_matrix(graph)
            except TypeError:
                raise TypeError("The type of the input graph is not recognised. `graph` must be "
                                "a scipy.sparse or numpy matrix, networkx or retworkx graph, or "
                                "stim.DetectorErrorModel.")
            self.load_from_check_matrix(graph, weights, error_probabilities,
                                        repetitions, timelike_weights, measurement_error_probabilities,
                                        **kwargs)

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
