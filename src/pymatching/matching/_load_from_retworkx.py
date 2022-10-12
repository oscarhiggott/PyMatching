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

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    import pymatching

import numpy as np
import retworkx as rx

from pymatching._cpp_pymatching import MatchingGraph


def load_from_retworkx(self: 'pymatching.Matching', graph: rx.PyGraph, *, min_num_fault_ids: int = None) -> None:
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
    min_num_fault_ids: int
        Sets the minimum number of fault ids in the matching graph. Let `max_id` be the maximum fault id assigned to
        any of the edges in the graph. Then setting this argument will ensure that
        `Matching.num_fault_ids=max(min_num_fault_ids, max_id)`. Note that `Matching.num_fault_ids` sets the length
        of the correction array output by `Matching.decode`.
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
    num_fault_ids = 0 if min_num_fault_ids is None else min_num_fault_ids
    g = MatchingGraph(num_nodes, num_fault_ids)
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
