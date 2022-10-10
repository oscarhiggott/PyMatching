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
import networkx as nx

from pymatching._cpp_pymatching import MatchingGraph


def load_from_networkx(self: 'pymatching.Matching', graph: nx.Graph) -> None:
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
