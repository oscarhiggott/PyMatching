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

from typing import Union, Set

import numpy as np


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
        The weight of the edge, which must be non-negative, by default 1.0
    error_probability: float, optional
        The probability that the edge is flipped. This is used by the `add_noise()` method
        to sample from the distribution defined by the matching graph (in which each edge
        is flipped independently with the corresponding `error_probability`). By default None
    merge_strategy: str, optional
        Which strategy to use if the edge (`node1`, `node2`) is already in the graph. The available options
        are "disallow", "independent", "smallest-weight", "first-only" and "last-only". "disallow" raises a
        `ValueError` if the edge (`node1`, `node2`) is already present. The "independent" strategy assumes that
        the existing edge (`node1`, `node2`) and the edge being added represent independent error mechanisms, and
        they are merged into a new edge with updated weights and error_probabilities accordingly (it is assumed
        that each weight represents the log-likelihood ratio log((1-p)/p) where p is the `error_probability` and
        where the natural logarithm is used. The fault_ids associated with the existing edge are kept only, since
        where the natural logarithm is used. The fault_ids associated with the existing edge are kept only, since
        the code has distance 2 if parallel edges have different fault_ids anyway). The "smallest-weight" strategy
        keeps only the new edge if it has a smaller weight than the existing edge, otherwise the graph is left
        unchanged. The "first-only" strategy keeps only the existing edge, and ignores the edge being added.
        The "last-only" strategy always keeps the edge being added, replacing the existing edge.
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
        The weight of the edge, which must be non-negative, by default 1.0
    error_probability: float, optional
        The probability that the edge is flipped. This is used by the `add_noise()` method
        to sample from the distribution defined by the matching graph (in which each edge
        is flipped independently with the corresponding `error_probability`). By default None
    merge_strategy: str, optional
        Which strategy to use if the edge (`node1`, `node2`) is already in the graph. The available options
        are "disallow", "independent", "smallest-weight", "first-only" and "last-only". "disallow" raises a
        `ValueError` if the edge (`node1`, `node2`) is already present. The "independent" strategy assumes that
        the existing edge (`node1`, `node2`) and the edge being added represent independent error mechanisms, and
        they are merged into a new edge with updated weights and error_probabilities accordingly (it is assumed
        that each weight represents the log-likelihood ratio log((1-p)/p) where p is the `error_probability` and
        where the natural logarithm is used. The fault_ids associated with the existing edge are kept only, since
        where the natural logarithm is used. The fault_ids associated with the existing edge are kept only, since
        the code has distance 2 if parallel edges have different fault_ids anyway). The "smallest-weight" strategy
        keeps only the new edge if it has a smaller weight than the existing edge, otherwise the graph is left
        unchanged. The "first-only" strategy keeps only the existing edge, and ignores the edge being added.
        The "last-only" strategy always keeps the edge being added, replacing the existing edge.
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
