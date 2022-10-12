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

from typing import Set, TYPE_CHECKING
if TYPE_CHECKING:
    import pymatching


def set_boundary_nodes(self: 'pymatching.Matching', nodes: Set[int]) -> None:
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


def set_min_num_fault_ids(self: 'pymatching.Matching', min_num_fault_ids: int) -> None:
    """
    Set the minimum number of fault ids in the matching graph.

    Let `max_id` be the maximum fault id assigned to any of the edges in the graph. Then setting
    `min_num_fault_ids` will ensure that `Matching.num_fault_ids=max(min_num_fault_ids, max_id)`.
    Note that `Matching.num_fault_ids` sets the length of the correction array output by `Matching.decode`.
    Parameters
    ----------
    min_num_fault_ids: int
        The required minimum number of fault ids in the matching graph

    """
    self._matching_graph.set_min_num_observables(min_num_fault_ids)


@property
def num_fault_ids(self: 'pymatching.Matching') -> int:
    """
    The number of fault IDs defined in the matching graph
    Returns
    -------
    int
        Number of fault IDs
    """
    return self._matching_graph.get_num_observables()


@property
def boundary(self: 'pymatching.Matching') -> Set[int]:
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
def num_nodes(self: 'pymatching.Matching') -> int:
    """
    The number of nodes in the matching graph
    Returns
    -------
    int
        The number of nodes
    """
    return self._matching_graph.get_num_nodes()


@property
def num_edges(self: 'pymatching.Matching') -> int:
    """
    The number of edges in the matching graph
    Returns
    -------
    int
        The number of edges
    """
    return self._matching_graph.get_num_edges()


@property
def num_detectors(self: 'pymatching.Matching') -> int:
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
