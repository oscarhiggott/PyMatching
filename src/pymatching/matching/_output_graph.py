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

import warnings
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    import pymatching

import networkx as nx
import retworkx as rx


def to_networkx(self: 'pymatching.Matching') -> nx.Graph:
    """Convert to NetworkX graph
    Returns a NetworkX graph corresponding to the matching graph. Each edge
    has attributes `fault_ids`, `weight` and `error_probability` and each node has
    the attribute `is_boundary`.
    Returns
    -------
    NetworkX.Graph
        NetworkX Graph corresponding to the matching graph
    """
    G = nx.Graph()
    G.add_edges_from(self.edges())
    boundary = self.boundary
    for i in G.nodes:
        is_boundary = i in boundary
        G.nodes[i]['is_boundary'] = is_boundary
    return G


def to_retworkx(self: 'pymatching.Matching') -> rx.PyGraph:
    """Convert to retworkx graph
    Returns a retworkx graph object corresponding to the matching graph. Each edge
    payload is a ``dict`` with keys `fault_ids`, `weight` and `error_probability` and
    each node has a ``dict`` payload with the key ``is_boundary`` and the value is
    a boolean.
    Returns
    -------
    retworkx.PyGraph
        retworkx graph corresponding to the matching graph
    """
    G = rx.PyGraph(multigraph=False)
    G.add_nodes_from([{} for _ in range(self.num_nodes)])
    G.extend_from_weighted_edge_list(self.edges())
    boundary = self.boundary
    for i in G.node_indices():
        is_boundary = i in boundary
        G[i]['is_boundary'] = is_boundary
    return G
