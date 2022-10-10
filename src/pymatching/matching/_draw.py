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

import matplotlib.cbook
import numpy as np
import networkx as nx


def draw(self: 'pymatching.Matching') -> None:
    """Draw the matching graph using matplotlib
    Draws the matching graph as a matplotlib graph. Stabiliser nodes are
    filled grey and boundary nodes are filled white. The line thickness of each
    edge is determined from its weight (with min and max thicknesses of 0.2 pts
    and 2 pts respectively).
    Note that you may need to call `plt.figure()` before and `plt.show()` after calling
    this function.
    """
    # Ignore matplotlib deprecation warnings from networkx.draw_networkx
    warnings.filterwarnings("ignore", category=matplotlib.cbook.mplDeprecation)
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
