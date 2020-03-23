import numpy as np
import networkx as nx
from scipy.sparse import csc_matrix, spmatrix
import numpy as np

from mwpm._cpp_mwpm import (all_pairs_shortest_path, 
                            decode, 
                            UnweightedStabiliserGraph,
                            WeightedStabiliserGraph)


def check_two_checks_per_qubit(H):
    if np.count_nonzero(H.indptr[1:]-H.indptr[0:-1]-2) != 0:
        raise ValueError("Parity check matrix does not have two "
                            "non-zero entries per column")


def syndrome_graph_from_check_matrix(H):
    H = csc_matrix(H)
    check_two_checks_per_qubit(H)
    H.sort_indices()
    G = nx.Graph()
    for i in range(len(H.indices)//2):
        G.add_edge(H.indices[2*i], H.indices[2*i+1], id=i)
    return G


class MWPM:
    def __init__(self, H, weights=None):
        if not isinstance(H, nx.Graph):
            H = csc_matrix(H)
            check_two_checks_per_qubit(H)
            H.sort_indices()
            self.num_stabilisers = H.shape[0]
            num_qubits = H.shape[1]
            if weights is None:
                self.stabiliser_graph = UnweightedStabiliserGraph(
                    H.indices
                )
            else:
                weights = np.asarray(weights)
                if weights.shape[0] != self.num_qubits:
                    raise ValueError("Weights array must have num_qubits elements")
                if np.any(weights >= 0.):
                    raise ValueError("All weights must be non-negative.")
                self.stabiliser_graph = WeightedStabiliserGraph(
                    H.indices,
                    weights
                )
        else:
            self.num_stabilisers = H.number_of_nodes()
            g = WeightedStabiliserGraph(self.num_stabilisers)
            for (u, v, attr) in H.edges(data=True):
                if ('qubit_id' not in attr) or ('weight' not in attr):
                    raise ValueError("Graph edges must have 'qubit_id' and 'weight' attributes."
                                f" Instead attr={attr}")
                g.add_edge(u, v, attr['qubit_id'], attr['weight'])
            g.compute_all_pairs_shortest_paths()
            self.stabiliser_graph = g
    
    @property
    def num_qubits(self):
        return self.stabiliser_graph.get_num_qubits()
    
    def decode(self, z):
        if len(z.shape) == 1 and z.shape[0] == self.num_stabilisers:
            defects = z.nonzero()[0]
        elif len(z.shape) == 2 and z.shape[0] == self.num_stabilisers:
            times, checks = z.T.nonzero()
            defects = times*self.num_stabilisers + checks
        else:
            raise ValueError(f"The shape ({z.shape}) of the syndrome vector z is not valid.")
        if (len(defects) % 2) != 0:
            raise ValueError(f"There must be an even number of nonzero syndrome elements.")
        return decode(self.stabiliser_graph, defects)
