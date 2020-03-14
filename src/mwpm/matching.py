import numpy as np
import networkx as nx
from scipy.sparse import csc_matrix

from mwpm._cpp_mwpm import all_pairs_shortest_path, decode, StabiliserGraph


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
    def __init__(self, H):
        self.H = csc_matrix(H)
        check_two_checks_per_qubit(self.H)
        self.H.sort_indices()
        self.stabiliser_graph = StabiliserGraph(self.H.indices, self.H.shape[0], self.H.shape[1])
    
    def decode(self, z):
        if len(z.shape) == 1 and z.shape[0] == self.H.shape[0]:
            defects = z.nonzero()[0]
        elif len(z.shape) == 2 and z.shape[0] == self.H.shape[0]:
            times, checks = z.T.nonzero()
            defects = times*self.H.shape[0] + checks
        else:
            raise ValueError(f"The shape ({z.shape}) of the syndrome vector z is not valid.")
        if (len(defects) % 2) != 0:
            raise ValueError(f"There must be an even number of nonzero syndrome elements.")
        return decode(self.stabiliser_graph, defects)
