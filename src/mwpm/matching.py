import numpy as np
from scipy.sparse import csc_matrix

from mwpm._cpp_mwpm import _stabiliser_graph, all_pairs_shortest_path, _decode


class MWPM:
    def __init__(self, H):
        self.H = csc_matrix(H)
        if np.count_nonzero(self.H.indptr[1:]-self.H.indptr[0:-1]-2) != 0:
            raise ValueError("Parity check matrix does not have two "
                             "non-zero entries per column")
        self.graph_data = _stabiliser_graph(self.H.indices, self.H.shape[0])
        self.shortest_paths = all_pairs_shortest_path(self.graph_data.g)
    
    def decode(self, z):
        defects = np.nonzero(z)[0]
        return _decode(self.shortest_paths, defects, 
                       self.graph_data.qubit, self.H.shape[1])
