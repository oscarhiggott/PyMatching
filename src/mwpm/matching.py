import numpy as np
from scipy.sparse import csc_matrix

from mwpm._cpp_mwpm import _stabiliser_graph, all_pairs_shortest_path, _decode, StabiliserGraph


class MWPM:
    def __init__(self, H):
        self.H = csc_matrix(H)
        if np.count_nonzero(self.H.indptr[1:]-self.H.indptr[0:-1]-2) != 0:
            raise ValueError("Parity check matrix does not have two "
                             "non-zero entries per column")
        self.stabiliser_graph = StabiliserGraph(self.H.indices, self.H.shape[0], self.H.shape[1])
    
    def decode(self, z):
        defects = np.nonzero(z)[0]
        print(f"Running decode")
        return _decode(self.stabiliser_graph, defects)
