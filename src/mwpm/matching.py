import numpy as np
from scipy.sparse import csc_matrix

from mwpm._cpp_mwpm import all_pairs_shortest_path, decode, StabiliserGraph


class MWPM:
    def __init__(self, H):
        self.H = csc_matrix(H)
        if np.count_nonzero(self.H.indptr[1:]-self.H.indptr[0:-1]-2) != 0:
            raise ValueError("Parity check matrix does not have two "
                             "non-zero entries per column")
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
