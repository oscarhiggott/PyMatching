import networkx as nx
import numpy as np

from mwpm import MWPM


def test_add_noise():
    p = 0.1
    N = 1000
    std = (p*(1-p)/N)**0.5
    g = nx.Graph()
    for i in range(N):
        g.add_edge(i, i+1, qubit_id=i, weight=-np.log(p), error_probability=p)
    m = MWPM(g)
    for i in range(5):
        noise, syndrome = m.add_noise()
        assert (sum(syndrome) % 2) == 0
        assert (p-5*std) * N < sum(noise) < (p+5*std) * N
        for i in range(1, N-1):
            assert syndrome[i] == (noise[i-1] + noise[i]) % 2
