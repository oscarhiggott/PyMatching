import networkx as nx
import numpy as np

from pymatching import Matching


def test_add_noise():
    p = 0.1
    N = 1000
    std = (p*(1-p)/N)**0.5
    g = nx.Graph()
    for i in range(N):
        g.add_edge(i, i+1, qubit_id=i, weight=-np.log(p), error_probability=p)
    m = Matching(g)
    for i in range(5):
        noise, syndrome = m.add_noise()
        assert (sum(syndrome) % 2) == 0
        assert (p-5*std) * N < sum(noise) < (p+5*std) * N
        for i in range(1, N-1):
            assert syndrome[i] == (noise[i-1] + noise[i]) % 2


def test_add_noise_with_boundary():
    g = nx.Graph()
    for i in range(11):
        g.add_edge(i, i+1, qubit_id=i, error_probability=(i+1) % 2)
    for i in range(5, 12):
        g.nodes()[i]['is_boundary'] = True
    m = Matching(g)
    noise, syndrome = m.add_noise()
    assert sum(syndrome) == 6
    assert np.array_equal(noise, (np.arange(11)+1) % 2)
    assert m.boundary == list(range(5, 12))
    assert np.array_equal(syndrome, np.array([1,1,1,1,1,1,0,0,0,0,0,0]))
