import numpy as np
import networkx as nx
from scipy.sparse import csc_matrix

from mwpm import PerfectMatching, Options


def stabiliser_graph(H):
    H = csc_matrix(H)
    if np.count_nonzero(H.indptr[1:]-H.indptr[0:-1]-2) != 0:
        raise ValueError("Parity check matrix does not have two "
                         "non-zero entries per column")
    G = nx.Graph()
    for i in range(len(H.indices)//2):
        G.add_edge(H.indices[2*i], H.indices[2*i+1],id=i)
    return G


def syndrome_graph(G, z):
    vs = np.nonzero(z)[0]
    S = nx.Graph()
    for i, v in enumerate(vs[0:-1]):
        length = nx.single_source_shortest_path_length(G, v)
        for other in vs[i+1:]:
            S.add_edge(v, other, weight=length[other])
    return S


def solve_matching(S, verbose=False):
    nodes = list(S)
    mapping = dict(zip(nodes,range(len(nodes))))
    inverse_mapping = dict(zip(range(len(nodes)), nodes))
    G = nx.relabel_nodes(S, mapping)
    pm = PerfectMatching(len(G), G.number_of_edges())
    for i, j, data in G.edges(data=True):
        pm.add_edge(i, j, data['weight'])
    pm.options.verbose=verbose
    pm.solve()
    matches = []
    for i in range(len(S)):
        j = pm.get_match(i)
        if i < j:
            matches.append((inverse_mapping[i],inverse_mapping[j]))
    return matches


def correction(G, matches):
    es = []
    for i, j in matches:
        path = nx.bidirectional_shortest_path(G, i, j)
        for k in range(len(path)-1):
            e = G.get_edge_data(path[k],path[k+1])['id']
            es.append(e)
    c = np.zeros(G.number_of_edges())
    for i in es:
        c[i]=1
    return c


def decode(H, z):
    G = stabiliser_graph(H)
    S = syndrome_graph(G, z)
    matches = solve_matching(S)
    return correction(G, matches)
