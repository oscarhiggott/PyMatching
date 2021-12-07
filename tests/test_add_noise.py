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

import networkx as nx
import numpy as np
from scipy.sparse import csr_matrix

from pymatching import Matching


def test_add_noise():
    p = 0.1
    N = 1000
    std = (p*(1-p)/N)**0.5
    g = nx.Graph()
    for i in range(N):
        g.add_edge(i, i+1, fault_ids=i, weight=-np.log(p), error_probability=p)
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
        g.add_edge(i, i+1, fault_ids=i, error_probability=(i+1) % 2)
    for i in range(5, 12):
        g.nodes()[i]['is_boundary'] = True
    m = Matching(g)
    noise, syndrome = m.add_noise()
    assert sum(syndrome) == 5
    assert np.array_equal(noise, (np.arange(11)+1) % 2)
    assert m.boundary == set(range(5, 12))
    assert np.array_equal(
        syndrome,
        np.array([1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0])
    )


def test_add_noise_without_error_probabilities_returns_none():
    m = Matching(csr_matrix(np.array([[1,1,0],[0,1,1]])))
    assert m.add_noise() is None
    m = Matching(csr_matrix(np.array([[1,1,0],[0,1,1]])),
             error_probabilities=np.array([0.5,0.7,-0.1]))
    assert m.add_noise() is None
