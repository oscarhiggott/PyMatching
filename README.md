# PyMatching 2

![Continuous Integration](https://github.com/oscarhiggott/PyMatching/workflows/ci/badge.svg)
[![codecov](https://codecov.io/gh/oscarhiggott/PyMatching/branch/master/graph/badge.svg)](https://codecov.io/gh/oscarhiggott/PyMatching)
[![docs](https://readthedocs.org/projects/pymatching/badge/?version=latest&style=plastic)](https://readthedocs.org/projects/pymatching/builds/)
[![PyPI version](https://badge.fury.io/py/PyMatching.svg)](https://badge.fury.io/py/PyMatching)
[![Unitary Fund](https://img.shields.io/badge/Supported%20By-UNITARY%20FUND-brightgreen.svg?style=for-the-badge)](http://unitary.fund)

PyMatching is a fast Python/C++ library for decoding quantum error correcting (QEC) codes using the Minimum Weight
Perfect Matching (MWPM) decoder.
Given the syndrome measurements from a quantum error correction circuit, the MWPM decoder finds the most probable set 
of errors, given the assumption that error mechanisms are _independent_, as well as _graphlike_ (each error causes 
either one or two detection events).
The MWPM decoder is the most popular decoder for decoding [surface codes](https://arxiv.org/abs/quant-ph/0110143), 
and can also be used to decode various other code families, including 
[subsystem codes](https://arxiv.org/abs/1207.1443), 
[honeycomb codes](https://quantum-journal.org/papers/q-2021-10-19-564/) and 
[2D hyperbolic codes](https://arxiv.org/abs/1506.04029).

Version 2 includes a new implementation of the blossom algorithm which is **100-1000x faster** than previous versions
of PyMatching.
PyMatching can be configured using arbitrary weighted graphs, with or without a boundary, and can be combined with 
Craig Gidney's [Stim](https://github.com/quantumlib/Stim) library to simulate and decode error correction circuits 
in the presence of circuit-level noise. The [sinter](https://pypi.org/project/sinter/) package combines Stim and 
PyMatching to perform fast, parallelised monte-carlo sampling of quantum error correction circuits.

Documentation for PyMatching can be found at: [pymatching.readthedocs.io](https://pymatching.readthedocs.io/en/stable/)

Our [paper](https://arxiv.org/abs/2303.15933) gives more background on the MWPM decoder and our implementation (sparse blossom) released in PyMatching v2.

To see how stim, sinter and pymatching can be used to estimate the threshold of an error correcting code with 
circuit-level noise, try out the [stim getting started notebook](https://github.com/quantumlib/Stim/blob/main/doc/getting_started.ipynb).

## The new >100x faster implementation for Version 2

Version 2 features a new implementation of the blossom algorithm, which I wrote with Craig Gidney.
Our new implementation, which we refer to as the _sparse blossom_ algorithm, can be seen as a generalisation of the 
blossom algorithm to handle the decoding problem relevant to QEC. 
We solve the problem of finding minimum-weight paths between detection events in a detector graph 
_directly_, which avoids the need to use costly all-to-all Dijkstra searches to find a MWPM in a derived 
graph using the original blossom algorithm.
The new version is also exact - unlike previous versions of PyMatching, no approximation is made.
See our [paper](https://arxiv.org/abs/2303.15933) for more details.

Our new implementation is **over 100x faster** than previous versions of PyMatching, and is 
**over 100,000x faster** than NetworkX (benchmarked with surface code circuits). At 0.1% circuit-noise, PyMatching can 
decode both X and Z basis measurements of surface code circuits up to distance 17 in under 1 microsecond per round 
of syndrome extraction on a single core. Furthermore, the runtime is roughly linear in the number 
of nodes in the graph.

The plot below compares the performance of PyMatching v2 with the previous 
version (v0.7) as well as with NetworkX for decoding surface code circuits with circuit-level depolarising noise. 
All decoders were run on a single core of an M1 processor, processing both the X and Z basis measurements.
The equations T=N^x in the legend (and plotted as dashed lines) are 
obtained from a fit to the same dataset for 
distance > 10, where N is the number of detectors (nodes) per round, and T is the decoding time per round.
See the [benchmarks](https://github.com/oscarhiggott/PyMatching/raw/master/benchmarks) folder in the repository 
for the data and stim circuits, as well as additional benchmarks.


![PyMatching new vs old vs NetworkX](https://github.com/oscarhiggott/PyMatching/raw/master/benchmarks/surface_codes/surface_code_rotated_memory_x_p_0.001_d_5_7_9_13_17_23_29_39_50_both_bases/pymatching_v0.7_vs_pymatching_v2_vs_networkx_timing_p=0.001_per_round_both_bases_decoded.png)


Sparse blossom is conceptually similar to the approach described in [this paper](https://arxiv.org/abs/1307.1740) 
by Austin Fowler, although our approach differs in many of the details (as explained in our [paper](https://arxiv.org/abs/2303.15933)).
There are even more similarities with the very nice independent work by Yue Wu, who recently released the 
[fusion-blossom](https://pypi.org/project/fusion-blossom/) library.
One of the differences with our approach is that fusion-blossom grows the exploratory regions of alternating trees 
in a similar way to how clusters are grown in Union-Find, whereas our approach instead progresses along a timeline, 
and uses a global priority queue to grow alternating trees.
Yue also has a paper coming soon, so stay tuned for that as well.

## Installation

The latest version of PyMatching can be downloaded and installed from [PyPI](https://pypi.org/project/PyMatching/) 
with the command:

```
pip install pymatching --upgrade
```


## Usage

PyMatching can load matching graphs from a check matrix, a `stim.DetectorErrorModel`, a `networkx.Graph`, a 
`rustworkx.PyGraph` or by adding edges individually with `pymatching.Matching.add_edge` and 
`pymatching.Matching.add_boundary_edge`.

### Decoding Stim circuits

PyMatching can be combined with [Stim](https://github.com/quantumlib/Stim). Generally, the easiest and fastest way to 
do this is using [sinter](https://pypi.org/project/stim/) (use v1.10.0 or later), which uses PyMatching and Stim to run 
parallelised monte carlo simulations of quantum error correction circuits.
However, in this section we will use Stim and PyMatching directly, to demonstrate how their Python APIs can be used.
To install stim, run `pip install stim --upgrade`.

First, we generate a stim circuit. Here, we use a surface code circuit included with stim:

```python
import numpy as np
import stim
import pymatching
circuit = stim.Circuit.generated("surface_code:rotated_memory_x", 
                                 distance=5, 
                                 rounds=5, 
                                 after_clifford_depolarization=0.005)
```

Next, we use stim to generate a `stim.DetectorErrorModel` (DEM), which is effectively a 
[Tanner graph](https://en.wikipedia.org/wiki/Tanner_graph) describing the circuit-level noise model.
By setting `decompose_errors=True`, stim decomposes all error mechanisms into _edge-like_ error 
mechanisms (which cause either one or two detection events).
This ensures that our DEM is graphlike, and can be loaded by pymatching:

```python
model = circuit.detector_error_model(decompose_errors=True)
matching = pymatching.Matching.from_detector_error_model(model)
```

Next, we will sample 1000 shots from the circuit. Each shot (a row of `shots`) contains the full syndrome (detector 
measurements), as well as the logical observable measurements, from simulating the noisy circuit:

```python
sampler = circuit.compile_detector_sampler()
syndrome, actual_observables = sampler.sample(shots=1000, separate_observables=True)
```

Now we can decode! We compare PyMatching's predictions of the logical observables with the actual observables sampled 
with stim, in order to count the number of mistakes and estimate the logical error rate:

```python
num_errors = 0
for i in range(syndrome.shape[0]):
    predicted_observables = matching.decode(syndrome[i, :])
    num_errors += not np.array_equal(actual_observables[i, :], predicted_observables)

print(num_errors)  # prints 8
```

As of PyMatching v2.1.0, you can use `matching.decode_batch` to decode a batch of shots instead.
Since `matching.decode_batch` iterates over the shots in C++, it's faster than iterating over calls 
to `matching.decode` in Python. The following cell is therefore a faster 
equivalent to the cell above:

```python
predicted_observables = matching.decode_batch(syndrome)
num_errors = np.sum(np.any(predicted_observables != actual_observables, axis=1))

print(num_errors)  # prints 8
```

### Loading from a parity check matrix

We can also load a `pymatching.Matching` object from a binary
[parity check matrix](https://en.wikipedia.org/wiki/Parity-check_matrix), another representation of a Tanner graph.
Each row in the parity check matrix `H` corresponds to a parity check, and each column corresponds to an 
error mechanism.
The element `H[i,j]` of `H` is 1 if parity check `i` is flipped by error mechanism `j`, and 0 otherwise.
To be used by PyMatching, the error mechanisms in `H` must be _graphlike_.
This means that each column must contain either one or two 1s (if a column has a single 1, it represents a half-edge 
connected to the boundary).

We can give each edge in the graph a weight, by providing PyMatching with a `weights` numpy array.
Element `weights[j]` of the `weights` array sets the edge weight for the edge corresponding to column `j` of `H`.
If the error mechanisms are treated as independent, then we typically want to set the weight of edge `j` to 
the log-likelihood ratio `log((1-p_j)/p_j)`, where `p_j` is the error probability associated with edge `j`.
With this setting, PyMatching will find the most probable set of error mechanisms, given the syndrome.

With PyMatching configured using `H` and `weights`, decoding a binary syndrome vector `syndrome` (a numpy array 
of length `H.shape[0]`) corresponds to finding a set of errors defined in a binary `predictions` vector 
satisfying `H@predictions % 2 == syndrome` while minimising the total solution weight `predictions@weights`.

In quantum error correction, rather than predicting which exact set of error mechanisms occurred, we typically want to 
predict the outcome of _logical observable_ measurements, which are the parities of error mechanisms.
These can be represented by a binary matrix `observables`. Similar to the check matrix, `observables[i,j]` is 1 if 
logical observable `i` is flipped by error mechanism `j`.
For example, suppose our syndrome `syndrome`, was the result of a set of errors `noise` (a binary array of 
length `H.shape[1]`), such that `syndrome = H@noise % 2`.
Our decoding is successful if `observables@noise % 2 == observables@predictions % 2`.

Putting this together, we can decode a distance 5 repetition code as follows:

```python
import numpy as np
from scipy.sparse import csc_matrix
import pymatching
H = csc_matrix([[1, 1, 0, 0, 0],
                 [0, 1, 1, 0, 0],
                 [0, 0, 1, 1, 0],
                 [0, 0, 0, 1, 1]])
weights = np.array([4, 3, 2, 3, 4])   # Set arbitrary weights for illustration
matching = pymatching.Matching(H, weights=weights)
prediction = matching.decode(np.array([0, 1, 0, 1]))
print(prediction)  # prints: [0 0 1 1 0]
# Optionally, we can return the weight as well:
prediction, solution_weight = matching.decode(np.array([0, 1, 0, 1]), return_weight=True)
print(prediction)  # prints: [0 0 1 1 0]
print(solution_weight)  # prints: 5.0
```

And in order to estimate the logical error rate for a physical error rate of 10%, we can sample 
as follows:

```python
import numpy as np
from scipy.sparse import csc_matrix
import pymatching
H = csc_matrix([[1, 1, 0, 0, 0],
                [0, 1, 1, 0, 0],
                [0, 0, 1, 1, 0],
                [0, 0, 0, 1, 1]])
observables = csc_matrix([[1, 0, 0, 0, 0]])
error_probability = 0.1
weights = np.ones(H.shape[1]) * np.log((1-error_probability)/error_probability)
matching = pymatching.Matching.from_check_matrix(H, weights=weights)
num_shots = 1000
num_errors = 0
for i in range(num_shots):
    noise = (np.random.random(H.shape[1]) < error_probability).astype(np.uint8)
    syndrome = H@noise % 2
    prediction = matching.decode(syndrome)
    predicted_observables = observables@prediction % 2
    actual_observables = observables@noise % 2
    num_errors += not np.array_equal(predicted_observables, actual_observables)
print(num_errors)  # prints 4
```

Note that we can also ask PyMatching to predict the logical observables directly, by supplying them 
to the `faults_matrix` argument when constructing the `pymatching.Matching` object. This allows the decoder to make 
some additional optimisations, that speed up the decoding procedure a bit. The following example uses this approach, 
and is equivalent to the example above:

```python
import numpy as np
from scipy.sparse import csc_matrix
import pymatching

H = csc_matrix([[1, 1, 0, 0, 0],
                [0, 1, 1, 0, 0],
                [0, 0, 1, 1, 0],
                [0, 0, 0, 1, 1]])
observables = csc_matrix([[1, 0, 0, 0, 0]])
error_probability = 0.1
weights = np.ones(H.shape[1]) * np.log((1-error_probability)/error_probability)
matching = pymatching.Matching.from_check_matrix(H, weights=weights, faults_matrix=observables)
num_shots = 1000
num_errors = 0
for i in range(num_shots):
    noise = (np.random.random(H.shape[1]) < error_probability).astype(np.uint8)
    syndrome = H@noise % 2
    predicted_observables = matching.decode(syndrome)
    actual_observables = observables@noise % 2
    num_errors += not np.array_equal(predicted_observables, actual_observables)

print(num_errors)  # prints 6
```

We'll make one more optimisation, which is to use `matching.decode_batch` to decode the batch of shots, rather than
iterating over calls to `matching.decode` in Python:

```python
import numpy as np
from scipy.sparse import csc_matrix
import pymatching

H = csc_matrix([[1, 1, 0, 0, 0],
                [0, 1, 1, 0, 0],
                [0, 0, 1, 1, 0],
                [0, 0, 0, 1, 1]])
observables = csc_matrix([[1, 0, 0, 0, 0]])
error_probability = 0.1
num_shots = 1000
weights = np.ones(H.shape[1]) * np.log((1-error_probability)/error_probability)
matching = pymatching.Matching.from_check_matrix(H, weights=weights, faults_matrix=observables)
noise = (np.random.random((num_shots, H.shape[1])) < error_probability).astype(np.uint8)
shots = (noise @ H.T) % 2
actual_observables = (noise @ observables.T) % 2
predicted_observables = matching.decode_batch(shots)
num_errors = np.sum(np.any(predicted_observables != actual_observables, axis=1))
print(num_errors)  # prints 6

```



Instead of using a check matrix, the Matching object can also be constructed using
the [`Matching.add_edge`](https://pymatching.readthedocs.io/en/stable/api.html#pymatching.matching.Matching.add_edge)
and 
[`Matching.add_boundary_edge`](https://pymatching.readthedocs.io/en/stable/api.html#pymatching.matching.Matching.add_boundary_edge) 
methods, or by loading from a NetworkX or rustworkx graph. 

For more details on how to use PyMatching,
see [the documentation](https://pymatching.readthedocs.io).

## Attribution

When using PyMatching please cite our [paper](https://arxiv.org/abs/2303.15933) on the sparse blossom algorithm (implemented in version 2):

```
@article{Higgott2025sparseblossom,
  doi = {10.22331/q-2025-01-20-1600},
  url = {https://doi.org/10.22331/q-2025-01-20-1600},
  title = {Sparse {B}lossom: correcting a million errors per core second with minimum-weight matching},
  author = {Higgott, Oscar and Gidney, Craig},
  journal = {{Quantum}},
  issn = {2521-327X},
  publisher = {{Verein zur F{\"{o}}rderung des Open Access Publizierens in den Quantenwissenschaften}},
  volume = {9},
  pages = {1600},
  month = jan,
  year = {2025}
}
```

Note: the previous PyMatching [paper](https://arxiv.org/abs/2105.13082) descibes the implementation in version 0.7 and 
earlier of PyMatching (not v2).

## Acknowledgements

We are grateful to the Google Quantum AI team for supporting the development of PyMatching v2. Earlier versions of 
PyMatching were supported by Unitary Fund and EPSRC.
