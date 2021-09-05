# PyMatching


![Continuous Integration](https://github.com/oscarhiggott/PyMatching/workflows/Continuous%20Integration/badge.svg)
[![codecov](https://codecov.io/gh/oscarhiggott/PyMatching/branch/master/graph/badge.svg)](https://codecov.io/gh/oscarhiggott/PyMatching)
[![docs](https://readthedocs.org/projects/pymatching/badge/?version=latest&style=plastic)](https://readthedocs.org/projects/pymatching/builds/)
[![PyPI version](https://badge.fury.io/py/PyMatching.svg)](https://badge.fury.io/py/PyMatching)
[![Unitary Fund](https://img.shields.io/badge/Supported%20By-UNITARY%20FUND-brightgreen.svg?style=for-the-badge)](http://unitary.fund)

PyMatching is a fast Python/C++ library for decoding quantum error correcting codes (QECC) using the Minimum Weight Perfect Matching (MWPM) decoder. PyMatching can decode codes for which each error generates a pair of syndrome defects (or only a single defect at a boundary). Codes that satisfy these properties include two-dimensional topological codes such as the [toric code](https://en.wikipedia.org/wiki/Toric_code), the [surface code](https://arxiv.org/abs/quant-ph/0110143) and [2D hyperbolic codes](https://arxiv.org/abs/1506.04029), amongst others. PyMatching can also be used as a subroutine to decode other codes, such as the 3D toric code and the [color code](https://arxiv.org/abs/1905.07393). PyMatching can handle boundaries, measurement errors and weighted edges in the matching graph. Since the core algorithms are written in C++, PyMatching is much faster than a pure Python NetworkX implementation.

Documentation for PyMatching can be found at: [pymatching.readthedocs.io](https://pymatching.readthedocs.io/en/stable/)

## Installation

PyMatching can be downloaded and installed from [PyPI](https://pypi.org/project/PyMatching/) with the command:
```
pip install pymatching
```
This is the recommended way to install PyMatching since pip will fetch the pre-compiled binaries, rather than building the C++ extension from source on your machine. Note that PyMatching requires Python 3.

If instead you would like to install PyMatching from source, clone the repository (using the `--recursive` flag to include the lib/pybind11 submodule) and then use `pip` to install:
```
git clone --recursive https://github.com/oscarhiggott/PyMatching.git
pip install -e ./PyMatching
```
The installation may take a few minutes since the C++ extension has to be compiled. If you'd also like to run the tests, first install [pytest](https://docs.pytest.org/en/stable/), and then run:
```
pytest ./PyMatching/tests ./PyMatching/src
```

## Usage

In order to decode a parity check matrix `H` (a `scipy.sparse` matrix) with syndrome vector `z` (a bitstring which is a numpy array of dtype int), first construct the `Matching` object after importing it:
```
from pymatching import Matching
m = Matching(H)
```
Now to decode, simply run:
```
c = m.decode(z)
```
which outputs a bitstring `c`, which is a numpy array of ints corresponding to the minimum-weight correction. Note that the `m` by `n` parity check matrix `H` should correspond to the Z (or X) stabilisers of a CSS code with `n` qubits, `m` Z (or X) stabilisers, and with either one or two non-zero entries per column.

To decode instead in the presence of measurement errors, each stabiliser measurement is repeated `L` times, and decoding then takes place over a 3D matching graph (see Section IV B of [this paper](https://arxiv.org/abs/quant-ph/0110143)), which can be constructed directly from the check matrix `H` using:
```
m = Matching(H, repetitions=L)
```
and then decoded from an `m` by `L` numpy array syndrome `z` using:
```
c = m.decode(z)
```

The Matching object can also be constructed from a NetworkX graph instead of a check matrix, and can handle weighted edges. For full details see [the documentation](https://pymatching.readthedocs.io/en/stable/usage.html).

## Performance

While all the functionality of PyMatching is available via the Python bindings, the core algorithms and data structures are implemented in C++, with the help of the [LEMON](https://lemon.cs.elte.hu/trac/lemon) and [Boost Graph](https://www.boost.org/doc/libs/1_74_0/libs/graph/doc/index.html) libraries. PyMatching also uses a local variant of the MWPM decoder (explained in the Appendix of [this paper](https://arxiv.org/abs/2010.09626)) that has a runtime that is approximately linear, rather than quadratic, in the number of nodes. As a result, PyMatching is orders of magnitude faster than a standard pure Python [NetworkX](https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.matching.max_weight_matching.html) implementation, as shown here for decoding the toric code under an independent noise model with p=0.05 and noiseless syndrome measurements:

<img src="https://raw.githubusercontent.com/oscarhiggott/PyMatching/master/docs/_static/pymatching_vs_networkx.png" width="400">

## Exact vs. local matching

PyMatching includes both the standard "exact" minimum-weight perfect matching decoder, as well as a close approximation of it, called _local matching_, which is much faster.
Local matching allows each node corresponding to a syndrome defect (-1 measurement) to be matched to one of the `num_neighbours` defects that are closest to it in the matching graph.
By default, PyMatching uses local matching with `num_neighbours=30`, but a different choice of `num_neighbours` can be set when decoding, e.g.:
```
c = m.decode(z, num_neighbours=40)
```
Note that by setting `num_neighbours=sum(z)`, local matching corresponds to exact matching.

Rather than setting `num_neighbours=sum(z)`, an alternative option for using exact matching is provided by setting `num_neighbours=None`. If this option is chosen, the shortest paths between all pairs of nodes in the matching graph are pre-computed and cached the first time `m.decode` is called, and then reused for later uses of `m.decode`. This differs from local matching, where shortest paths are computed on the fly.
As a result, setting `num_neighbours=None` is more memory intensive than local matching, with the required memory scaling quadratically with the number of nodes in the matching graph, however for exact matching it is faster than setting `num_neighbours=sum(z)`.

For typical decoding problems, local matching is an extremely close approximation of exact matching even for small `num_neighbours`. The following graph shows the threshold of local matching for the toric code with noisy syndrome measurements (a 3D matching graph), as a function of `num_neighbours`. For `num_neighbours>=16`, the local matching threshold is consistent with the 2.92\% threshold found with exact matching:

<img src="https://raw.githubusercontent.com/oscarhiggott/PyMatching/master/docs/_static/toric_noisy_syndromes_threshold_vs_num_neighbours.png" width=300>

The runtime of local matching scales linearly with `num_neighbours`, as shown by the following graph, generated using an L=20 toric code:

<img src="https://raw.githubusercontent.com/oscarhiggott/PyMatching/master/docs/_static/toric_timing_analysis_num_neighbours_vs_t.png" width=300>

A more detailed description and analysis of local matching can be found in the PyMatching [paper](https://arxiv.org/abs/2105.13082).

Note that PyMatching used `num_neighbours=20` as a default for v0.3.1 and earlier.

## Attribution

When using PyMatching for research, please cite:
```
@article{higgott2021pymatching,
  title={{PyMatching}: A Python package for decoding quantum codes with minimum-weight perfect matching},
  author={Higgott, Oscar},
  journal={arXiv preprint arXiv:2105.13082},
  year={2021}
}
```

Please also consider citing the [LEMON](https://dl.acm.org/doi/10.1016/j.entcs.2011.06.003) and Boost Graph libraries.
