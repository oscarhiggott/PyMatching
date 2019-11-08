# MWPM

A library for decoding quantum error correcting codes (QECC) using the Minimum Weight Perfect Matching decoder.

## Installation

Clone the repository (using the `--recursive` flag to include the lib/pybind11 submodule). Make sure you have cmake installed, and then use `pip` to install:
```
git clone --recursive https://github.com/oscarhiggott/MWPM.git
pip install -e ./MWPM
```

## Usage

In order to decode a parity check matrix `H` (a scipy.sparse matrix with) with syndrome vector `z` (a bitstring which is a numpy array of dtype int), first construct the `MWPM` object after importing it:
```
from mwpm import MWPM
m = MWPM(H)
```
This precomputes the all-pairs shortest paths in the stabiliser graph corresponding to the parity check matrix H. Now to decode, simply run:
```
c = m.decode(z)
```
which outputs a bitstring `c`, which is a numpy array of dtype int. Note that the `mxn` parity check matrix `H` should correspond to the `Z` (or `X`) stabilisers of a CSS QECC with `n` qubits and `m` `Z` (or `X`) stabilisers.
