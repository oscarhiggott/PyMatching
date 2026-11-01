"""fastgap: PyMatching decodes with a certified lower bound on the complementary gap.

fastgap runs the ordinary, untruncated PyMatching (sparse blossom) decode and then computes a
lower bound on the complementary gap ``W^c - W*`` from the decoder's own final dual state, with a
dense Dijkstra over reduced costs priced by landmark potentials. It never
runs a second decode and never sweeps the empty spacetime volume.

Example::

    import stim, fastgap
    circuit = stim.Circuit.generated("surface_code:rotated_memory_x", distance=5, rounds=5,
                                     after_clifford_depolarization=0.005)
    dem = circuit.detector_error_model(decompose_errors=True)
    idx = fastgap.GapIndex.from_dem(dem)
    dec = fastgap.GapDecoder(idx)
    shots = circuit.compile_detector_sampler().sample(1000, bit_packed=True)
    res = dec.decode_batch(shots)
    res.gap_lb      # certified lower bounds on the gap, in nats
"""

from fastgap import _cpp_fastgap
from fastgap._cpp_fastgap import BoundarylessLogicalError, InvariantViolation, ParallelEdgeConflictError
from fastgap.decoder import BatchResult, GapDecoder, GapResult, pack_shots
from fastgap.exact import ExactBatch, brute_force_gap, exact_gap_batch
from fastgap.index import GapIndex
from fastgap.landmarks import InteriorGraph, LandmarkSource, LandmarkStrategy

__version__ = "0.1.0"

__all__ = [
    "BatchResult",
    "BoundarylessLogicalError",
    "ExactBatch",
    "GapDecoder",
    "GapIndex",
    "GapResult",
    "InteriorGraph",
    "InvariantViolation",
    "LandmarkSource",
    "LandmarkStrategy",
    "ParallelEdgeConflictError",
    "brute_force_gap",
    "exact_gap_batch",
    "pack_shots",
]
