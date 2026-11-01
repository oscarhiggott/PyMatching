"""GapDecoder: MWPM decode plus a certified lower bound on the complementary gap."""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np

from fastgap import _cpp_fastgap
from fastgap.index import GapIndex


def pack_shots(shots: np.ndarray, num_detectors: int, bit_packed: bool) -> np.ndarray:
    """Returns a C-contiguous (num_shots, ceil(num_detectors / 8)) uint8 bit-packed array."""
    shots = np.asarray(shots)
    if shots.ndim == 1:
        shots = shots.reshape(1, -1)
    if bit_packed:
        return np.ascontiguousarray(shots, dtype=np.uint8)
    if shots.shape[1] != num_detectors:
        raise ValueError(f"unpacked shots must have {num_detectors} columns, got {shots.shape[1]}")
    return np.ascontiguousarray(np.packbits(shots.astype(np.uint8), axis=1, bitorder="little"))


def _threshold_int(idx: GapIndex, threshold: Optional[float]) -> int:
    # gap_int >= ceil(tau * nc) implies gap >= tau, so post-selection decisions stay sound.
    if threshold is None:
        return -1
    if threshold < 0:
        raise ValueError("threshold must be non-negative")
    return int(math.ceil(float(threshold) * idx.normalising_constant))


def _prediction_bits(masks: np.ndarray, num_observables: int) -> np.ndarray:
    masks = np.asarray(masks, dtype=np.uint64)
    bits = (masks[..., None] >> np.arange(num_observables, dtype=np.uint64)) & np.uint64(1)
    return bits.astype(np.uint8)


@dataclass
class GapResult:
    """Result of one shot. Weights and gaps are in user units (nats); ``*_int`` fields are in
    PyMatching's internal integer units, which is what soundness comparisons must use."""

    prediction: np.ndarray
    weight: float
    gap_lb: float
    gap_lb_db: float
    gap_ub: float
    censored: bool
    walk_simple: bool
    all_hops_exact: bool
    blossom_pass: bool
    num_defects: int
    num_blossoms: int
    num_hops: int
    hop_slack: float
    t_decode_ns: int
    t_gap_ns: int
    weight_int: int
    gap_lb_int: float
    gap_ub_int: float
    walk: List[Tuple[int, int, int, bool]] = field(default_factory=list)

    @property
    def gap_ub_db(self) -> float:
        return GapIndex.to_db(self.gap_ub)


@dataclass
class BatchResult:
    """Arrays of per-shot results (same fields as :class:`GapResult`)."""

    prediction: np.ndarray
    weight: np.ndarray
    gap_lb: np.ndarray
    gap_lb_db: np.ndarray
    gap_ub: np.ndarray
    censored: np.ndarray
    walk_simple: np.ndarray
    all_hops_exact: np.ndarray
    blossom_pass: np.ndarray
    num_defects: np.ndarray
    num_blossoms: np.ndarray
    num_hops: np.ndarray
    hop_slack: np.ndarray
    t_decode_ns: np.ndarray
    t_gap_ns: np.ndarray
    weight_int: np.ndarray
    gap_lb_int: np.ndarray
    gap_ub_int: np.ndarray

    def __len__(self) -> int:
        return len(self.weight)

    @property
    def gap_ub_db(self) -> np.ndarray:
        return GapIndex.to_db(self.gap_ub)

    def classify(self, threshold: float) -> np.ndarray:
        """Post-selection at threshold tau (CLAUDE.md §3.5): +1 keep (gap_lb >= tau), -1 discard
        (gap_ub < tau), 0 uncertain."""
        out = np.zeros(len(self), dtype=np.int8)
        out[self.gap_lb >= threshold] = 1
        out[(out == 0) & (self.gap_ub < threshold)] = -1
        return out


class GapDecoder:
    """Decodes shots with PyMatching's sparse blossom and returns a certified lower bound on the
    complementary gap from the decoder's final dual state (no second decode).

    Args:
        index: the :class:`GapIndex` for the DEM.
        num_threads: number of threads that cooperate on *each* shot's gap search, including the
            calling thread (which runs the search; the other ``num_threads - 1`` prefetch price
            rows). Created once here, never per shot. Shots are always decoded one after another, so this
            lowers per-shot latency rather than raising throughput. Shots with fewer than
            ``parallel_grain`` defects run on one thread. Use at most the size of the fastest
            core cluster (5 on an M5 Pro); beyond that latency gets worse.
        pin_threads: pin the threads to cores where the OS allows (Linux; on macOS they only ask
            for the user-interactive QoS class).
        check: run the debug invariants (dual feasibility, relay tightness, d_hat <= dist) on
            every shot. Slow; meant for tests.
    """

    def __init__(self, index: GapIndex, num_threads: int = 1, pin_threads: bool = False, check: bool = False):
        self.index = index
        self.check = bool(check)
        self._core = _cpp_fastgap.GapDecoderCore(index._core, int(num_threads), bool(pin_threads))

    @property
    def num_threads(self) -> int:
        return int(self._core.num_threads)

    @property
    def parallel_grain(self) -> int:
        """Shots with fewer defects than this use one thread (results never depend on it)."""
        return int(self._core.parallel_grain)

    @parallel_grain.setter
    def parallel_grain(self, value: int) -> None:
        self._core.parallel_grain = int(value)

    def decode(
        self,
        syndrome=None,
        threshold: Optional[float] = None,
        upper_bound: bool = False,
        *,
        detection_events=None,
        pricing: str = "index",
        keep_walk: bool = False,
    ) -> GapResult:
        """Decodes one shot.

        Args:
            syndrome: dense binary array with one entry per detector (like
                ``pymatching.Matching.decode``). Alternatively pass ``detection_events``.
            threshold: early-termination threshold tau in user units (CLAUDE.md §3.6). The
                search stops once every open label is >= tau and reports ``gap_lb = tau`` with
                ``censored=True``.
            upper_bound: also compute ``gap_ub`` (CLAUDE.md §3.5).
            pricing: ``"index"`` (landmark bound: the reported bound) or ``"exact"`` (exact
                reduced costs via Dijkstra: the slow reference Delta_sigma).
            keep_walk: return the hops of the walk in ``walk`` as (from, to, price_int, exact),
                with -1 for L and -2 for R.
        """
        if detection_events is None:
            if syndrome is None:
                raise ValueError("pass a syndrome or detection_events")
            syndrome = np.asarray(syndrome).reshape(-1)
            if syndrome.shape[0] != self.index.num_detectors:
                raise ValueError(f"syndrome must have {self.index.num_detectors} entries, got {syndrome.shape[0]}")
            detection_events = np.flatnonzero(syndrome)
        events = np.asarray(detection_events, dtype=np.uint64).reshape(-1)
        r = self._core.decode(
            events,
            threshold=_threshold_int(self.index, threshold),
            upper_bound=bool(upper_bound),
            pricing=pricing,
            check=self.check,
            keep_walk=bool(keep_walk),
        )
        idx = self.index
        gap_lb = idx.to_user_units(int(r["gap_lb"]))
        gap_ub = idx.to_user_units(int(r["gap_ub"]))
        slack = int(r["hop_slack"])
        return GapResult(
            prediction=_prediction_bits(np.uint64(r["prediction_mask"]), idx.num_observables),
            weight=int(r["w_star"]) / idx.normalising_constant,
            gap_lb=gap_lb,
            gap_lb_db=GapIndex.to_db(gap_lb),
            gap_ub=gap_ub,
            censored=bool(r["censored"]),
            walk_simple=bool(r["walk_simple"]),
            all_hops_exact=bool(r["all_hops_exact"]),
            blossom_pass=bool(r["blossom_pass"]),
            num_defects=int(r["num_defects"]),
            num_blossoms=int(r["num_blossoms"]),
            num_hops=int(r["num_hops"]),
            hop_slack=np.nan if slack < 0 else slack / idx.normalising_constant,
            t_decode_ns=int(r["t_decode_ns"]),
            t_gap_ns=int(r["t_gap_ns"]),
            weight_int=int(r["w_star"]),
            gap_lb_int=np.inf if int(r["gap_lb"]) < 0 else int(r["gap_lb"]),
            gap_ub_int=np.inf if int(r["gap_ub"]) < 0 else int(r["gap_ub"]),
            walk=[(int(a), int(b), int(p), bool(e)) for a, b, p, e in r["walk"]],
        )

    def decode_batch(
        self,
        shots,
        threshold: Optional[float] = None,
        upper_bound: bool = False,
        *,
        bit_packed_shots: bool = True,
        pricing: str = "index",
    ) -> BatchResult:
        """Decodes many shots. ``shots`` is (num_shots, ceil(num_detectors/8)) bit-packed uint8
        (little-endian bit order, as from ``stim``'s ``bit_packed=True``), or set
        ``bit_packed_shots=False`` for a dense (num_shots, num_detectors) array."""
        idx = self.index
        packed = pack_shots(shots, idx.num_detectors, bit_packed_shots)
        r = self._core.decode_batch(
            packed,
            threshold=_threshold_int(idx, threshold),
            upper_bound=bool(upper_bound),
            pricing=pricing,
            check=self.check,
        )
        gap_lb = idx.to_user_units(r["gap_lb"])
        gap_ub = idx.to_user_units(r["gap_ub"])
        slack = r["hop_slack"].astype(np.float64)
        slack = np.where(slack < 0, np.nan, slack / idx.normalising_constant)
        return BatchResult(
            prediction=_prediction_bits(r["prediction_mask"], idx.num_observables),
            weight=r["w_star"] / idx.normalising_constant,
            gap_lb=gap_lb,
            gap_lb_db=GapIndex.to_db(gap_lb),
            gap_ub=gap_ub,
            censored=r["censored"],
            walk_simple=r["walk_simple"],
            all_hops_exact=r["all_hops_exact"],
            blossom_pass=r["blossom_pass"],
            num_defects=r["num_defects"],
            num_blossoms=r["num_blossoms"],
            num_hops=r["num_hops"],
            hop_slack=slack,
            t_decode_ns=r["t_decode_ns"],
            t_gap_ns=r["t_gap_ns"],
            weight_int=r["w_star"],
            gap_lb_int=np.where(r["gap_lb"] < 0, np.inf, r["gap_lb"].astype(np.float64)),
            gap_ub_int=np.where(r["gap_ub"] < 0, np.inf, r["gap_ub"].astype(np.float64)),
        )
