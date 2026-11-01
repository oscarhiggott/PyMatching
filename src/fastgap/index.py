"""GapIndex: per-DEM precomputation (canonical graph, boundary and landmark potentials)."""

from __future__ import annotations

import json
import logging
from typing import Any, Dict, List, Optional, Union

import numpy as np

from fastgap import _cpp_fastgap
from fastgap.landmarks import (
    InteriorGraph,
    LandmarkSource,
    LandmarkStrategy,
    dedupe,
    resolve,
    strategy_params,
)

LOG10_E_TIMES_10 = 10.0 * np.log10(np.e)

logger = logging.getLogger("fastgap")


def _coords_of(dem) -> Optional[Dict[int, List[float]]]:
    getter = getattr(dem, "get_detector_coordinates", None)
    if getter is None:
        return None
    try:
        return {int(k): list(v) for k, v in getter().items()}
    except Exception:  # pragma: no cover - coordinates are optional
        return None


class GapIndex:
    """Precomputed data for one (DEM, observable) pair (CLAUDE.md §5.3).

    Build with :meth:`from_dem` (or :meth:`from_gap_dem` for an external reference DEM). The
    index stores the gauge bits, the boundary potentials ``h_L``/``h_R``, the landmark potentials
    and the landmark strategy config with its hash. Cache indices with :meth:`save`/:meth:`load`.
    """

    def __init__(self, core):
        self._core = core
        self._exact_cache: Dict[int, Any] = {}
        self._graph: Optional[InteriorGraph] = None

    # ------------------------------------------------------------------ construction

    @classmethod
    def from_dem(
        cls,
        dem,
        observable: int = 0,
        strategy: Union[str, LandmarkStrategy] = "auto",
        num_landmarks: int = 16,
        meta: Optional[Dict[str, Any]] = None,
    ) -> "GapIndex":
        """Builds an index from a graphlike stim DEM (``decompose_errors=True``).

        Args:
            dem: a ``stim.DetectorErrorModel`` (or its text).
            observable: the observable index k whose gap is computed.
            strategy: landmark strategy name (see ``fastgap.landmarks.STRATEGIES``), ``"auto"``,
                or a ``LandmarkStrategy`` object.
            num_landmarks: landmark budget K (L and R come on top).
            meta: extra information for the strategy.

        Raises:
            fastgap.BoundarylessLogicalError: observable k has a logical that never touches a
                boundary (e.g. a toroidal code).
            fastgap.ParallelEdgeConflictError: two mechanisms on the same detectors differ in
                the gauge-fixed observable (a local weight-2 logical).
        """
        text = dem if isinstance(dem, str) else str(dem)
        core = _cpp_fastgap.GapIndexCore.create_from_dem(text, int(observable))
        coords = None if isinstance(dem, str) else _coords_of(dem)
        meta = dict(meta or {})
        if not isinstance(dem, str):
            meta.setdefault("dem", dem)
        meta.setdefault("observable", int(observable))
        meta.setdefault("num_landmarks", int(num_landmarks))
        idx = cls(core)
        idx._place_landmarks(strategy, coords, meta, num_landmarks)
        return idx

    @classmethod
    def from_gap_dem(
        cls,
        gap_dem,
        obs_detector: Optional[int] = None,
        strategy: Union[str, LandmarkStrategy] = "auto",
        num_landmarks: int = 16,
        meta: Optional[Dict[str, Any]] = None,
    ) -> "GapIndex":
        """Builds an index from a reference DEM whose observable is already a detector.

        This is the path for Gidney's desaturation ``gap_dem`` (CLAUDE.md §7.3): fastgap's graph
        is built from exactly the edges PyMatching keeps for ``gap_dem``. Edges to the
        observable detector become L boundary edges, ordinary boundary edges become R edges, and
        the observable detector itself is dropped (its events are ignored). The edge set, weights
        and normalising constant are then checked against the reference graph.
        """
        import pymatching

        matching = pymatching.Matching.from_detector_error_model(gap_dem)
        num_detectors = int(gap_dem.num_detectors)
        if obs_detector is None:
            obs_detector = num_detectors - 1
        edges = []
        for u, v, attrs in matching.edges():
            w = float(attrs["weight"])
            p = float(attrs.get("error_probability", -1.0))
            if v is None:
                if u == obs_detector:
                    edges.append((-1, -1, False, w, p))
                else:
                    edges.append((u, -1, False, w, p))
            elif u == obs_detector or v == obs_detector:
                other = v if u == obs_detector else u
                edges.append((other, -1, True, w, p))
            else:
                edges.append((u, v, False, w, p))
        core = _cpp_fastgap.GapIndexCore.create_from_reference_edges(num_detectors, int(obs_detector), edges)
        mismatch = core.compare_with_reference_dem(str(gap_dem))
        if mismatch:
            raise AssertionError(f"fastgap graph does not match the reference gap_dem graph: {mismatch}")
        coords = _coords_of(gap_dem)
        if coords is not None:
            coords.pop(int(obs_detector), None)
        meta = dict(meta or {})
        meta.setdefault("num_landmarks", int(num_landmarks))
        idx = cls(core)
        idx._place_landmarks(strategy, coords, meta, num_landmarks)
        return idx

    def _place_landmarks(self, strategy, coords, meta, k):
        graph = self.graph
        strat = resolve(strategy, graph, coords, meta)
        sources: List[LandmarkSource] = dedupe(list(strat.select(graph, coords, meta, int(k))))[: int(k)]
        sources = [LandmarkSource(frozenset(int(x) for x in s.nodes), str(s.label)) for s in sources]
        payload = [(s.label, sorted(s.nodes)) for s in sources]
        config = json.dumps(
            {
                "strategy": getattr(strat, "name", type(strat).__name__),
                "params": strategy_params(strat),
                "num_landmarks": int(k),
                "observable": int(self.observable),
                "sources": payload,
            },
            sort_keys=True,
        )
        self._core.set_landmarks(payload, config)
        mem = self.memory_report()
        logger.info(
            "fastgap index: %d detectors, strategy %s (%d landmarks, hash %s); potentials %.0f B/detector",
            mem["num_detectors"], self.strategy_name, mem["num_landmarks"], self.strategy_hash,
            mem["potential_bytes_per_detector"],
        )

    # ------------------------------------------------------------------ persistence

    def save(self, path) -> None:
        self._core.save(str(path))

    @classmethod
    def load(cls, path) -> "GapIndex":
        return cls(_cpp_fastgap.GapIndexCore.load(str(path)))

    # ------------------------------------------------------------------ properties

    @property
    def graph(self) -> InteriorGraph:
        if self._graph is None:
            self._graph = InteriorGraph(self._core)
        return self._graph

    @property
    def num_detectors(self) -> int:
        return int(self._core.num_detectors)

    @property
    def num_observables(self) -> int:
        return int(self._core.num_observables)

    @property
    def observable(self) -> int:
        return int(self._core.observable)

    @property
    def normalising_constant(self) -> float:
        """Internal units per user unit (``2c`` in CLAUDE.md §2.1)."""
        return float(self._core.normalising_constant)

    @property
    def d_lr(self) -> float:
        """Weight of the lightest logical (user units)."""
        return self.to_user_units(self._core.d_lr)

    @property
    def d_lr_int(self) -> int:
        return int(self._core.d_lr)

    @property
    def strategy_config(self) -> Dict[str, Any]:
        return json.loads(self._core.strategy_config)

    @property
    def strategy_name(self) -> str:
        return self.strategy_config.get("strategy", "")

    @property
    def strategy_hash(self) -> str:
        return f"{int(self._core.strategy_hash):016x}"

    @property
    def landmarks(self) -> List[LandmarkSource]:
        return [LandmarkSource(frozenset(int(x) for x in nodes), label) for label, nodes in self._core.landmarks]

    @property
    def gauge(self) -> np.ndarray:
        return self._core.gauge

    @property
    def potential_bytes(self) -> int:
        return int(self._core.potential_bytes)

    def memory_report(self) -> Dict[str, float]:
        """Potential memory (CLAUDE.md §7.4)."""
        n = max(1, self.num_detectors)
        return {
            "num_detectors": self.num_detectors,
            "num_landmarks": len(self.landmarks),
            "potential_bytes": self.potential_bytes,
            "potential_bytes_per_detector": self.potential_bytes / n,
        }

    def canonical_dem(self):
        """The canonical (gauge-fixed) DEM as a ``stim.DetectorErrorModel``."""
        import stim

        return stim.DetectorErrorModel(self._core.canonical_dem_text())

    def to_user_units(self, x):
        """Internal integers to user units; DIST_INF (or -1 from batch arrays) maps to inf."""
        arr = np.asarray(x, dtype=np.float64)
        inf_mask = (arr < 0) | (arr >= _cpp_fastgap.DIST_INF // 2)
        out = np.where(inf_mask, np.inf, arr / self.normalising_constant)
        return float(out) if out.ndim == 0 else out

    @staticmethod
    def to_db(x):
        """Gap in decibels: 10 log10(e) * gap (CLAUDE.md §2.1)."""
        return np.asarray(x, dtype=np.float64) * LOG10_E_TIMES_10 if np.ndim(x) else float(x) * LOG10_E_TIMES_10

    def __repr__(self) -> str:
        return (
            f"fastgap.GapIndex(num_detectors={self.num_detectors}, observable={self.observable}, "
            f"strategy={self.strategy_name!r}, num_landmarks={len(self.landmarks)}, "
            f"hash={self.strategy_hash})"
        )
