import numpy as np
import pymatching
import pytest


def test_decode_reweight():
    # Simple graph: 0 -- 1 -- 2
    # Edge (0, 1) weight 2
    # Edge (1, 2) weight 2
    m = pymatching.Matching()
    m.add_edge(0, 1, fault_ids={0}, weight=2)
    m.add_edge(1, 2, fault_ids={1}, weight=2)

    # decode([1, 0, 1]) -> detection events 0, 2
    # Should match 0 to 2 via 1. Weight 4.
    res, weight = m.decode(np.array([1, 0, 1]), return_weight=True)
    assert weight == 4.0

    # Reweight (0, 1) to 5. New weight 5+2 = 7.
    reweights = np.array([[0, 1, 5.0]])
    res, weight = m.decode(
        np.array([1, 0, 1]), return_weight=True, edge_reweights=reweights
    )
    assert weight == 7.0

    # Check weights restored
    res, weight = m.decode(np.array([1, 0, 1]), return_weight=True)
    assert weight == 4.0


def test_decode_reweight_boundary():
    m = pymatching.Matching()
    m.add_boundary_edge(0, fault_ids={0}, weight=2)
    m.add_edge(0, 1, fault_ids={1}, weight=3)

    # decode([1, 0]) -> event at 0. Matches to boundary (weight 2).
    res, weight = m.decode(np.array([1, 0]), return_weight=True)
    assert weight == 2.0

    # Reweight boundary edge to 5.
    reweights = np.array([[0, -1, 5.0]])
    res, weight = m.decode(
        np.array([1, 0]), return_weight=True, edge_reweights=reweights
    )
    assert weight == 5.0

    # Restored
    res, weight = m.decode(np.array([1, 0]), return_weight=True)
    assert weight == 2.0


def test_decode_batch_reweight():
    m = pymatching.Matching()
    m.add_edge(0, 1, fault_ids={0}, weight=2)
    m.add_edge(1, 2, fault_ids={1}, weight=2)

    shots = np.array([[1, 0, 1], [1, 0, 1]], dtype=np.uint8)

    # Shot 0: reweight (0, 1) to 5. Expected weight 7.
    # Shot 1: no reweight. Expected weight 4.

    reweights = [np.array([[0, 1, 5.0]]), None]

    preds, weights = m.decode_batch(
        shots, return_weights=True, edge_reweights=reweights
    )
    assert weights[0] == 7.0
    assert weights[1] == 4.0

    # Check restored
    preds, weights = m.decode_batch(shots, return_weights=True)
    assert weights[0] == 4.0
    assert weights[1] == 4.0


def test_decode_batch_reweight_all_same():
    m = pymatching.Matching()
    m.add_edge(0, 1, fault_ids={0}, weight=2)

    shots = np.array([[1, 1], [1, 1]], dtype=np.uint8)
    # Reweight to 5
    rw = np.array([[0, 1, 5.0]])
    reweights = [rw, rw]

    preds, weights = m.decode_batch(
        shots, return_weights=True, edge_reweights=reweights
    )
    assert weights[0] == 5.0
    assert weights[1] == 5.0

    preds, weights = m.decode_batch(shots, return_weights=True)
    assert weights[0] == 2.0
