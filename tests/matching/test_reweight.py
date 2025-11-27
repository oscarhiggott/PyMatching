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
    res, weight = m.decode(np.array([1, 0, 1]), return_weight=True, edge_reweights=reweights)
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
    res, weight = m.decode(np.array([1, 0]), return_weight=True, edge_reweights=reweights)
    assert weight == 5.0
    
    # Restored
    res, weight = m.decode(np.array([1, 0]), return_weight=True)
    assert weight == 2.0

def test_decode_batch_reweight():
    m = pymatching.Matching()
    m.add_edge(0, 1, fault_ids={0}, weight=2)
    m.add_edge(1, 2, fault_ids={1}, weight=2)
    
    shots = np.array([
        [1, 0, 1],
        [1, 0, 1]
    ], dtype=np.uint8)
    
    # Shot 0: reweight (0, 1) to 5. Expected weight 7.
    # Shot 1: no reweight. Expected weight 4.
    
    reweights = [
        np.array([[0, 1, 5.0]]),
        None
    ]
    
    preds, weights = m.decode_batch(shots, return_weights=True, edge_reweights=reweights)
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
     
     preds, weights = m.decode_batch(shots, return_weights=True, edge_reweights=reweights)
     assert weights[0] == 5.0
     assert weights[1] == 5.0
     
     preds, weights = m.decode_batch(shots, return_weights=True)
     assert weights[0] == 2.0

def test_decode_reweight_large_observables():
    # If num_observables > 64, the search graph should be present even if enable_correlations=False.
    m = pymatching.Matching()
    # Add enough edges with unique fault_ids to exceed 64 observables
    for i in range(70):
        m.add_edge(i, i+1, fault_ids={i}, weight=1)
    
    assert m.num_fault_ids >= 70
    
    # Decode a simple case: error on edge (0, 1)
    # Expected weight 1.
    syndrome = np.zeros(m.num_nodes, dtype=np.uint8)
    syndrome[0] = 1
    syndrome[1] = 1
    
    res, weight = m.decode(syndrome, return_weight=True)
    assert weight == 1.0
    
    # Reweight edge (0, 1) to 10.
    reweights = np.array([[0, 1, 10.0]])
    res, weight = m.decode(syndrome, return_weight=True, edge_reweights=reweights)
    assert weight == 10.0
    
    # Verify the weight is restored
    res, weight = m.decode(syndrome, return_weight=True)
    assert weight == 1.0

def test_reweight_sign_flip_raises_error():
    m = pymatching.Matching()
    m.add_edge(0, 1, weight=2)
    m.add_edge(1, 2, weight=-3)
    
    # Positive to negative (flip)
    with pytest.raises(ValueError, match="sign flip not allowed"):
        m.decode(np.array([1, 0, 1]), edge_reweights=np.array([[0, 1, -5.0]]))
        
    # Negative to positive (flip)
    with pytest.raises(ValueError, match="sign flip not allowed"):
        m.decode(np.array([1, 0, 1]), edge_reweights=np.array([[1, 2, 3.0]]))

def test_reweight_negative_to_negative():
    # Graph: 0 -- 1 -- 2
    # (0, 1) weight 5
    # (1, 2) weight -3.  
    # Solution for detection events at 0, 2.
    # Standard matching: 0 matches to 2 via 1. Path: (0,1), (1,2).
    # Cost: 5 + (-3) = 2.
    # Note: Negative weight -3 means edge (1,2) is pre-flipped.
    # Events at 0, 2 means syndrome is 1 at 0, 1 at 2.
    # If (1,2) is pre-flipped, it causes events at 1, 2.
    # Observed syndrome: 0:1, 1:0, 2:1.
    # Adjusted syndrome (xor with negative weight syndrome):
    # 0:1, 1:1, 2:0.
    # Now we match 0 and 1. Path (0, 1) cost 5.
    # Total cost = 5 + (-3) = 2.
    
    m = pymatching.Matching()
    m.add_edge(0, 1, weight=5)
    m.add_edge(1, 2, weight=-3)
    
    # Check baseline
    res, weight = m.decode(np.array([1, 0, 1]), return_weight=True)
    assert weight == 2.0
    
    # Reweight (1, 2) to -10.
    # New cost calculation:
    # Path (0, 1) cost 5.
    # Total cost = 5 + (-10) = -5.
    reweights = np.array([[1, 2, -10.0]])
    res, weight = m.decode(np.array([1, 0, 1]), return_weight=True, edge_reweights=reweights)
    assert weight == -5.0
    
    # Verify restoration
    res, weight = m.decode(np.array([1, 0, 1]), return_weight=True)
    assert weight == 2.0
