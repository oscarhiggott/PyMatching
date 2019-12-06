import numpy as np

from quantumcode import load_css_from_hdf5
from mwpm.slow_matching import decode


def test_slow_decode():
    fn = "css_toric_[[18,2,3]]_rank_deficient"
    css = load_css_from_hdf5(fn)
    H = css.x_stabilisers
    z = np.zeros(H.shape[0])
    vs = [0,4,3,6]
    for v in vs:
        z[v] = 1
    c = decode(H, z)
    expected = np.array([0., 0., 0., 0., 0., 0., 
                         1., 0., 0., 0., 0., 0., 
                         0., 0., 0., 1., 0., 0.])
    assert(np.array_equal(c, expected))
