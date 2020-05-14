from pymatching import rand_float


def test_rand_float():
    N = 1000
    s = sum(rand_float(0.,1.) for i in range(N))
    assert 430 < s < 570
