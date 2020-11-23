

def test_has_version():
    from pymatching._version import __version__
    assert len(__version__) > 0
    assert isinstance(__version__, str)
