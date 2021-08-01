import doctest

import pymatching


def test_cpp_docstrings():
    doctest.testmod(pymatching._cpp_mwpm, raise_on_error=True)
