Installation
===============

PyMatching can be downloaded and installed from `PyPI <https://pypi.org/project/PyMatching/>`_ with the command::

    pip install pymatching

This is the recommended way to install PyMatching since pip will fetch the pre-compiled binaries, rather than building the C++ extension from source on your machine. 
Note that PyMatching requires Python 3.

If instead you would like to install PyMatching from source, clone the repository (using the `--recursive` flag to include the lib/pybind11 submodule) and then use `pip` to install::

    git clone --recursive https://github.com/oscarhiggott/PyMatching.git
    pip install -e ./PyMatching

The installation may take a few minutes since the C++ extension has to be compiled. If you'd also like to run the tests, first install `pytest <https://docs.pytest.org/en/stable/>`_, and then run::

    pytest ./PyMatching/tests
