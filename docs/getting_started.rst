Getting Started
===============

Installing PyMatching
---------------------

PyMatching requires Python 3.6.0 or later. To install PyMatching, 
first clone the repository (you must use the ``--recursive`` flag)::

    git clone --recursive https://github.com/oscarhiggott/PyMatching.git

and then install using ``pip``::

    pip install -e ./PyMatching


Usage
-----

Much of the functionality of PyMatching is available through the 
:py:class:`pymatching.matching.Matching` class, which can be imported in Python with::

    from pymatching import Matching