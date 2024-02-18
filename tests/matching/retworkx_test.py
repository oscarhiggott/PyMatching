import pytest

from pymatching import Matching


def test_load_from_retworkx_deprecated():
    rx = pytest.importorskip("rustworkx")
    with pytest.deprecated_call():
        g = rx.PyGraph()
        g.add_nodes_from([{} for _ in range(3)])
        g.add_edge(0, 1, dict(fault_ids=0))
        g.add_edge(0, 2, dict(fault_ids=1))
        g.add_edge(1, 2, dict(fault_ids=2))
        m = Matching()
        m.load_from_retworkx(g)


def test_to_retworkx_deprecated():
    _ = pytest.importorskip("rustworkx")
    with pytest.deprecated_call():
        m = Matching()
        m.add_edge(0, 1, {0})
        m.add_edge(1, 2, {1})
        m.to_retworkx()
