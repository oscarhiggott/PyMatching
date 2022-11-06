import os
import sys
from contextlib import contextmanager

import pymatching
from pymatching import cli
from pymatching._cli_argv import cli_argv

from .config import DATA_DIR


@contextmanager
def three_errors_data():
    out_fn = os.path.join(DATA_DIR, "three_errors_predictions.dets")
    if os.path.isfile(out_fn):
        os.remove(out_fn)
    assert not os.path.exists(out_fn)
    args = [
        "predict",
        "--dem", os.path.join(DATA_DIR, "three_errors.dem"),
        "--in", os.path.join(DATA_DIR, "three_errors.dets"),
        "--in_format", "dets",
        "--out", out_fn,
        "--out_format", "dets",
    ]
    yield args
    assert os.path.isfile(out_fn)
    with open(out_fn) as f:
        assert f.read() == """shot
shot L0
shot L2
shot L1
"""
    os.remove(out_fn)


def test_cli():
    with three_errors_data() as args:
        cli(command_line_args=args)


def test_protected_cli():
    with three_errors_data() as args:
        pymatching._cpp_pymatching.main(command_line_args=args)


def test_cli_argv():
    from unittest.mock import patch
    with three_errors_data() as args:
        with patch.object(sys, 'argv', ["cli"] + args):
            cli_argv()
