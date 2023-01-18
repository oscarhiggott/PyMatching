import os
import sys
from contextlib import contextmanager

import pymatching
from pymatching import cli
from pymatching._cli_argv import cli_argv

from .config import DATA_DIR


@contextmanager
def three_errors_data(in_fmt: str):
    out_fn = os.path.join(DATA_DIR, "three_errors_predictions.dets")
    if os.path.isfile(out_fn):
        os.remove(out_fn)
    assert not os.path.exists(out_fn)
    args = [
        "predict",
        "--dem", os.path.join(DATA_DIR, "three_errors.dem"),
        "--in", os.path.join(DATA_DIR, "three_errors." + in_fmt),
        "--in_format", in_fmt,
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
    with three_errors_data("dets") as args:
        cli(command_line_args=args)


def test_protected_cli():
    with three_errors_data("dets") as args:
        pymatching._cpp_pymatching.main(command_line_args=args)


def test_protected_cli_b8_in():
    with three_errors_data("b8") as args:
        pymatching._cpp_pymatching.main(command_line_args=args)


def test_cli_argv():
    from unittest.mock import patch
    with three_errors_data("dets") as args:
        with patch.object(sys, 'argv', ["cli"] + args):
            cli_argv()


def test_load_surface_code_b8_cli():
    dem_path = os.path.join(DATA_DIR, "surface_code_rotated_memory_x_13_0.01.dem")
    dets_b8_in_path = os.path.join(DATA_DIR, "surface_code_rotated_memory_x_13_0.01_1000_shots.b8")
    out_fn = os.path.join(DATA_DIR, "surface_code_rotated_memory_x_13_0.01_1000_shots_temp_predictions.b8")

    pymatching._cpp_pymatching.main(command_line_args=[
        "predict",
        "--dem", dem_path,
        "--in", dets_b8_in_path,
        "--in_format", "b8",
        "--out", out_fn,
        "--out_format", "b8",
        "--in_includes_appended_observables"
    ])
    os.remove(out_fn)
