import sys
from pathlib import Path
from typing import List
from unittest.mock import patch

import pytest

import pymatching
from pymatching import cli
from pymatching._cli_argv import cli_argv


def predict_args(dem_file: Path, input_file: Path, output_file: Path) -> List[str]:
    return [
        "predict",
        "--dem", str(dem_file),
        "--in", str(input_file),
        "--in_format", input_file.suffix[1:],
        "--out", str(output_file),
        "--out_format", "dets",
    ]


@pytest.mark.parametrize("cli_function", [cli, pymatching._cpp_pymatching.main])
@pytest.mark.parametrize("input_format", ["dets", "b8"])
def test_calling_cli_creates_expected_file(tmp_path, data_dir, cli_function, input_format):
    output_file = tmp_path / "three_errors_predictions.dets"
    cli_function(command_line_args=predict_args(
        data_dir / "three_errors.dem",
        data_dir / f"three_errors.{input_format}",
        output_file))
    assert output_file.is_file()
    with open(output_file) as prediction_file:
        assert prediction_file.readlines() == ["shot\n", "shot L0\n", "shot L2\n", "shot L1\n"]


@pytest.mark.parametrize("input_format", ["dets", "b8"])
def test_patching_cli_argv_creates_expected_file(tmp_path, data_dir, input_format):
    output_file = tmp_path / "three_errors_predictions.dets"
    args = predict_args(
        data_dir / "three_errors.dem",
        data_dir / f"three_errors.{input_format}",
        output_file)
    with patch.object(sys, "argv", ["cli"] + args):
        cli_argv()
        assert output_file.is_file()
        with open(output_file) as prediction_file:
            assert prediction_file.readlines() == ["shot\n", "shot L0\n", "shot L2\n", "shot L1\n"]


def test_load_surface_code_b8_cli(tmp_path, data_dir):
    dem_path = data_dir / "surface_code_rotated_memory_x_13_0.01.dem"
    dets_b8_in_path = data_dir / "surface_code_rotated_memory_x_13_0.01_1000_shots.b8"
    out_fn = tmp_path / "surface_code_rotated_memory_x_13_0.01_1000_shots_temp_predictions.b8"

    pymatching._cpp_pymatching.main(command_line_args=[
        "predict",
        "--dem", str(dem_path),
        "--in", str(dets_b8_in_path),
        "--in_format", "b8",
        "--out", str(out_fn),
        "--out_format", "b8",
        "--in_includes_appended_observables"
    ])
