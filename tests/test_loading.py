import os
import sys
from pathlib import Path

import pytest

from yt_idefix.api import PlutoVtkDataset

HERE = Path(__file__).parent
DATADIR = HERE / "data"


ds_path = DATADIR.joinpath("pluto_sod", "data.0001.vtk")


def test_load_from_str():
    PlutoVtkDataset(str(ds_path))


def test_load_from_path():
    PlutoVtkDataset(ds_path)


def test_load_from_parent_str():
    # https://github.com/neutrinoceros/yt_idefix/issues/88
    os.chdir(ds_path.parent)
    fn = os.path.join("..", ds_path.parent.name, ds_path.name)
    PlutoVtkDataset(fn)


@pytest.mark.skipif(
    sys.version_info < (3, 9) or not ds_path.is_relative_to(Path.home()),
    reason=(
        "$HOME isn't a parent to the test dataset, "
        "or Python is too old (< 3.9) for us to test that condition easily."
    ),
    # see https://docs.python.org/3/library/pathlib.html#pathlib.PurePath.is_relative_to
)
def test_load_from_home_str():
    # https://github.com/neutrinoceros/yt_idefix/issues/91
    fn = os.path.join("~", ds_path.relative_to(Path.home()))
    PlutoVtkDataset(fn)
