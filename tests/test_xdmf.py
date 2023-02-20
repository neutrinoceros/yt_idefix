from pathlib import Path

import pytest

import yt
from yt_idefix.api import PlutoXdmfDataset

pytest.importorskip("h5py")


HERE = Path(__file__).parent
DATADIR = HERE / "data"
ds_path_plutoXDMF = DATADIR.joinpath("pluto_isentropic_vortex", "data.0010.flt.h5")


def test_class_validation(xdmf_file):
    file = xdmf_file
    assert PlutoXdmfDataset._is_valid(str(file["path"]))


def test_load_class(xdmf_file):
    file = xdmf_file
    PlutoXdmfDataset(str(file["path"]))


# TODO: make this a pytest-mpl test
def test_slice_plot(xdmf_file):
    file = xdmf_file
    ds = yt.load(file["path"], geometry=file["geometry"], unit_system="code")
    yt.SlicePlot(ds, normal=(0, 0, 1), fields=("gas", "density"))


def test_projection_plot(xdmf_file):
    file = xdmf_file
    ds = yt.load(file["path"], geometry=file["geometry"], unit_system="code")
    yt.ProjectionPlot(ds, normal=(0, 0, 1), fields=("gas", "density"))


def test_load_magic(xdmf_file):
    ds = yt.load(xdmf_file["path"], geometry=xdmf_file["geometry"])
    assert isinstance(ds, PlutoXdmfDataset)
