import os
from pathlib import Path

import pytest
from packaging.version import Version

import yt
from yt.extensions.idefix.api import IdefixVtkDataset


def myload(fn, *args, **kwargs):
    # drop in replacement for yt.load, that compensate for the limitations of loading vtk files
    fn = os.fspath(fn)
    if fn.endswith(".vtk"):
        return IdefixVtkDataset(fn, *args, **kwargs)
    else:
        return yt.load(fn, *args, **kwargs)


YT_VERSION = Version(yt.__version__)

DATA_DIR = Path(__file__).parent / "data"
VTK_FILES = (
    # KHI in cartesian coordinates, geom included in datafile, no streching
    (DATA_DIR.joinpath("khi", "data.0100.vtk"), None),
    # contributed by Gaylor; this dataset is vertically stretched
    # (DATA_DIR.joinpath("planet-disk", "data.0034.vtk"), "polar"),
    # contributed by Marc; this dataset is stretched in the radial direction
    (DATA_DIR.joinpath("cataclysmic-disk", "data.0086.vtk"), "polar"),
)


@pytest.mark.parametrize("fn", [_[0] for _ in VTK_FILES])
def test_validation(fn):
    assert IdefixVtkDataset._is_valid(fn)


@pytest.mark.parametrize(("fn", "geometry"), VTK_FILES)
def test_data_access(fn, geometry):
    ds = myload(fn, geometry=geometry)
    ad = ds.all_data()
    ad["gas", "density"]


@pytest.mark.parametrize(("fn", "geometry"), VTK_FILES)
def test_slice_plot(fn, geometry):
    ds = myload(fn, geometry=geometry)
    yt.SlicePlot(ds, "z", ("gas", "density"))


@pytest.mark.parametrize(("fn", "geometry"), VTK_FILES)
def test_load(fn, geometry):
    ds = myload(fn, geometry=geometry)
    if geometry is None:
        assert ds.geometry in ("cartesian", "polar", "spherical")
    else:
        assert ds.geometry == geometry


# TODO: make this a pytest-mpl test
@pytest.mark.parametrize(("fn", "geometry"), VTK_FILES)
def test_plot(fn, geometry):
    ds = myload(fn, geometry=geometry, unit_system="code")
    yt.SlicePlot(ds, "z", ("gas", "density"))
    # p.save("/tmp/")


@pytest.mark.xfail(
    YT_VERSION < Version("4.1"),
    reason="all .vtk files are detected as Athena-compatible, "
    "hence dataformat is considered ambiguous before y 4.1",
)
@pytest.mark.parametrize(("fn", "geometry"), VTK_FILES)
def test_load_magic(fn, geometry):
    ds = yt.load(fn, geometry=geometry)
    assert isinstance(ds, IdefixVtkDataset)
