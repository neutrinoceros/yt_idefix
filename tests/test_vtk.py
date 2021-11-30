from pathlib import Path

import pytest
from packaging.version import Version

import yt
import yt_idefix
from yt.extensions.idefix.api import IdefixVtkDataset

YT_VERSION = Version(yt.__version__)

DATA_DIR = Path(__file__).parent / "data"
VTK_FILES = (
    # KHI in cartesian coordinates, geom included in datafile, no streching
    (DATA_DIR.joinpath("khi", "data.0100.vtk"), None),
    # contributed by Gaylor; this dataset is vertically stretched
    # (DATA_DIR.joinpath("planet-disk", "data.0034.vtk"), "polar"),
    # contributed by Marc; this dataset is stretched in the radial direction
    (DATA_DIR.joinpath("cataclysmic-disk", "data.0086.vtk"), "polar"),
    # small spherical 3D test case
    (DATA_DIR.joinpath("FargoMHDSpherical", "data.0010.vtk"), "spherical"),
)


@pytest.mark.parametrize("fn", [_[0] for _ in VTK_FILES])
def test_validation(fn):
    assert IdefixVtkDataset._is_valid(fn)


@pytest.mark.parametrize(("fn", "geometry"), VTK_FILES)
def test_data_access(fn, geometry):
    ds = yt_idefix.load(fn, geometry=geometry)
    ad = ds.all_data()
    ad["gas", "density"]


@pytest.mark.parametrize(("fn", "geometry"), VTK_FILES)
def test_load(fn, geometry):
    ds = yt_idefix.load(fn, geometry=geometry)
    if geometry is None:
        assert ds.geometry in ("cartesian", "polar", "spherical")
    else:
        assert ds.geometry == geometry


# TODO: make this a pytest-mpl test
@pytest.mark.parametrize(("fn", "geometry"), VTK_FILES)
def test_slice_plot(fn, geometry):
    ds = yt_idefix.load(fn, geometry=geometry, unit_system="code")
    if YT_VERSION <= Version("4.0.1"):
        if ds.geometry == "spherical":
            normal = "phi"
        else:
            normal = "z"
        yt.SlicePlot(ds, normal=normal, fields=("gas", "density"))
    else:
        # this should work but it's broken with yt 4.0.1
        # it is fixed in https://github.com/yt-project/yt/pull/3489
        yt.SlicePlot(ds, normal=(0, 0, 1), fields=("gas", "density"))


@pytest.mark.xfail(
    YT_VERSION < Version("4.0.2"),
    reason="all .vtk files are detected as Athena-compatible, "
    "hence dataformat is considered ambiguous before yt 4.0.2",
)
@pytest.mark.parametrize(("fn", "geometry"), VTK_FILES)
def test_load_magic(fn, geometry):
    ds = yt.load(fn, geometry=geometry)
    assert isinstance(ds, IdefixVtkDataset)


@pytest.mark.parametrize(("fn", "geometry"), VTK_FILES)
def test_load_stretched(fn, geometry):
    yt_idefix.load_stretched(fn, geometry=geometry)
