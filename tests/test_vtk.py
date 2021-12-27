import os
from pathlib import Path
from typing import Any, Dict

import pytest
from packaging.version import Version
from unyt import assert_allclose_units, unyt_quantity

import yt
import yt_idefix
from yt_idefix.api import IdefixVtkDataset, PlutoVtkDataset

DatabaseT = Dict[str, Dict[str, Any]]

YT_VERSION = Version(yt.__version__)

DATA_DIR = Path(__file__).parent / "data"

IDEFIX_VTK_FILES: DatabaseT = {
    # KHI in cartesian coordinates, geom included in datafile, no streching
    "khi": {
        "path": DATA_DIR.joinpath("khi", "data.0100.vtk"),
        "has_geometry": True,
        "geometry": "cartesian",
        "dimensionality": 2,
        "kind": "idefix",
    },
    # contributed by Gaylor Wafflard Fernandez; this dataset stretched in z and R
    "planet_disk": {
        "path": DATA_DIR.joinpath("planet-disk", "data.0034.vtk"),
        "has_geometry": False,
        "geometry": "polar",
        "dimensionality": 3,
        "kind": "idefix",
    },
    # contributed by Marc Van den Bossche; this dataset is stretched in the radial direction
    "cataclysmic_disk": {
        "path": DATA_DIR.joinpath("cataclysmic-disk", "data.0086.vtk"),
        "has_geometry": False,
        "geometry": "polar",
        "dimensionality": 2,
        "kind": "idefix",
    },
    # small spherical 3D test case
    "spherical_disk": {
        "path": DATA_DIR.joinpath("FargoMHDSpherical", "data.0010.vtk"),
        "has_geometry": True,
        "geometry": "spherical",
        "dimensionality": 3,
        "kind": "idefix",
    },
}

PLUTO_VTK_FILES: DatabaseT = {
    "sedov": {
        "path": DATA_DIR.joinpath("pluto_sedov", "data.0018.vtk"),
        "has_geometry": True,
        "geometry": "cartesian",
        "dimensionality": 3,
        "kind": "pluto",
        "units": {
            "length": unyt_quantity(3.08567758e21, "cm"),
            "velocity": unyt_quantity(1e5, "cm/s"),
            "mass": unyt_quantity(4.91416082e40, "g"),
            "time": unyt_quantity(3.08567758e16, "s"),
            "magnetic": unyt_quantity(4.58462477e-07, "G"),
        },
    },
    "sod": {
        "has_geometry": True,
        "path": DATA_DIR.joinpath("pluto_sod", "data.0001.vtk"),
        "geometry": "cartesian",
        "dimensionality": 1,
        "kind": "pluto",
        "current_time": None,
    },
}

VTK_FILES = {**IDEFIX_VTK_FILES, **PLUTO_VTK_FILES}

# subset of files in which geometry isn't parsable and needs to be provided by the user
VTK_FILES_NO_GEOMETRY = {
    k: v for k, v in VTK_FILES.items() if v["has_geometry"] is False
}

# and the complementary
VTK_FILES_WITH_GEOMETRY = {
    k: v for k, v in VTK_FILES.items() if v["has_geometry"] is True
}


@pytest.mark.parametrize("file", VTK_FILES.values(), ids=VTK_FILES.keys())
def test_class_validation(file):
    cls = {
        "idefix": IdefixVtkDataset,
        "pluto": PlutoVtkDataset,
    }[file["kind"]]
    assert cls._is_valid(file["path"])


@pytest.mark.parametrize("file", PLUTO_VTK_FILES.values(), ids=PLUTO_VTK_FILES.keys())
def test_parse_pluto_metadata(file):
    ds = yt_idefix.load(file["path"])
    assert os.path.samefile(
        ds._definitions_header, file["path"].parent / "definitions.h"
    )
    assert ds.geometry == file["geometry"]
    if "current_time" in file:
        if file["current_time"] is None:
            assert ds.current_time == -ds.time_unit
        else:
            assert ds.current_time == file["current_time"]
    else:
        assert ds.current_time >= 0


@pytest.mark.parametrize("file", PLUTO_VTK_FILES.values(), ids=PLUTO_VTK_FILES.keys())
def test_pluto_units(file):
    ds = yt_idefix.load(file["path"])
    if "units" not in file:
        return
    for u, expected in file["units"].items():
        assert_allclose_units(getattr(ds, f"{u}_unit"), expected)


def test_pluto_wrong_definitions_header():
    with pytest.raises(
        FileNotFoundError,
        match=(
            "Header file definitions2.h couldn't be found. "
            "The 'geometry' keyword argument must be specified."
        ),
    ):
        yt_idefix.load(
            PLUTO_VTK_FILES["sedov"]["path"], definitions_header="definitions2.h"
        )


def test_pluto_wrong_definitions_header_with_geometry():
    with pytest.warns(
        UserWarning,
        match=(
            "Header file definitions2.h couldn't be found. "
            "The code units are set to be 1.0 in cgs by default."
        ),
    ):
        yt_idefix.load(
            PLUTO_VTK_FILES["sedov"]["path"],
            definitions_header="definitions2.h",
            geometry="cartesian",
        )


@pytest.mark.parametrize("file", VTK_FILES.values(), ids=VTK_FILES.keys())
def test_data_access(file):
    ds = yt_idefix.load(file["path"], geometry=file["geometry"])
    ad = ds.all_data()
    ad["gas", "density"]


@pytest.mark.parametrize("file", VTK_FILES.values(), ids=VTK_FILES.keys())
def test_load_user_passed_geometry(file):
    ds = yt_idefix.load(file["path"], geometry=file["geometry"])
    assert ds.geometry == file["geometry"]
    assert ds.dimensionality == file["dimensionality"]


@pytest.mark.parametrize(
    "file", VTK_FILES_NO_GEOMETRY.values(), ids=VTK_FILES_NO_GEOMETRY.keys()
)
def test_load_missing_geometry_metadata(file):
    with pytest.raises(
        ValueError,
        match=(
            "Geometry couldn't be parsed from file. "
            "The 'geometry' keyword argument must be specified."
        ),
    ):
        yt_idefix.load(file["path"])


@pytest.mark.parametrize(
    "file", VTK_FILES_WITH_GEOMETRY.values(), ids=VTK_FILES_WITH_GEOMETRY.keys()
)
def test_load_parseable_geometry_metadata(file):
    ds = yt_idefix.load(file["path"])
    assert ds.geometry == file["geometry"]


# TODO: make this a pytest-mpl test
@pytest.mark.parametrize("file", VTK_FILES.values(), ids=VTK_FILES.keys())
def test_slice_plot(file):
    ds = yt_idefix.load(file["path"], geometry=file["geometry"], unit_system="code")
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
@pytest.mark.parametrize("file", VTK_FILES.values(), ids=VTK_FILES.keys())
def test_load_magic(file):
    ds = yt.load(file["path"], geometry=file["geometry"])
    assert isinstance(ds, IdefixVtkDataset)


@pytest.mark.parametrize("file", VTK_FILES.values(), ids=VTK_FILES.keys())
def test_load_stretched(file):
    yt_idefix.load_stretched(file["path"], geometry=file["geometry"])
