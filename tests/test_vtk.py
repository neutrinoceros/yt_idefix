import os

import pytest
from packaging.version import Version
from unyt import assert_allclose_units

import yt
import yt_idefix
from yt_idefix.api import IdefixVtkDataset, PlutoVtkDataset

YT_VERSION = Version(yt.__version__)


def test_class_validation(vtk_file):
    file = vtk_file
    cls = {
        "idefix": IdefixVtkDataset,
        "pluto": PlutoVtkDataset,
    }[file["kind"]]
    assert cls._is_valid(file["path"])


def test_parse_pluto_metadata(pluto_vtk_file):
    file = pluto_vtk_file
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


def test_pluto_units(pluto_vtk_file):
    file = pluto_vtk_file
    ds = yt_idefix.load(file["path"])
    if "units" not in file:
        return
    for u, expected in file["units"].items():
        assert_allclose_units(getattr(ds, f"{u}_unit"), expected)


def test_pluto_wrong_definitions_header(pluto_vtk_file):
    with pytest.raises(
        FileNotFoundError,
        match=(
            "Header file definitions2.h couldn't be found. "
            "The 'geometry' keyword argument must be specified."
        ),
    ):
        yt_idefix.load(pluto_vtk_file["path"], definitions_header="definitions2.h")


def test_pluto_wrong_definitions_header_with_geometry(pluto_vtk_file):
    with pytest.warns(
        UserWarning,
        match=(
            "Header file definitions2.h couldn't be found. "
            "The code units are set to be 1.0 in cgs by default."
        ),
    ):
        yt_idefix.load(
            pluto_vtk_file["path"],
            definitions_header="definitions2.h",
            geometry="cartesian",
        )


def test_data_access(vtk_file):
    file = vtk_file
    ds = yt_idefix.load(file["path"], geometry=file["geometry"])
    ad = ds.all_data()
    ad["gas", "density"]


def test_load_user_passed_geometry(vtk_file):
    file = vtk_file
    ds = yt_idefix.load(file["path"], geometry=file["geometry"])
    assert ds.geometry == file["geometry"]
    assert ds.dimensionality == file["dimensionality"]


def test_load_missing_geometry_metadata(vtk_file_no_geom):
    with pytest.raises(
        ValueError,
        match=(
            "Geometry couldn't be parsed from file. "
            "The 'geometry' keyword argument must be specified."
        ),
    ):
        yt_idefix.load(vtk_file_no_geom["path"])


def test_load_parseable_geometry_metadata(vtk_file_with_geom):
    file = vtk_file_with_geom
    ds = yt_idefix.load(file["path"])
    assert ds.geometry == file["geometry"]


def test_derived_field(vtk_file):
    file = vtk_file
    ds = yt_idefix.load(file["path"], geometry=file["geometry"])
    # this derived field is compatible for all current test data
    # it'll be replaced once a better universal field is defined internally.

    def mom_den(field, data):
        return data["idefix-vtk", "RHO"] * data["idefix-vtk", "VX1"]

    ds.add_field(
        name=("gas", "momentum_density"),
        function=mom_den,
        units="g*cm**-2/s",
        sampling_type="cell",
    )
    ad = ds.all_data()
    test = ad["gas", "momentum_density"]
    expect = ad["idefix-vtk", "RHO"] * ad["idefix-vtk", "VX1"]
    assert_allclose_units(test, expect)


# TODO: make this a pytest-mpl test
def test_slice_plot(vtk_file):
    file = vtk_file
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
def test_load_magic(vtk_file):
    ds = yt.load(vtk_file["path"], geometry=vtk_file["geometry"])
    assert isinstance(ds, IdefixVtkDataset)


def test_load_stretched(vtk_file):
    yt_idefix.load_stretched(vtk_file["path"], geometry=vtk_file["geometry"])
