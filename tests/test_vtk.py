import os
import re

import pytest
from more_itertools import distinct_combinations
from unyt import Unit, assert_allclose_units

import yt
from yt_idefix.api import IdefixVtkDataset, PlutoVtkDataset
from yt_idefix.data_structures import VtkMixin

# A sample list of units for test.
# The first three values are chosen randomly
# and others are calculated correspondingly.
SAMPLE_UNITS = {
    "time_unit": (2.0, "s"),
    "length_unit": (4.0, "cm"),
    "mass_unit": (5.0, "kg"),
    "density_unit": (0.078125, "kg/cm**3"),
    "velocity_unit": (2.0, "cm/s"),
    "magnetic_unit": (62.66570686577499, "gauss"),
}


def test_class_validation(vtk_file):
    file = vtk_file
    cls = {
        "idefix": IdefixVtkDataset,
        "pluto": PlutoVtkDataset,
    }[file["kind"]]
    assert cls._is_valid(str(file["path"]))


def test_parse_pluto_metadata(pluto_vtk_file):
    file = pluto_vtk_file
    ds = yt.load(file["path"])
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


def test_units(vtk_file_with_units):
    file = vtk_file_with_units
    ds = yt.load(file["path"])
    for u, expected in file["units"].items():
        assert_allclose_units(getattr(ds, f"{u}_unit"), expected)


def test_pluto_complete_units_override(pluto_vtk_file):
    su = SAMPLE_UNITS
    unit_combs = distinct_combinations(su, 3)
    for comb in unit_combs:
        uo = {}
        for unit in comb:
            uo[unit] = su[unit]
        if set(uo) in PlutoVtkDataset.invalid_unit_combinations:
            continue
        ds = yt.load(
            pluto_vtk_file["path"],
            geometry=pluto_vtk_file["geometry"],
            units_override=uo,
        )
        for attr, value in su.items():
            assert_allclose_units(getattr(ds, attr), ds.quan(*value))


def test_pluto_one_unit_override(vtk_file_with_units):
    file = vtk_file_with_units
    # Pluto's length_unit and density_unit will be combined with
    uo = {"time_unit": (2.0, "yr")}
    ds = yt.load(file["path"], geometry=file["geometry"], units_override=uo)
    # Pluto's length_unit should be preserved
    assert_allclose_units(ds.length_unit, file["units"]["length"])
    # Pluto's velocity_unit should be overrided
    expect_velocity = ds.length_unit / ds.quan(*uo["time_unit"])
    assert_allclose_units(ds.velocity_unit, expect_velocity)


def test_pluto_two_units_override(vtk_file_with_units):
    file = vtk_file_with_units
    # Pluto's length_unit will be combined with
    uo = {"time_unit": (2.0, "yr"), "density_unit": (32.0, "g/cm**3")}
    ds = yt.load(file["path"], geometry=file["geometry"], units_override=uo)
    # Pluto's length_unit should be preserved
    assert_allclose_units(ds.length_unit, file["units"]["length"])
    # Pluto's velocity_unit should be overrided
    expect_velocity = ds.length_unit / ds.quan(*uo["time_unit"])
    expect_mass = ds.quan(*uo["density_unit"]) * ds.length_unit**3
    assert_allclose_units(ds.velocity_unit, expect_velocity)
    assert_allclose_units(ds.mass_unit, expect_mass)


def test_pluto_invalid_units_override(pluto_vtk_file):
    for uo in PlutoVtkDataset.invalid_unit_combinations:
        with pytest.raises(ValueError, match=r".* cannot derive all units\n.*"):
            yt.load(
                pluto_vtk_file["path"],
                geometry=pluto_vtk_file["geometry"],
                units_override=uo,
            )


def test_pluto_temperature_unit_override(pluto_vtk_file):
    with pytest.raises(
        ValueError,
        match=(
            "Temperature is not allowed in units_override, "
            "since it's always in Kelvin in PLUTO"
        ),
    ):
        yt.load(
            pluto_vtk_file["path"],
            geometry=pluto_vtk_file["geometry"],
            units_override={"temperature_unit": (2.0, "K")},
        )


def test_pluto_over_units_override(pluto_vtk_file):
    with pytest.raises(
        ValueError,
        match=(
            re.escape(
                "More than 3 degrees of freedom were specified "
                "in units_override (6 given)"
            )
        ),
    ):
        yt.load(
            pluto_vtk_file["path"],
            geometry=pluto_vtk_file["geometry"],
            units_override=SAMPLE_UNITS,
        )


def test_pluto_wrong_definitions_header(pluto_vtk_file):
    with pytest.raises(
        FileNotFoundError,
        match=(r".*No such file or directory: '.*definitions2\.h'"),
    ):
        yt.load(pluto_vtk_file["path"], definitions_header="definitions2.h")


def test_code_time(vtk_file_with_units):
    ds = yt.load(vtk_file_with_units["path"])
    code_time = Unit("code_time", registry=ds.unit_registry)
    assert_allclose_units(ds.time_unit, 1.0 * code_time)
    assert_allclose_units(ds.current_time.in_cgs(), ds.current_time.value * code_time)
    quantity = ds.quan(1.0, "code_time")
    assert_allclose_units(quantity.in_cgs(), quantity.value * code_time)


def test_data_access(vtk_file):
    file = vtk_file
    ds = yt.load(file["path"], geometry=file["geometry"])
    ad = ds.all_data()
    ad["gas", "density"]


def test_load_user_passed_geometry(vtk_file):
    file = vtk_file
    ds = yt.load(file["path"], geometry=file["geometry"])
    assert ds.geometry == file["geometry"]
    assert ds.dimensionality == file["dimensionality"]


def test_load_missing_geometry_metadata(vtk_file_no_geom):
    with pytest.raises(
        ValueError,
        match=(
            "Geometry couldn't be parsed from disk. "
            "The 'geometry' keyword argument must be specified."
        ),
    ):
        yt.load(vtk_file_no_geom["path"])


def test_load_parseable_geometry_metadata(vtk_file_with_geom):
    file = vtk_file_with_geom
    ds = yt.load(file["path"])
    assert ds.geometry == file["geometry"]


def test_derived_field(vtk_file):
    file = vtk_file
    ds = yt.load(file["path"], geometry=file["geometry"])
    # this derived field is compatible for all current test data
    # it'll be replaced once a better universal field is defined internally.

    def mom_den(field, data):
        return data["RHO"] * data["VX1"]

    ds.add_field(
        name=("gas", "momentum_density"),
        function=mom_den,
        units="g*cm**-2/s",
        sampling_type="cell",
    )
    ad = ds.all_data()
    test = ad["gas", "momentum_density"]
    expect = ad["RHO"] * ad["VX1"]
    assert_allclose_units(test, expect)


# TODO: make this a pytest-mpl test
def test_slice_plot(vtk_file):
    file = vtk_file
    ds = yt.load(file["path"], geometry=file["geometry"], unit_system="code")
    yt.SlicePlot(ds, normal=(0, 0, 1), fields=("gas", "density"))


def test_projection_plot(vtk_file):
    file = vtk_file
    ds = yt.load(file["path"], geometry=file["geometry"], unit_system="code")
    yt.ProjectionPlot(ds, normal=(0, 0, 1), fields=("gas", "density"))


def test_load_magic(vtk_file):
    ds = yt.load(vtk_file["path"], geometry=vtk_file["geometry"])
    assert isinstance(ds, VtkMixin)
