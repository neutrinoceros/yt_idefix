from pathlib import Path

import numpy as np
from numpy.testing import assert_equal

import yt
import yt_idefix  # noqa: F401
from yt_idefix._io.dmp_io import read_idefix_dmpfile
from yt_idefix.api import IdefixDmpDataset

DATA_DIR = Path(__file__).parent / "data"
idefix_khi = DATA_DIR.joinpath("khi", "dump.0100.dmp")


def test_load():
    ds = yt.load(idefix_khi)
    assert isinstance(ds, IdefixDmpDataset)
    assert ds.dimensionality == 2


def test_region():
    ds = yt.load(idefix_khi)
    ds.r[:]


def test_fields():
    ds = yt.load(idefix_khi)
    expected = [
        ("idefix-dmp", "Vc-RHO"),
        ("idefix-dmp", "Vc-VX1"),
        ("idefix-dmp", "Vc-VX2"),
    ]
    assert ds.field_list == expected

    assert ("gas", "density") in ds.derived_field_list
    assert ("gas", "velocity_x") in ds.derived_field_list
    assert ("gas", "velocity_y") in ds.derived_field_list


def test_get_data():
    ds = yt.load(idefix_khi)
    ds.r[:]["density"]


def test_read_dmp():
    fn = idefix_khi
    fprops, fdata = read_idefix_dmpfile(fn)

    expected_fields = {
        "x1",
        "xl1",
        "xr1",
        "x2",
        "xl2",
        "xr2",
        "x3",
        "xl3",
        "xr3",
        "Vc-RHO",
        "Vc-VX1",
        "Vc-VX2",
        "time",
        "dt",
        "vtkFileNumber",
        "dumpFileNumber",
        "geometry",
        "periodicity",
    }
    detected_fields = set(fprops.keys())
    assert expected_fields.issubset(detected_fields)

    expected_shape = [1024, 256, 1]
    detected_shape = np.concatenate([fprops[k][-1] for k in ("x1", "x2", "x3")])
    assert_equal(detected_shape, expected_shape)
    for field_name, data in fprops.items():
        if not field_name.startswith("V"):
            continue
        _dtype, _ndim, dim = data
        assert_equal(dim, expected_shape)
