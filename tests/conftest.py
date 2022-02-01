from __future__ import annotations

import os
import re
import sys
from pathlib import Path
from typing import Any

import pytest
import unyt as un
import yaml


def pytest_configure(config):
    if sys.version_info >= (3, 10):
        # from numpy, still visible in 1.22.x
        config.addinivalue_line(
            "filterwarnings",
            (
                "ignore:The distutils.sysconfig module is deprecated, use sysconfig instead:DeprecationWarning"
            ),
        )


def parse_quantity(s) -> un.unyt_quantity:
    # FUTURE: use unyt.unyt_quantity.from_string instead
    # https://github.com/yt-project/unyt/pull/191

    # This is partially adapted from the following SO thread
    # https://stackoverflow.com/questions/41668588/regex-to-match-scientific-notation
    _NUMB_PATTERN = r"^[+/-]?((?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)|\d*\.?\d+|\d+\.?\d*|nan\s|inf\s)"  # noqa: E501
    # *all* greek letters are considered valid unit string elements.
    # This may be an overshoot. We rely on unyt.Unit to do the actual validation
    _UNIT_PATTERN = r"([α-ωΑ-Ωa-zA-Z]+(\*\*([+/-]?[0-9]+)|[*/])?)+"
    _QUAN_PATTERN = rf"{_NUMB_PATTERN}\s*{_UNIT_PATTERN}"
    _NUMB_REGEXP = re.compile(_NUMB_PATTERN)
    _UNIT_REGEXP = re.compile(_UNIT_PATTERN)
    _QUAN_REGEXP = re.compile(_QUAN_PATTERN)

    v = s.strip()
    match = re.fullmatch(_NUMB_PATTERN, v)
    if match is not None:
        return float(match.group()) * un.Unit()
    if not re.match(_QUAN_REGEXP, v):
        raise ValueError(f"Received invalid quantity expression '{s}'.")
    res = re.search(_NUMB_REGEXP, v)
    if res is None:
        raise ValueError
    num = res.group()
    res = re.search(_UNIT_REGEXP, v[res.span()[1] :])
    if res is None:
        raise ValueError
    unit = res.group()
    return float(num) * un.Unit(unit)


DATA_DIR = Path(__file__).parent / "data"

VTK_FILES: dict[str, dict[str, Any]] = {}


def load_meta(pdir, meta_file):
    with open(meta_file) as fh:
        metadata = yaml.load(fh, yaml.SafeLoader)

    metadata["attrs"]["path"] = pdir / metadata["attrs"]["path"]
    if "units" in metadata["attrs"]:
        keys = list(metadata["attrs"]["units"].keys())
        for u in keys:
            metadata["attrs"]["units"][u] = parse_quantity(
                metadata["attrs"]["units"][u]
            )
    return metadata


for ddir in os.listdir(DATA_DIR):
    pdir = DATA_DIR / ddir
    meta_file = pdir / "meta.yaml"
    if not pdir.is_dir():
        continue
    if not meta_file.is_file():
        continue

    metadata = load_meta(pdir, meta_file)
    VTK_FILES.update({metadata["id"]: metadata["attrs"]})


@pytest.fixture(params=VTK_FILES.values(), ids=VTK_FILES.keys(), scope="session")
def vtk_file(request):
    return request.param


# useful subsets
VTK_FILES_NO_GEOMETRY = {
    k: v for k, v in VTK_FILES.items() if v["has_geometry"] is False
}


@pytest.fixture(
    params=VTK_FILES_NO_GEOMETRY.values(),
    ids=VTK_FILES_NO_GEOMETRY.keys(),
    scope="session",
)
def vtk_file_no_geom(request):
    return request.param


VTK_FILES_WITH_GEOMETRY = {
    k: v for k, v in VTK_FILES.items() if v["has_geometry"] is True
}


@pytest.fixture(
    params=VTK_FILES_WITH_GEOMETRY.values(),
    ids=VTK_FILES_WITH_GEOMETRY.keys(),
    scope="session",
)
def vtk_file_with_geom(request):
    return request.param


IDEFIX_VTK_FILES = {k: v for k, v in VTK_FILES.items() if v["kind"] == "idefix"}


@pytest.fixture(
    params=IDEFIX_VTK_FILES.values(), ids=IDEFIX_VTK_FILES.keys(), scope="session"
)
def idefix_vtk_file(request):
    return request.param


PLUTO_VTK_FILES = {k: v for k, v in VTK_FILES.items() if v["kind"] == "pluto"}


@pytest.fixture(
    params=PLUTO_VTK_FILES.values(), ids=PLUTO_VTK_FILES.keys(), scope="session"
)
def pluto_vtk_file(request):
    return request.param
