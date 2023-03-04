from __future__ import annotations

import os
import sys
from importlib.metadata import version
from pathlib import Path
from typing import Any

import pytest
import unyt as un
import yaml
from packaging.version import Version

NUMPY_VERSION = Version(version("numpy"))


def pytest_configure(config):
    if sys.version_info >= (3, 10) and NUMPY_VERSION < Version("1.23"):
        config.addinivalue_line(
            "filterwarnings",
            (
                "ignore:The distutils.sysconfig module is deprecated, use sysconfig instead:DeprecationWarning"
            ),
        )


DATA_DIR = Path(__file__).parent / "data"

VTK_FILES: dict[str, dict[str, Any]] = {}


def load_meta(pdir: Path, meta_file: Path) -> list[dict[str:Any]]:
    with open(meta_file) as fh:
        metadata = yaml.load(fh, yaml.SafeLoader)

    if "units" in metadata["attrs"]:
        keys = list(metadata["attrs"]["units"].keys())
        for u in keys:
            metadata["attrs"]["units"][u] = un.unyt_quantity.from_string(
                metadata["attrs"]["units"][u]
            )

    attrs = metadata["attrs"]
    paths = metadata["paths"]
    retv = []
    for p in paths:
        pp = pdir / p
        retv.append(
            {
                "id": f"{metadata['id']}{pp.suffix.replace('.', '_')}",
                "attrs": {**attrs, **{"path": pp}},
            }
        )
    return retv


for ddir in os.listdir(DATA_DIR):
    if not (pdir := DATA_DIR / ddir).is_dir():
        continue

    if not (meta_file := pdir / "meta.yaml").is_file():
        continue

    datasets = load_meta(pdir, meta_file)
    for ds in datasets:
        if ds["attrs"]["path"].suffix == ".vtk":
            VTK_FILES.update({ds["id"]: ds["attrs"]})
        else:
            raise ValueError(f"Failed to determine data type for {ds['path']}")


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


VTK_FILES_WITH_UNITS = {k: v for k, v in VTK_FILES.items() if v["has_units"] is True}


@pytest.fixture(
    params=VTK_FILES_WITH_UNITS.values(),
    ids=VTK_FILES_WITH_UNITS.keys(),
    scope="session",
)
def vtk_file_with_units(request):
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
