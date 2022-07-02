from __future__ import annotations

import os
import warnings

import numpy as np
from packaging.version import Version

import yt
from yt.utilities.exceptions import YTUnidentifiedDataType
from yt_idefix._io import vtk_io
from yt_idefix.data_structures import IdefixVtkDataset, PlutoVtkDataset

__all__ = ["load", "load_stretched"]


YT_VERSION = Version(yt.__version__)


class VisibleDeprecationWarning(Warning):
    pass


class IdefixVtkDatasetSeries(yt.DatasetSeries):
    _dataset_cls = IdefixVtkDataset


def load(fn, *args, **kwargs):
    msg = (
        "yt_idefix.load is a diminished wrapper for to yt.load, "
        "and is not necessary with yt > 4.0.2\nplease "
    )

    if YT_VERSION <= Version("4.0.2"):
        msg += "upgrade yt and "
    msg += "use yt.load instead"
    msg += (
        "yt_idefix.load was deprecated in version 0.11.0 "
        "and will be removed completely in a future version"
    )
    warnings.warn(
        msg,
        category=VisibleDeprecationWarning,
        stacklevel=2,
    )

    fn = os.path.expanduser(fn)
    if not os.path.exists(fn):
        raise FileNotFoundError(fn)

    if YT_VERSION < Version("4.0.2") and fn.endswith(".vtk"):
        if any(wildcard in fn for wildcard in "[]?!*"):
            return IdefixVtkDatasetSeries(fn, *args, **kwargs)
        elif IdefixVtkDataset._is_valid(fn, *args, **kwargs):
            return IdefixVtkDataset(fn, *args, **kwargs)
        elif PlutoVtkDataset._is_valid(fn, *args, **kwargs):
            return PlutoVtkDataset(fn, *args, **kwargs)
        else:
            raise YTUnidentifiedDataType(fn, *args, **kwargs)
    else:
        return yt.load(fn, *args, **kwargs)


def load_stretched(fn, *, geometry: str | None = None, **kwargs):
    """
    A small, specialized wrapper around yt.load_hexahedral_mesh
    Compatible with streched grids but comes at a significant cost as of yt 4.0
    - limited support for common operations such as projections
    - fields are not properly dimensioned
    - no lazy-loading (all data has to reside in memory)
    - only supports vtk outputs (not dumps)
    """
    warnings.warn(
        "yt_idefix.load_strecthed is deprecated "
        "and will be removed completely in a future version. "
        "Please use yt.load instead.",
        category=VisibleDeprecationWarning,
        stacklevel=2,
    )

    # brute force validation
    if not (
        IdefixVtkDataset._is_valid(fn, geometry=geometry)
        or PlutoVtkDataset._is_valid(fn, geometry=geometry)
    ):
        raise TypeError(
            "yt_idefix.load_stretched only supports Idefix and Pluto vtk files"
        )

    # actual parsing
    with open(fn, "rb") as fh:
        md = vtk_io.read_metadata(fh)
        coords = vtk_io.read_grid_coordinates(
            fh, geometry=geometry or md.get("geometry")
        ).padded()
        field_offset_index = vtk_io.read_field_offset_index(
            fh, shape=coords.array_shape
        )

    if geometry is None:
        geometry = md["geometry"]

    data: dict[str, np.ndarray] = {}
    with open(fn, "rb") as fh:
        for name, offset in field_offset_index.items():
            data[name] = vtk_io.read_single_field(
                fh, shape=coords.array_shape, offset=offset, skip_data=False
            )

    coordinates, connectivity = yt.hexahedral_connectivity(
        coords.x1, coords.x2, coords.x3
    )

    bbox = np.array(
        [
            [coords.x1[0], coords.x1[-1]],
            [coords.x2[0], coords.x2[-1]],
            [coords.x3[0], coords.x3[-1]],
        ]
    )

    ds = yt.load_hexahedral_mesh(
        data,
        coordinates=coordinates,
        connectivity=connectivity,
        bbox=bbox,
        geometry=geometry,
        **kwargs,
    )
    try:
        ds.current_time = md["time"] * ds.time_unit
    except KeyError:
        # old Idefix vtk and Pluto vtks don't contain a "time" field
        pass
    return ds
