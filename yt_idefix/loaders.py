from __future__ import annotations

import warnings

import numpy as np

import yt
from yt_idefix._io import vtk_io
from yt_idefix.data_structures import IdefixVtkDataset, PlutoVtkDataset

__all__ = ["load", "load_stretched"]


class VisibleDeprecationWarning(Warning):
    pass


def load(fn, *args, **kwargs):
    warnings.warn(
        "yt_idefix.load is now a strict alias to yt.load, "
        "please use yt.load directly. "
        "yt_idefix.load was deprecated in version 0.11.0 "
        "and will be removed no later than version 0.13.0",
        category=VisibleDeprecationWarning,
        stacklevel=2,
    )
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
        md = vtk_io.read_metadata(fh, geometry=geometry)
        coords = vtk_io.read_grid_coordinates(fh, md).padded()
        shape = md["array_shape"]
        field_offset_index = vtk_io.read_field_offset_index(fh, shape=shape)

    if geometry is None:
        geometry = md["geometry"]

    data: dict[str, np.ndarray] = {}
    with open(fn, "rb") as fh:
        for name, offset in field_offset_index.items():
            data[name] = vtk_io.read_single_field(
                fh, shape=shape, offset=offset, skip_data=False
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

    return yt.load_hexahedral_mesh(
        data,
        coordinates=coordinates,
        connectivity=connectivity,
        bbox=bbox,
        geometry=geometry,
        **kwargs,
    )
