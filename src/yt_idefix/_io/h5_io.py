from __future__ import annotations

import os

import numpy as np

from yt.utilities.on_demand_imports import _h5py as h5py

from .commons import Coordinates, Shape, get_native_coordinates_from_cartesian

KNOWN_GEOMETRIES: dict[int, str] = {
    0: "cartesian",
    1: "polar",
    2: "spherical",
    3: "cylindrical",
}


def read_grid_coordinates(
    filename: str | os.PathLike[str],
    *,
    geometry: str | None = None,
) -> Coordinates:
    # Return cell edges coordinates
    fh = h5py.File(filename, "r")
    if geometry not in (valid_geometries := tuple(KNOWN_GEOMETRIES.values())):
        raise ValueError(
            f"Got unknown geometry {geometry!r}, expected one of {valid_geometries}"
        )

    nodesX = np.array(fh["/node_coords/X"]).astype("=f8", copy=False)
    nodesY = np.array(fh["/node_coords/Y"]).astype("=f8", copy=False)
    nodesZ = np.array(fh["/node_coords/Z"]).astype("=f8", copy=False)

    # this is reversed compared the vtk implementation in vtk_io.py
    shape = Shape(*(nodesX.shape))
    coords: list[np.ndarray] = []
    # now assuming that fh is positioned at the end of metadata
    if geometry in ("cartesian", "cylindrical") or (
        geometry in ("polar", "spherical") and nodesX.ndim == 1
    ):
        if nodesX.ndim == 1:
            # Default is assumed and not parsed from pluto.ini/grid.out (not present in data)
            nodesY = np.array([0.0, 1.0])
            nodesZ = np.array([0.0, 1.0])
            if geometry == "spherical":
                if np.fabs(np.ptp(nodesX)) < 1e-8:
                    nodesX = np.hstack(([nodesX[1] - nodesX[0]], nodesX))
                    nodesX /= np.sin(0.5)
            elif geometry == "polar":
                if np.fabs(np.ptp(nodesX)) < 1e-8:
                    nodesX = np.hstack(([nodesX[1] - nodesX[0]], nodesX))
                    nodesX /= np.cos(0.5)
            array_shape = Shape(shape[0], 1, 1).to_cell_centered()
        elif nodesX.ndim == 2:
            nodesX = nodesX[0, :]
            nodesY = nodesY[:, 0]
            # Default is assumed and not parsed from pluto.ini/grid.out (not present in data)
            if geometry == "cartesian":
                nodesZ = np.array([0.0, 1.0])
            else:
                nodesZ = np.array([0.0, 2 * np.pi])
            array_shape = Shape(*reversed(shape[:-1])).to_cell_centered()
        else:
            nodesX = nodesX[0, 0, :]
            nodesY = nodesY[0, :, 0]
            nodesZ = nodesZ[:, 0, 0]
            array_shape = Shape(*reversed(shape)).to_cell_centered()
        coords = [nodesX, nodesY, nodesZ]
    elif geometry in ("polar", "spherical") and nodesX.ndim > 1:
        if nodesX.ndim == 2:
            nodesX = np.expand_dims(nodesX, axis=0)
            nodesY = np.expand_dims(nodesY, axis=0)
            nodesZ = np.expand_dims(nodesZ, axis=0)
            array_shape = Shape(*reversed(shape[:-1])).to_cell_centered()
        else:
            array_shape = Shape(*reversed(shape)).to_cell_centered()

        xcart = np.transpose(nodesX, axes=(2, 1, 0))
        ycart = np.transpose(nodesY, axes=(2, 1, 0))
        zcart = np.transpose(nodesZ, axes=(2, 1, 0))

        coords = get_native_coordinates_from_cartesian(xcart, ycart, zcart, geometry)
    fh.close()
    return Coordinates(coords[0], coords[1], coords[2], array_shape)
