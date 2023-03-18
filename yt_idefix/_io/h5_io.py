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

    nodesX = fh["/node_coords/X"].astype("=f8")
    nodesY = fh["/node_coords/Y"].astype("=f8")
    nodesZ = fh["/node_coords/Z"].astype("=f8")
    dimensions = len(np.array(nodesX).shape)
    shape = Shape(
        *np.array(nodesX).shape
    )  # this is reversed compared the vtk implementation in vtk_io.py
    coords: list[np.ndarray] = []
    # now assuming that fh is positioned at the end of metadata
    if geometry in ("cartesian", "cylindrical") or (
        geometry in ("polar", "spherical") and dimensions == 1
    ):
        if dimensions == 1:
            nodesX = np.array(nodesX)
            # Default is assumed and not parsed from pluto.ini/grid.out (not present in data)
            nodesY = np.array([0.0, 1.0])
            nodesZ = np.array([0.0, 1.0])
            if geometry == "spherical":
                if np.fabs(np.max(nodesX) - np.min(nodesX)) < 1e-8:
                    nodesX = fh["/cell_coords/X"].astype("=f8")
                    nodesX = np.hstack(([nodesX[1] - nodesX[0]], np.array(nodesX)))
                    nodesX = np.array(nodesX) / np.sin(0.5)
            elif geometry == "polar":
                if np.fabs(np.max(nodesX) - np.min(nodesX)) < 1e-8:
                    nodesX = fh["/cell_coords/X"].astype("=f8")
                    nodesX = np.hstack(([nodesX[1] - nodesX[0]], np.array(nodesX)))
                    nodesX = np.array(nodesX) / np.cos(0.5)
            array_shape = Shape(shape[0], 1, 1).to_cell_centered()
        elif dimensions == 2:
            nodesX = np.array(nodesX[0, :])
            nodesY = np.array(nodesY[:, 0])
            # Default is assumed and not parsed from pluto.ini/grid.out (not present in data)
            if geometry == "cartesian":
                nodesZ = np.array([0.0, 1.0])
            else:
                nodesZ = np.array([0.0, 2 * np.pi])
            array_shape = Shape(*reversed(shape[:-1])).to_cell_centered()
        else:
            nodesX = np.array(nodesX)[0, 0, :]
            nodesY = np.array(nodesY)[0, :, 0]
            nodesZ = np.array(nodesZ)[:, 0, 0]
            array_shape = Shape(*reversed(shape)).to_cell_centered()
        coords = [nodesX, nodesY, nodesZ]
    elif geometry in ("polar", "spherical") and dimensions > 1:
        if dimensions == 2:
            nodesX = np.array([nodesX])
            nodesY = np.array([nodesY])
            nodesZ = np.array([nodesZ])
            array_shape = Shape(*reversed(shape[:-1])).to_cell_centered()
        else:
            array_shape = Shape(*reversed(shape)).to_cell_centered()

        xcart = np.transpose(nodesX, axes=(2, 1, 0))
        ycart = np.transpose(nodesY, axes=(2, 1, 0))
        zcart = np.transpose(nodesZ, axes=(2, 1, 0))

        coords = get_native_coordinates_from_cartesian(xcart, ycart, zcart, geometry)
    fh.close()
    return Coordinates(coords[0], coords[1], coords[2], array_shape)
