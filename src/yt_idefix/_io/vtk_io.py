from __future__ import annotations

import struct
import warnings
from typing import Any, BinaryIO, Literal, overload

import numpy as np

from .commons import Coordinates, Shape, get_native_coordinates_from_cartesian

KNOWN_GEOMETRIES: dict[int, str] = {
    0: "cartesian",
    1: "polar",
    2: "spherical",
    3: "cylindrical",
}


def read_header(filename: str) -> str:
    with open(filename, "rb") as fh:
        return "".join(fh.readline(256).decode() for _ in range(2))


@overload
def read_single_field(
    fh: BinaryIO,
    *,
    shape: tuple[int, int, int],
    offset: int | None = None,
    skip_data: Literal[False],
) -> np.ndarray: ...


@overload
def read_single_field(
    fh: BinaryIO,
    *,
    shape: tuple[int, int, int],
    offset: int | None = None,
    skip_data: Literal[True],
) -> None: ...


def read_single_field(
    fh,
    *,
    shape,
    offset=None,
    skip_data=False,
):
    count = np.prod(shape)
    if offset is not None and fh.tell() != offset:
        fh.seek(offset)
    if skip_data:
        fh.seek(count * np.dtype("f").itemsize, 1)
        data = None
    else:
        data = np.fromfile(fh, ">f", count=count)
        data.shape = shape[::-1]
        data = data.T
    return data


def read_shape(s: str) -> Shape:
    # read a specific line containing nx, ny, nz

    assert s.startswith("DIMENSIONS")
    raw: list[str] = s.split()[1:]
    if len(raw) != 3:
        raise RuntimeError
    return Shape(*(int(_) for _ in raw))


def parse_shape(s: str, md: dict[str, Any]) -> None:
    md["shape"] = read_shape(s)


# this may not be kept in the following form
def read_metadata(fh: BinaryIO) -> dict[str, Any]:
    fh.seek(0)
    # skip over the first 4 lines which normally contains
    # VTK DataFile Version x.x
    # <Comments>
    # BINARY
    # DATASET RECTILINEAR_GRID or STRUCTURED_GRID
    for _ in range(4):
        next(fh)

    metadata: dict[str, Any] = {}
    line = next(fh).decode()  # DIMENSIONS NX NY NZ or FIELD
    if line.startswith("FIELD"):
        # Idefix >= 0.8
        nfield = int(line.split()[2])
        for _ in range(nfield):
            d = next(fh).decode()
            if d.startswith("GEOMETRY"):
                geom_flag: int = struct.unpack(">i", fh.read(4))[0]
                geometry_from_data = KNOWN_GEOMETRIES.get(geom_flag)
                if geometry_from_data is None:
                    warnings.warn(
                        f"Unknown geometry enum value {geom_flag}, please report this.",
                        stacklevel=2,
                    )
                metadata["geometry"] = geometry_from_data
            elif d.startswith("TIME"):
                metadata["time"] = struct.unpack(">f", fh.read(4))[0]
            elif d.startswith("PERIODICITY"):
                metadata["periodicity"] = tuple(
                    np.fromfile(fh, dtype=">i4", count=3).astype(bool)
                )
            else:
                warnings.warn(f"Found unknown field {d!r}", stacklevel=2)
            next(fh)  # skip extra linefeed (empty line)
        parse_shape(next(fh).decode(), metadata)

    elif line.startswith("DIMENSIONS"):
        parse_shape(line, metadata)

    else:
        raise RuntimeError(f"Failed to parse {line!r}")

    return metadata


def read_grid_coordinates(
    fh: BinaryIO,
    *,
    geometry: str | None = None,
) -> Coordinates:
    # Return cell edges coordinates

    md = read_metadata(fh)

    geometry = md.get("geometry", geometry)
    if geometry not in (valid_geometries := tuple(KNOWN_GEOMETRIES.values())):
        raise ValueError(
            f"Got unknown geometry {geometry!r}, expected one of {valid_geometries}"
        )

    shape = md["shape"]
    coords: list[np.ndarray] = []
    # now assuming that fh is positioned at the end of metadata
    if geometry in ("cartesian", "cylindrical"):
        # In Idefix, cylindrical geometry is only meant to be used in 2D,
        # so the grid structure is effectively cartesian (R, z)
        for nx in shape:
            next(fh)
            coords.append(np.fromfile(fh, dtype=">f", count=nx))
            next(fh)
        line = next(fh).decode()
        point_type, npoints = (
            t(_) for t, _ in zip((str, int), line.split(), strict=True)
        )
        next(fh)
        if point_type == "CELL_DATA":
            array_shape = shape.to_cell_centered()
        else:
            if point_type == "POINT_DATA":
                # raw data is cell center coords, extrapolation would be needed here
                # However it's probably never going to be worth it since this branch
                # corresponds to Idefix < 0.8 and should never be hit in practice.
                # Raising a warning should suffice
                warnings.warn(
                    "point_type=POINT_DATA case is not fully supported. "
                    "Domain edges will be slightly off.",
                    stacklevel=2,
                )
            else:
                warnings.warn(
                    f"Got unexpected value point_type={point_type!r}. "
                    "Results are not guaranteed.",
                    stacklevel=2,
                )
            array_shape = shape

    else:
        assert geometry in ("polar", "spherical")
        rshape = Shape(*reversed(shape))
        npoints = int(next(fh).decode().split()[1])  # POINTS NXNYNZ float
        assert shape.size == npoints
        points = np.fromfile(fh, dtype=">f", count=3 * npoints)
        next(fh)

        xcart = points[::3]
        xcart.shape = rshape
        xcart = xcart.T

        ycart = points[1::3]
        ycart.shape = rshape
        ycart = ycart.T

        zcart = points[2::3]
        zcart.shape = rshape
        zcart = zcart.T

        coords = get_native_coordinates_from_cartesian(xcart, ycart, zcart, geometry)

        data_type = next(fh).decode().split()[0]  # CELL_DATA (NX-1)(NY-1)(NZ-1)
        next(fh)

        array_shape = shape.to_cell_centered()
        assert data_type == "CELL_DATA"

    def warn_invalid(arr):
        bulk_msg = (
            "This may result from uninitialized data being written to disk, "
            "as is a known bug in PLUTO up to version 4.4.patch2 "
            "in 1D spherical simulations. "
            "If you are seeing this warning in a different situation, "
            "please report this."
        )
        if any(np.isnan(arr)):
            warnings.warn(
                f"NaNs found in coordinate array, behaviour is undefined. {bulk_msg}",
                stacklevel=8,
            )
        elif not np.all(arr[:-1] <= arr[1:]):
            warnings.warn(
                f"Coordinate array is not sorted, behaviour is undefined. {bulk_msg}",
                stacklevel=8,
            )
        return arr

    return Coordinates(
        warn_invalid(coords[0]),
        warn_invalid(coords[1]),
        warn_invalid(coords[2]),
        array_shape,
    )


def read_field_offset_index(
    fh: BinaryIO, shape: Shape, *, default_field_list: list[str]
) -> dict[str, int]:
    # assuming fh is correctly positioned (read_grid_coordinates must be called first)
    retv: dict[str, int] = {}

    while True:
        line = fh.readline()
        if len(line) < 2:
            break
        s = line.decode()
        datatype, varname, dtype = s.split()

        # some versions of Pluto define field names in lower case
        # so we normalize standard output field names to upper case
        # to avoid duplicating data in PlutoFields.known_other_fields
        if varname.upper() in default_field_list:
            varname = varname.upper()

        if datatype == "SCALARS":
            next(fh)
            retv[varname] = fh.tell()
            read_single_field(fh, shape=shape, skip_data=True)
        elif datatype == "VECTORS":
            for axis in "XYZ":
                vname = f"{varname}_{axis}"
                retv[vname] = fh.tell()
                read_single_field(fh, shape=shape, skip_data=True)
        else:
            raise RuntimeError(f"Unknown datatype {datatype!r}")
        fh.readline()
    return retv
