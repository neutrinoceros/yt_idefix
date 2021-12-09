from __future__ import annotations

import struct
import warnings
from typing import Any, BinaryIO

import numpy as np

from .commons import ByteSize, Coordinates, Shape

KNOWN_GEOMETRIES = {0: "cartesian", 1: "polar", 2: "spherical", 3: "cylindrical"}


def read_header(filename: str) -> str:
    with open(filename, "rb") as fh:
        return "".join(fh.readline(256).decode() for _ in range(2))


def read_single_field(
    fh: BinaryIO,
    offset: int | None = None,
    *,
    shape: tuple[int, int, int],
    skip_data: bool = False,
) -> np.ndarray | None:
    count = np.prod(shape)
    if offset is not None and fh.tell() != offset:
        fh.seek(offset)
    if skip_data:
        fh.seek(count * ByteSize.FLOAT, 1)
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
def read_metadata(fh: BinaryIO, *, geometry: str | None = None) -> dict[str, Any]:

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
                        f"Unknown geometry enum value {geom_flag}, please report this."
                    )
                elif geometry not in (None, geometry_from_data):
                    warnings.warn(
                        f"Got inconsistent geometries:\n"
                        f" - {geometry_from_data!r} (from file)\n"
                        f" - {geometry!r} (from user)\n"
                        "Ignoring user input"
                    )
                metadata["geometry"] = geometry_from_data
            elif d.startswith("TIME"):
                metadata["time"] = struct.unpack(">f", fh.read(4))[0]
            elif d.startswith("PERIODICITY"):
                metadata["periodicity"] = tuple(
                    np.fromfile(fh, dtype=">i4", count=3).astype(bool)
                )
            else:
                warnings.warn(f"Found unknown field {d!r}")
            next(fh)  # skip extra linefeed (empty line)
        parse_shape(next(fh).decode(), metadata)

    elif line.startswith("DIMENSIONS"):
        parse_shape(line, metadata)

    else:
        raise RuntimeError(f"Failed to parse {line!r}")

    if "geometry" not in metadata:
        # Idefix < 0.8, or PLUTO
        if geometry is None:
            raise ValueError(
                "Geometry couldn't be parsed from file. "
                "The 'geometry' keyword argument must be specified."
            )
        metadata["geometry"] = geometry

    return metadata


def read_grid_coordinates(
    fh: BinaryIO, md: dict[str, Any] | None = None
) -> Coordinates:
    # Return cell edges coordinates

    if md is None:
        fh.seek(0)
        md = read_metadata(fh)

    geometry = md["geometry"]
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
        point_type, npoints = (t(_) for t, _ in zip((str, int), line.split()))
        next(fh)
        if point_type == "CELL_DATA":
            md["array_shape"] = shape.to_cell_centered()
        else:
            if point_type == "POINT_DATA":
                # raw data is cell center coords, extrapolation would be needed here
                # However it's probably never going to be worth it since this branch
                # corresponds to Idefix < 0.8 and should never be hit in practice.
                # Raising a warning should suffice
                warnings.warn(
                    "point_type=POINT_DATA case is not fully supported. "
                    "Domain edges will be slightly off."
                )
            else:
                warnings.warn(
                    f"Got unexpected value point_type={point_type!r}. "
                    "Results are not guaranteed."
                )
            md["array_shape"] = shape

    elif geometry in ("polar", "spherical"):
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

        # Reconstruct the polar coordinate system
        if geometry == "polar":
            # on disk coordinates are cell face coordinates.
            # We need to convert from cartesian to polar,
            # but no interpolation is needed.
            r = np.sqrt(xcart[:, 0, 0] ** 2 + ycart[:, 0, 0] ** 2)
            theta = np.unwrap(np.arctan2(ycart[0, :, 0], xcart[0, :, 0]))
            z = zcart[0, 0, :]

            data_type = next(fh).decode().split()[0]  # CELL_DATA (NX-1)(NY-1)(NZ-1)
            assert data_type == "CELL_DATA"
            next(fh)

            # manually changing phase origin (theta) to match
            # results from Idefix's pytools
            coords = [r, theta + np.pi, z]
        elif geometry == "spherical":
            # Reconstruct the spherical coordinate system
            if shape.n3 == 1:
                r = np.sqrt(xcart[:, 0, 0] ** 2 + ycart[:, 0, 0] ** 2)
                phi = np.unwrap(
                    np.arctan2(zcart[0, shape.n2 // 2, :], xcart[0, shape.n2 // 2, :])
                )
                theta = np.arccos(
                    ycart[0, :, 0] / np.sqrt(xcart[0, :, 0] ** 2 + ycart[0, :, 0] ** 2)
                )
            else:
                r = np.sqrt(
                    xcart[:, 0, 0] ** 2 + ycart[:, 0, 0] ** 2 + zcart[:, 0, 0] ** 2
                )
                phi = np.unwrap(
                    np.arctan2(
                        ycart[shape.n1 // 2, shape.n2 // 2, :],
                        xcart[shape.n1 // 2, shape.n2 // 2, :],
                    )
                )
                theta = np.arccos(
                    zcart[0, :, 0]
                    / np.sqrt(
                        xcart[0, :, 0] ** 2 + ycart[0, :, 0] ** 2 + zcart[0, :, 0] ** 2
                    )
                )
            coords = [r, theta, phi]

            data_type = next(fh).decode().split()[0]  # CELL_DATA (NX-1)(NY-1)(NZ-1)
            if data_type != "CELL_DATA":
                raise RuntimeError
            next(fh)
        else:
            raise RuntimeError("This should be logically impossible.")

        md["array_shape"] = shape.to_cell_centered()

    else:
        raise RuntimeError(f"Found unknown geometry {geometry!r}")

    return Coordinates(*coords)


def read_field_offset_index(fh: BinaryIO, shape: Shape) -> dict[str, int]:
    # assuming fh is correctly positioned (read_grid_coordinates must be called first)
    retv: dict[str, int] = {}

    while True:
        line = fh.readline()
        if len(line) < 2:
            break
        s = line.decode()
        datatype, varname, dtype = s.split()

        # some versions of Pluto define field names in lower case
        # so we normalize to upper case to avoid duplicating data
        # in IdefixVtkFieldInfo.known_other_fields
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
