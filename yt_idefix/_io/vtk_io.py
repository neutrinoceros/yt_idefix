from __future__ import annotations

import struct
import warnings
from typing import Any, BinaryIO

import numpy as np

from .commons import ByteSize, Coordinates, Shape

KNOWN_GEOMETRIES = {0: "cartesian", 1: "polar", 2: "spherical"}


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
def read_metadata(fh: BinaryIO) -> dict[str, Any]:

    fh.seek(0)
    # skip over the first 4 lines which normally contains
    # VTK DataFile Version x.x
    # <Comments>
    # BINARY
    # DATASET RECTILINEAR_GRID or STRUCTURED_GRID
    for _ in range(4):
        next(fh)

    metadata = {}
    line = next(fh).decode()  # DIMENSIONS NX NY NZ or FIELD
    if line.startswith("FIELD"):
        # Idefix >= 0.8
        nfield = int(line.split()[2])
        for _ in range(nfield):
            d = next(fh).decode()
            if d.startswith("GEOMETRY"):
                geom_flag: int = struct.unpack(">i", fh.read(4))[0]
                metadata["geometry"] = KNOWN_GEOMETRIES.get(geom_flag, "unknown")
            elif d.startswith("TIME"):
                metadata["time"] = struct.unpack(">f", fh.read(4))[0]
            else:
                warnings.warn(f"Found unknown field {d!r}")
            next(fh)  # skip extra linefeed (empty line)
        parse_shape(next(fh).decode(), metadata)

    elif line.startswith("DIMENSIONS"):
        # Idefix < 0.8
        parse_shape(line, metadata)
        warnings.warn(
            "yt_idefix has limited support for vtk files with idefix < 0.8 . "
            "Assuming a cartesian geometry."
        )
        metadata["geometry"] = "cartesian"

    else:
        raise RuntimeError(f"Failed to parse {line!r}")

    return metadata


def read_grid_coordinates(
    fh: BinaryIO, md: dict[str, Any] | None = None
) -> Coordinates:
    if md is None:
        fh.seek(0)
        md = read_metadata(fh)

    coords: list[np.ndarray] = []
    # now assuming that fh is positioned at the end of metadata
    if md["geometry"] == "cartesian":
        for key in ("n1", "n2", "n3"):
            next(fh)
            coords.append(np.fromfile(fh, dtype=">f", count=getattr(md["shape"], key)))
            next(fh)
        line = next(fh).decode()
        point_type, npoints = (t(_) for t, _ in zip((str, int), line.split()))
        next(fh)
        if point_type == "CELL_DATA":
            # raw data is cell face coords
            # we'd need to interpolate cell centers
            warnings.warn(
                "point_type = CELL_DATA is not fully supported yet. "
                "This will be treated as POINT_DATA instead"
            )
            md["array_shape"] = md["shape"].to_cell_centered()
        else:
            md["array_shape"] = md["shape"]

    else:
        raise NotImplementedError("Non cartesian coordinates are not implemented yet")

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
