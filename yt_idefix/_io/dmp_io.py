from __future__ import annotations

import re
import struct
import warnings
from enum import IntEnum
from typing import BinaryIO

import numpy as np

from .commons import ByteSize, IdefixFieldProperties, IdefixMetadata

KNOWN_GEOMETRIES: dict[int, str] = {
    1: "cartesian",
    2: "cylindrial",
    3: "polar",
    4: "spherical",
}


class CharCount(IntEnum):
    # hardcoded in idefix
    HEADER = 128
    NAME = 16


# emulating C++
# enum DataType {DoubleType, SingleType, IntegerType};
DTYPES = {0: "d", 1: "f", 2: "i"}
DTYPES_2_NUMPY = {"d": "=f8", "f": "=f4", "i": "=i4"}
DTYPES_2_SIZE = {"i": ByteSize.INT, "f": ByteSize.FLOAT, "d": ByteSize.DOUBLE}


def read_null_terminated_string(fh: BinaryIO, maxsize: int = CharCount.NAME):
    """Read maxsize * ByteSize.CHAR bytes, but only parse non-null characters."""
    b = fh.read(maxsize * ByteSize.CHAR)
    s = b.decode("utf-8", errors="backslashreplace")
    s = s.split("\x00", maxsplit=1)[0]
    return s


def read_next_field_properties(fh: BinaryIO):
    """Emulate Idefix's OutputDump::ReadNextFieldProperty"""
    field_name = read_null_terminated_string(fh)

    fmt = "=i"
    dtype = DTYPES[struct.unpack(fmt, fh.read(struct.calcsize(fmt)))[0]]
    ndim = struct.unpack(fmt, fh.read(struct.calcsize(fmt)))[0]
    if ndim > 3:
        raise ValueError(ndim)
    fmt = f"={ndim}i"
    dim = np.array(struct.unpack(fmt, fh.read(struct.calcsize(fmt))))
    return field_name, dtype, ndim, dim


def read_chunk(
    fh: BinaryIO,
    ndim: int,
    dim: list[int],
    dtype: str,
    *,
    is_scalar: bool = False,
    skip_data: bool = False,
) -> np.ndarray | None:
    # NOTE: ret type is only dependent on skip_data...
    # this could be better expressed in the type annotationation but it would make
    # more sense to just refactor this function to avoid the boolean trap, so I'll keep wonky
    # type hints for now
    assert ndim == len(dim)
    count = np.product(dim)
    size = count * DTYPES_2_SIZE[dtype]
    if skip_data:
        fh.seek(size, 1)
        return None

    # note: this reversal may not be desirable in general
    if is_scalar:
        fmt = f"={count}{dtype}"
        retv = struct.unpack(fmt, fh.read(size))[0]
        return retv
    else:
        data = np.fromfile(fh, DTYPES_2_NUMPY[dtype], count=count)
        data.shape = dim[::-1]
        return data.T


def read_serial(
    fh: BinaryIO, ndim: int, dim: list[int], dtype: str, *, is_scalar: bool = False
) -> np.ndarray | None:
    """Emulate Idefix's OutputDump::ReadSerial"""
    assert ndim == 1  # corresponds to an error raised in IDEFIX
    return read_chunk(fh, ndim=ndim, dim=dim, dtype=dtype, is_scalar=is_scalar)


def read_distributed(
    fh: BinaryIO, dim: list[int], *, skip_data: bool = False
) -> np.ndarray | None:
    """Emulate Idefix's OutputDump::ReadDistributed"""
    # note: OutputDump::ReadDistributed only read doubles
    # because chucks written in integers are small enough
    # that parallelization is counter productive.
    # This a design choice on idefix's size.
    return read_chunk(fh, ndim=len(dim), dim=dim, dtype="d", skip_data=skip_data)


# The following functions are originally designed for yt


def read_header(filename: str) -> str:
    with open(filename, "rb") as fh:
        header = read_null_terminated_string(fh, maxsize=CharCount.HEADER)
    return header


def get_field_offset_index(fh: BinaryIO) -> dict[str, int]:
    """
    Go over a dumpfile, parse bytes offsets associated with each field.
    Returns
    -------
    field_index: mapping (field name -> offset)
    """
    field_index = {}

    # skip header
    fh.seek(CharCount.HEADER * ByteSize.CHAR)
    # skip grid properties
    for _ in range(9):
        _field_name, dtype, ndim, dim = read_next_field_properties(fh)
        read_serial(fh, ndim, dim, dtype)

    while True:
        offset = fh.tell()
        field_name, dtype, ndim, dim = read_next_field_properties(fh)
        if not re.match("^V[cs]-", field_name):
            break
        field_index[field_name] = offset
        read_distributed(fh, dim, skip_data=True)

    return field_index


def read_single_field(fh: BinaryIO, field_offset: int) -> np.ndarray:
    """
    Returns
    -------
    data: 3D np.ndarray with dtype float64
    """
    fh.seek(field_offset)
    field_name, dtype, ndim, dim = read_next_field_properties(fh)
    data = read_distributed(fh, dim)
    return data


def read_idefix_dmpfile(
    filename: str, skip_data: bool = False
) -> tuple[IdefixFieldProperties, IdefixMetadata]:
    with open(filename, "rb") as fh:
        return read_idefix_dump_from_buffer(fh, skip_data)


def read_idefix_dump_from_buffer(
    fh: BinaryIO, skip_data: bool = False
) -> tuple[IdefixFieldProperties, IdefixMetadata]:

    # skip header
    fh.seek(CharCount.HEADER * ByteSize.CHAR)

    fprops = {}
    fdata = {}
    for _ in range(9):
        # read grid properties
        # (cell centers, left and right edges in 3D -> 9 arrays)
        field_name, dtype, ndim, dim = read_next_field_properties(fh)
        data = read_serial(fh, ndim, dim, dtype)
        fprops[field_name] = dtype, ndim, dim
        fdata[field_name] = data

    field_name, dtype, ndim, dim = read_next_field_properties(fh)
    while field_name != "eof":
        # note that this could likely be implemented using a call to
        # `iter` with a sentinel value, to the condition that read_next_field_properties
        # would be splitted into 2 parts (I don't the sentinel pattern works with tuples)
        fprops[field_name] = dtype, ndim, dim
        if field_name.startswith("Vc-") or field_name.startswith("Vs-"):
            data = read_distributed(fh, dim, skip_data=skip_data)
        else:
            is_scalar = ndim == 1 and dim[0] == 1
            is_scalar &= field_name not in ("x1", "x2", "x3")
            data = read_serial(fh, ndim, dim, dtype, is_scalar=is_scalar)
        fdata[field_name] = data
        field_name, dtype, ndim, dim = read_next_field_properties(fh)

    if "geometry" not in fdata:
        warnings.warn("Could not determine geometry, assuming cartesian.")
        fdata["geometry"] = "cartesian"
    elif isinstance(fdata["geometry"], int):
        fdata["geometry"] = KNOWN_GEOMETRIES[fdata["geometry"]]
    else:
        raise RuntimeError(
            f"Got a geometry flag with unexpected type {type(fdata['geometry'])}; "
            f"got {fdata['geometry']}"
        )

    return fprops, fdata
