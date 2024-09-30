from __future__ import annotations

import re
import struct
import sys
import warnings
from enum import IntEnum
from typing import BinaryIO, Literal, cast, overload

import numpy as np

from .commons import Dim, IdefixFieldProperties, IdefixMetadata, Prec

if sys.version_info >= (3, 11):
    from typing import assert_never
else:
    from typing_extensions import assert_never


KNOWN_GEOMETRIES: dict[int, str] = {
    1: "cartesian",
    2: "cylindrical",
    3: "polar",
    4: "spherical",
}


class CharCount(IntEnum):
    # hardcoded in idefix
    HEADER = 128
    NAME = 16


HEADER_REGEXP = re.compile(
    r"Idefix (?P<version>[\w\.-]+) Dump Data" r"((?P<byteorder>(little|big)) endian)?"
)

ByteOrder = Literal["little", "big", "native"]


def parse_byteorder(fh: BinaryIO) -> ByteOrder:
    header = read_header(fh)
    match = HEADER_REGEXP.match(header)
    if match is None:
        warnings.warn(
            f"failed to parse dump header {header!r}",
            stacklevel=2,
            category=RuntimeWarning,
        )
        return "native"

    if (res := match.group("byteorder")) in ("little", "big"):
        return cast(Literal["little", "big"], res)
    else:
        # early versions of Idefix didn't include byteorder in dump headers,
        # fallback to native for backward compatibility
        return "native"


def byteorder2alignment(byteorder: ByteOrder) -> Literal["<", ">", "="]:
    if byteorder == "little":
        return "<"
    elif byteorder == "big":
        return ">"
    elif byteorder == "native":
        return "="
    else:
        assert_never(byteorder)


# emulating C++
# enum DataType {DoubleType, SingleType, IntegerType};
DTYPES: dict[int, Prec] = {0: "d", 1: "f", 2: "i", 3: "?"}
DTYPES_2_NUMPY: dict[Prec, str] = {"d": "f8", "f": "f4", "i": "i4"}


def read_null_terminated_string(fh: BinaryIO, maxsize: int = CharCount.NAME) -> str:
    """Read maxsize bytes, but only parse non-null characters."""
    b = fh.read(maxsize)
    s = b.decode("utf-8", errors="backslashreplace")
    s = s.split("\x00", maxsplit=1)[0]
    return s


def read_next_field_properties(
    fh: BinaryIO,
    *,
    byteorder: ByteOrder,
) -> tuple[str, Prec, Dim, np.ndarray]:
    """Emulate Idefix's OutputDump::ReadNextFieldProperty"""
    field_name = read_null_terminated_string(fh)

    alignment = byteorder2alignment(byteorder)

    fmt = f"{alignment}i"
    int_dtype = struct.unpack(fmt, fh.read(struct.calcsize(fmt)))[0]
    dtype = DTYPES[int_dtype]
    ndim = struct.unpack(fmt, fh.read(struct.calcsize(fmt)))[0]
    if not isinstance(ndim, int):
        raise TypeError(ndim)
    if not (1 <= ndim <= 3):
        raise ValueError(ndim)
    ndim = cast(Dim, ndim)
    fmt = f"{alignment}{ndim}i"
    dim = np.array(struct.unpack(fmt, fh.read(struct.calcsize(fmt))))
    return field_name, dtype, ndim, dim


@overload
def read_chunk(
    fh: BinaryIO,
    ndim: int,
    dim: list[int],
    dtype: str,
    *,
    byteorder: ByteOrder,
    is_scalar: bool,
    skip_data: Literal[True],
) -> None: ...


@overload
def read_chunk(
    fh: BinaryIO,
    ndim: int,
    dim: list[int],
    dtype: str,
    *,
    byteorder: ByteOrder,
    is_scalar: Literal[True],
    skip_data: Literal[False],
) -> float: ...


@overload
def read_chunk(
    fh: BinaryIO,
    ndim: int,
    dim: list[int],
    dtype: str,
    *,
    byteorder: ByteOrder,
    is_scalar: Literal[False],
    skip_data: Literal[False],
) -> np.ndarray: ...


@overload
def read_chunk(
    fh: BinaryIO,
    ndim: int,
    dim: list[int],
    dtype: str,
    *,
    byteorder: ByteOrder,
    is_scalar: bool,
    skip_data: bool,
) -> float | np.ndarray | None: ...


def read_chunk(
    fh,
    ndim,
    dim,
    dtype,
    *,
    byteorder: ByteOrder,
    skip_data=False,
    is_scalar=False,
):
    # NOTE: ret type is only dependent on skip_data...
    # this could be better expressed in the type annotations but it would make
    # more sense to just refactor this function to avoid the boolean trap, so I'll keep wonky
    # type hints for now
    assert ndim == len(dim)
    count = np.prod(dim)
    size = count * np.dtype(dtype).itemsize
    if skip_data:
        fh.seek(size, 1)
        return None

    alignment = byteorder2alignment(byteorder)

    # note: this reversal may not be desirable in general
    if is_scalar:
        fmt = f"{alignment}{count}{dtype}"
        retv = struct.unpack(fmt, fh.read(size))[0]
        return retv
    else:
        data = np.fromfile(fh, alignment + DTYPES_2_NUMPY[dtype], count=count)
        data.shape = dim[::-1]
        return data.astype("=" + DTYPES_2_NUMPY[dtype], copy=False).T


@overload
def read_serial(
    fh: BinaryIO,
    ndim: int,
    dim: np.ndarray,
    dtype: str,
    *,
    byteorder: ByteOrder,
    is_scalar: Literal[True],
) -> float: ...


@overload
def read_serial(
    fh: BinaryIO,
    ndim: int,
    dim: np.ndarray,
    dtype: str,
    *,
    byteorder: ByteOrder,
    is_scalar: Literal[False],
) -> np.ndarray: ...


@overload
def read_serial(
    fh: BinaryIO,
    ndim: int,
    dim: np.ndarray,
    dtype: str,
    *,
    byteorder: ByteOrder,
    is_scalar: bool = False,
) -> float | np.ndarray: ...


def read_serial(fh, ndim, dim, dtype, *, byteorder, is_scalar=False):
    """Emulate Idefix's OutputDump::ReadSerial"""
    assert ndim == 1  # corresponds to an error raised in IDEFIX
    return read_chunk(
        fh,
        ndim=ndim,
        dim=dim,
        dtype=dtype,
        byteorder=byteorder,
        is_scalar=is_scalar,
        skip_data=False,
    )


@overload
def read_distributed(
    fh: BinaryIO,
    dim: np.ndarray,
    *,
    byteorder: ByteOrder,
    dtype: str,
    skip_data: Literal[False],
) -> np.ndarray: ...


@overload
def read_distributed(
    fh: BinaryIO,
    dim: np.ndarray,
    *,
    byteorder: ByteOrder,
    dtype: str,
    skip_data: Literal[True],
) -> None: ...


@overload
def read_distributed(
    fh: BinaryIO,
    dim: np.ndarray,
    *,
    byteorder: ByteOrder,
    dtype: str,
    skip_data: bool = False,
) -> np.ndarray | None: ...


def read_distributed(fh, dim, *, byteorder, dtype, skip_data):
    """Emulate Idefix's OutputDump::ReadDistributed"""
    # note: OutputDump::ReadDistributed only reads doubles
    # because chunks written as integers are small enough
    # that parallelization is counter productive.
    # This a design choice on idefix's size.
    return read_chunk(
        fh,
        ndim=len(dim),
        dim=dim,
        byteorder=byteorder,
        dtype=dtype,
        skip_data=skip_data,
    )


# The following functions are originally designed for yt


def read_header(source: str | BinaryIO, /) -> str:
    if isinstance(source, str):
        with open(source, "rb") as fh:
            return read_header(fh)
    else:
        return read_null_terminated_string(source, maxsize=CharCount.HEADER)


def get_field_offset_index(fh: BinaryIO) -> dict[str, int]:
    """
    Go over a dumpfile, parse bytes offsets associated with each field.
    Returns
    -------
    field_index: mapping (field name -> offset)
    """
    field_index = {}

    fh.seek(0)
    byteorder = parse_byteorder(fh)

    # skip grid properties
    for _ in range(9):
        _field_name, dtype, ndim, dim = read_next_field_properties(
            fh, byteorder=byteorder
        )
        read_serial(fh, ndim, dim, dtype, byteorder=byteorder)

    while True:
        offset = fh.tell()
        field_name, dtype, ndim, dim = read_next_field_properties(
            fh, byteorder=byteorder
        )
        if not re.match("^V[cs]-", field_name):
            break
        field_index[field_name] = offset
        read_distributed(fh, dim, dtype=dtype, byteorder=byteorder, skip_data=True)

    return field_index


def read_single_field(
    fh: BinaryIO, field_offset: int, *, byteorder: ByteOrder
) -> np.ndarray:
    """
    Returns
    -------
    data: 3D np.ndarray with dtype float64
    """
    fh.seek(field_offset)
    field_name, dtype, ndim, dim = read_next_field_properties(fh, byteorder=byteorder)
    data = read_distributed(fh, dim, dtype=dtype, byteorder=byteorder, skip_data=False)
    return data


def read_idefix_dmpfile(
    filename: str, skip_data: bool = False
) -> tuple[IdefixFieldProperties, IdefixMetadata]:
    with open(filename, "rb") as fh:
        return read_idefix_dump_from_buffer(fh, skip_data)


def read_idefix_dump_from_buffer(
    fh: BinaryIO, skip_data: bool = False
) -> tuple[IdefixFieldProperties, IdefixMetadata]:
    fh.seek(0)
    byteorder = parse_byteorder(fh)

    data: float | np.ndarray | None
    fprops: IdefixFieldProperties = {}
    fdata: IdefixMetadata = {}
    for _ in range(9):
        # read grid properties
        # (cell centers, left and right edges in 3D -> 9 arrays)
        field_name, dtype, ndim, dim = read_next_field_properties(
            fh, byteorder=byteorder
        )
        data = read_serial(fh, ndim, dim, dtype, byteorder=byteorder)
        fprops[field_name] = dtype, ndim, dim
        fdata[field_name] = data

    field_name, dtype, ndim, dim = read_next_field_properties(fh, byteorder=byteorder)
    while field_name != "eof":
        # note that this could likely be implemented using a call to
        # `iter` with a sentinel value, to the condition that read_next_field_properties
        # would be splitted into 2 parts (I don't the sentinel pattern works with tuples)
        fprops[field_name] = dtype, ndim, dim
        if field_name.startswith(("Vc-", "Vs-")):
            data = read_distributed(
                fh, dim, dtype=dtype, byteorder=byteorder, skip_data=skip_data
            )
        else:
            is_scalar = ndim == 1 and dim[0] == 1
            is_scalar &= field_name not in ("x1", "x2", "x3")
            data = read_serial(
                fh, ndim, dim, dtype, byteorder=byteorder, is_scalar=is_scalar
            )
        fdata[field_name] = data
        field_name, dtype, ndim, dim = read_next_field_properties(
            fh, byteorder=byteorder
        )

    if isinstance(fdata["geometry"], int):
        fdata["geometry"] = KNOWN_GEOMETRIES[fdata["geometry"]]
    else:
        raise RuntimeError(
            f"Got a geometry flag with unexpected type {type(fdata['geometry'])}; "
            f"got {fdata['geometry']}"
        )

    return fprops, fdata
