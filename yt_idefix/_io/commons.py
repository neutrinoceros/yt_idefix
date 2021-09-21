from __future__ import annotations

import sys
from enum import IntEnum
from typing import Any, Dict, NamedTuple, Tuple

import numpy as np

if sys.version_info >= (3, 8):
    from functools import cached_property
    from typing import Literal

    Prec = Literal["d", "f", "i"]
    Dim = Literal[1, 2, 3]
else:
    cached_property = property
    Prec = str
    Dim = int


class Shape(NamedTuple):
    n1: int = 1
    n2: int = 1
    n3: int = 1

    @cached_property
    def size(self):
        return self.n1 * self.n2 * self.n3

    def to_cell_centered(self) -> Shape:
        vals: list[int] = [1, 1, 1]
        for i, attr in enumerate(self):
            vals[i] = max(1, attr - 1)
        return Shape(*vals)


class Coordinates(NamedTuple):
    x1: np.ndarray
    x2: np.ndarray
    x3: np.ndarray

    @cached_property
    def shape(self) -> Shape:
        return Shape(len(self.x1), len(self.x2), len(self.x3))


# map field name to numpy array init data:
# precision (-> datatype), dimensionality, [nx, ny, nz]
# the np.ndarray is assumed to contain *dim* elements
IdefixFieldProperties = Dict[str, Tuple[Prec, Dim, np.ndarray]]

# Map various str keys to scalars and arrays
IdefixMetadata = Dict[str, Any]


class ByteSize(IntEnum):
    """Constant byte size for data types"""

    CHAR = 1
    INT = 4
    FLOAT = 4
    DOUBLE = 8
