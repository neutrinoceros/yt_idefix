from __future__ import annotations

from typing import Any, Literal, NamedTuple

import numpy as np

Prec = Literal["d", "f", "i", "?"]
Dim = Literal[1, 2, 3]


class Shape(NamedTuple):
    n1: int = 1
    n2: int = 1
    n3: int = 1

    @property
    def size(self):
        return self.n1 * self.n2 * self.n3

    def to_cell_centered(self) -> Shape:
        vals: list[int] = [1, 1, 1]
        for i, attr in enumerate(self):
            vals[i] = max(1, attr - 1)
        return Shape(*vals)


class Coordinates(NamedTuple):
    # Store 3 1D coordinates arrays and one 'array_shape'
    x1: np.ndarray
    x2: np.ndarray
    x3: np.ndarray
    array_shape: Shape  # the 3D shape that 1D arrays should be broadcasted to

    @property
    def shape(self) -> Shape:
        return Shape(len(self.x1), len(self.x2), len(self.x3))

    @property
    def arrays(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        return self.x1, self.x2, self.x3

    def padded(self) -> Coordinates:
        arrs = [_.copy() for _ in self.arrays]
        for i, arr in enumerate(self.arrays):
            if arr.size == 1:
                arrs[i] = np.array((arr[0], arr[0] + 1))
        return Coordinates(arrs[0], arrs[1], arrs[2], self.array_shape)


def get_native_coordinates_from_cartesian(
    xcart: np.ndarray,
    ycart: np.ndarray,
    zcart: np.ndarray,
    geometry: str,
) -> list[np.ndarray]:
    shape = Shape(*xcart.shape)

    # Reconstruct the polar coordinate system
    if geometry == "polar":
        # on disk coordinates are cell face coordinates.
        # We need to convert from cartesian to polar,
        # but no interpolation is needed.
        r = np.sqrt(xcart[:, 0, 0] ** 2 + ycart[:, 0, 0] ** 2)
        theta = np.unwrap(np.arctan2(ycart[0, :, 0], xcart[0, :, 0]))
        z = zcart[0, 0, :]

        coords = [r, theta, z]
    elif geometry == "spherical":
        # Reconstruct the spherical coordinate system
        if shape.n3 == 1:
            r = np.sqrt(xcart[:, 0, 0] ** 2 + zcart[:, 0, 0] ** 2)
            phi = np.unwrap(
                np.arctan2(ycart[0, shape.n2 // 2, :], xcart[0, shape.n2 // 2, :])
            )
            theta = np.arccos(
                zcart[0, :, 0] / np.sqrt(xcart[0, :, 0] ** 2 + zcart[0, :, 0] ** 2)
            )
        else:
            r = np.sqrt(xcart[:, 0, 0] ** 2 + ycart[:, 0, 0] ** 2 + zcart[:, 0, 0] ** 2)
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
    else:
        raise NotImplementedError(
            f"This kind of geometry: {geometry} is not supported yet!"
        )
    return coords


# map field name to numpy array init data:
# precision (-> datatype), dimensionality, [nx, ny, nz]
# the np.ndarray is assumed to contain *dim* elements
IdefixFieldProperties = dict[str, tuple[Prec, Dim, np.ndarray]]

# Map various str keys to scalars and arrays
IdefixMetadata = dict[str, Any]
