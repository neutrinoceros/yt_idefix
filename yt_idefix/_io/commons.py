import sys
from typing import Any, Dict, Tuple

import numpy as np

if sys.version_info >= (3, 8):
    from typing import Literal

    Prec = Literal["d", "f", "i"]
    Dim = Literal[1, 2, 3]
else:
    Prec = str
    Dim = int

# map field name to numpy array init data:
# precision (-> datatype), dimensionality, [nx, ny, nz]
# the np.ndarray is assumed to contain *dim* elements
IdefixFieldProperties = Dict[str, Tuple[Prec, Dim, np.ndarray]]

# Map various str keys to scalars and arrays
IdefixMetadata = Dict[str, Any]
