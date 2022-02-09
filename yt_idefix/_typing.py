# this module defines some type constants that could be upstreamed
# so we don't want to use `from __future__ import annotation` here as long as yt doesn't
from typing import Tuple, Union

from unyt import unyt_quantity

# an intentionally restrictive list of types that can
# be passes to ds.quan (which is a proxy for unyt.unyt_quantity.__init__)
UnitLike = Union[unyt_quantity, Tuple[float, str]]
