from abc import ABC, abstractmethod
from typing import BinaryIO, Tuple, cast

import numpy as np

from yt.utilities.io_handler import BaseIOHandler

from ._io import dmp_io, vtk_io

try:
    from yt.utilities.io_handler import BaseParticleIOHandler
except ImportError:
    from ._vendors.base_particle_io_handler import (
        BaseParticleIOHandler,  # type: ignore [no-redef]
    )


class SingleGridIO(BaseIOHandler, ABC):
    _particle_reader = False

    def __init__(self, ds):
        BaseIOHandler.__init__(self, ds)
        self.ds = ds
        self._data_file = ds.parameter_filename

    def _read_fluid_selection(self, chunks, selector, fields, size):
        data = {}

        for field in fields:
            data[field] = np.empty(size, dtype="float64")

        with open(self._data_file, mode="rb") as fh:
            ind = 0
            for chunk in chunks:
                for grid in chunk.objs:
                    nd = 0
                    for field in fields:
                        ftype, fname = field
                        foffset = grid._index._field_offsets[fname]
                        values = self._read_single_field(fh, foffset)
                        nd = grid.select(selector, values, data[field], ind)
                    ind += nd
        return data

    def _read_chunk_data(self, chunk, fields):
        # This reads the data from a single chunk without doing any selection,
        # and is only used for caching data that might be used by multiple
        # different selectors later. For instance, this can speed up ghost zone
        # computation.
        raise NotImplementedError

    # The following methods are frontend-specific
    @abstractmethod
    def _read_single_field(self, fh: BinaryIO, offset: int) -> np.ndarray:
        pass


class PlutoVtkIO(SingleGridIO):
    _dataset_type = "pluto-vtk"

    def _read_single_field(self, fh: BinaryIO, offset: int) -> np.ndarray:
        shape = cast(Tuple[int, int, int], tuple(self.ds.domain_dimensions))
        return vtk_io.read_single_field(fh, shape=shape, offset=offset, skip_data=False)


class IdefixVtkIO(PlutoVtkIO, BaseParticleIOHandler):
    _dataset_type = "idefix-vtk"

    def _read_particle_coords(self, chunks, ptf):
        # This needs to *yield* a series of tuples of (ptype, (x, y, z)).
        # chunks is a list of chunks, and ptf is a dict where the keys are
        # ptypes and the values are lists of fields.
        raise NotImplementedError("Particles are not currently supported for Idefix")

    def _read_particle_fields(self, chunks, ptf, selector):
        # idefix doesn't have particles (yet)
        raise NotImplementedError("Particles are not currently supported for Idefix")


class IdefixDmpIO(SingleGridIO, BaseParticleIOHandler):
    _dataset_type = "idefix-dmp"

    def _read_single_field(self, fh: BinaryIO, offset: int) -> np.ndarray:
        return dmp_io.read_single_field(fh, offset)

    def _read_particle_coords(self, chunks, ptf):
        # This needs to *yield* a series of tuples of (ptype, (x, y, z)).
        # chunks is a list of chunks, and ptf is a dict where the keys are
        # ptypes and the values are lists of fields.
        raise NotImplementedError("Particles are not currently supported for Idefix")

    def _read_particle_fields(self, chunks, ptf, selector):
        # idefix doesn't have particles (yet)
        raise NotImplementedError("Particles are not currently supported for Idefix")
