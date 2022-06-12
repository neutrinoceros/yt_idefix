from abc import ABC, abstractmethod
from typing import BinaryIO, DefaultDict, List, Tuple, cast

import numpy as np
from packaging.version import Version

import yt
from yt.utilities.io_handler import BaseIOHandler

from ._io import dmp_io, vtk_io

YT_VERSION = Version(yt.__version__)
if YT_VERSION >= Version("4.1"):
    from yt.utilities.io_handler import BaseParticleIOHandler
else:
    # this class is vendored from yt 4.1dev
    class BaseParticleIOHandler(BaseIOHandler):  # type:ignore [no-redef]
        def _sorted_chunk_iterator(self, chunks):
            chunks = list(chunks)
            data_files = set()
            for chunk in chunks:
                for obj in chunk.objs:
                    data_files.update(obj.data_files)
            yield from sorted(data_files, key=lambda x: (x.filename, x.start))

        def _count_particles_chunks(
            self,
            psize: DefaultDict[str, int],
            chunks,
            ptf: DefaultDict[str, List[str]],
            selector,
        ) -> DefaultDict[str, int]:
            if getattr(selector, "is_all_data", False):
                for data_file in self._sorted_chunk_iterator(chunks):
                    for ptype in ptf.keys():
                        psize[ptype] += data_file.total_particles[ptype]
            else:
                # we must apply the selector and count the result
                psize = self._count_selected_particles(psize, chunks, ptf, selector)
            return psize


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
        # This will read chunks and yield the results.
        ptype = "dust"

        if "PX1" not in grid._index._field_offsets:
            return

        # TODO: define particle count
        if particle_count == 0:
            return

        cl = self.ds.quan(1, "code_length")
        with open(self._data_file, mode="rb") as fh:
            foffset = grid._index._field_offsets["PX1"]
            x = self._read_single_field(fh, foffset) * cl
            if self.ds.dimensionality >= 2:
                foffset = grid._index._field_offsets["PX2"]
                y = self._read_single_field(fh, foffset) * cl
            else:
                y = np.full_like(x, self.ds.domain_center[1])
            if self.ds.dimensionality == 3:
                foffset = grid._index._field_offsets["PX3"]
                z = self._read_single_field(fh, foffset) * cl
            else:
                z = np.full_like(x, self.ds.domain_center[2])

        yield ptype, (x, y, z), 0.0

    def _read_particle_fields(self, chunks, ptf, selector):
        # idefix doesn't have particles (yet)
        raise NotImplementedError("Particles fields are not implemented yet")
