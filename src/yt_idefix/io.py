import re
from abc import ABC, abstractmethod
from typing import BinaryIO, cast

import numpy as np

from yt.utilities.io_handler import BaseIOHandler, BaseParticleIOHandler
from yt.utilities.on_demand_imports import _h5py as h5py

from ._io import dmp_io, vtk_io


class SingleGridIO(BaseIOHandler, ABC):
    _particle_reader = False

    def _read_fluid_selection(self, chunks, selector, fields, size):
        data = {}

        for field in fields:
            data[field] = np.empty(size, dtype="float64")

        with open(self.ds.filename, mode="rb") as fh:
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
        shape = cast(tuple[int, int, int], tuple(self.ds.domain_dimensions))
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

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._byteorder = dmp_io.parse_byteorder(self.ds.filename)

    def _read_single_field(self, fh: BinaryIO, offset: int) -> np.ndarray:
        return dmp_io.read_single_field(fh, offset, byteorder=self._byteorder)

    def _read_particle_coords(self, chunks, ptf):
        # This needs to *yield* a series of tuples of (ptype, (x, y, z)).
        # chunks is a list of chunks, and ptf is a dict where the keys are
        # ptypes and the values are lists of fields.
        raise NotImplementedError("Particles are not currently supported for Idefix")

    def _read_particle_fields(self, chunks, ptf, selector):
        # idefix doesn't have particles (yet)
        raise NotImplementedError("Particles are not currently supported for Idefix")


class PlutoXdmfIOHandler(BaseIOHandler):
    _particle_reader = False
    _dataset_type = "pluto-xdmf"

    def _read_particle_coords(self, chunks, ptf):
        raise NotImplementedError("Reading of particle fields not yet implemented!")

    def _read_particle_fields(self, chunks, ptf, selector):
        raise NotImplementedError("Reading of particle fields not yet implemented!")

    def _read_fluid_selection(self, chunks, selector, fields, size):
        """
        Filenames are data.<snapnum>.<dbl/flt>.h5
        <snapnum> needs to be parse from the filename.
        <snapnum> is the corresponding entry in the <dbl/flt>.h5.out file
        Example <dbl/flt>.h5.out file:
            0 0.000000e+00 1.000000e-04 0 single_file little rho vx1 vx2 vx3 prs tr1 tr2 tr3 Temp ndens PbykB mach
            1 2.498181e+00 3.500985e-03 747 single_file little rho vx1 vx2 vx3 prs tr1 tr2 tr3 Temp ndens PbykB mach
            2 4.998045e+00 3.400969e-03 1458 single_file little rho vx1 vx2 vx3 prs tr1 tr2 tr3 Temp ndens PbykB mach
            3 7.497932e+00 3.386245e-03 2186 single_file little rho vx1 vx2 vx3 prs tr1 tr2 tr3 Temp ndens PbykB mach
        """
        if (match := re.search(r"\d{4}", self.ds.filename)) is not None:
            entry = int(match.group())
        else:
            raise RuntimeError(f"Failed to parse output number from {self.ds.filename}")

        data = {field: np.empty(size, dtype="float64") for field in fields}

        with h5py.File(self.ds.filename, "r") as fh:
            ind = 0
            for chunk in chunks:
                for grid in chunk.objs:
                    nd = 0
                    for field in fields:
                        _, fname = field
                        position = (
                            f"/Timestep_{entry}/vars/{self.ds._field_name_map[fname]}"
                        )
                        field_data = fh[position][:].astype("=f8")

                        while field_data.ndim < 3:
                            field_data = np.expand_dims(field_data, axis=0)

                        # X3 X2 X1 orderding of fields in PLUTO needs to rearranged to X1 X2 X3 order in yt.
                        values = np.transpose(field_data, axes=(2, 1, 0))
                        nd = grid.select(selector, values, data[field], ind)
                    ind += nd
        return data

    def _read_chunk_data(self, chunk, fields):
        raise NotImplementedError("No multipart file dumps happen in PLUTO till date!")
