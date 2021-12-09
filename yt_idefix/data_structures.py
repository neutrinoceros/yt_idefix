from __future__ import annotations

import os
import re
import warnings
import weakref
from abc import ABC, abstractmethod

import inifix
import numpy as np

from yt.data_objects.index_subobjects.grid_patch import AMRGridPatch
from yt.data_objects.static_output import Dataset
from yt.funcs import setdefaultattr
from yt.geometry.grid_geometry_handler import GridIndex

from ._io import dmp_io, vtk_io
from ._io.commons import IdefixFieldProperties, IdefixMetadata
from .fields import IdefixDmpFieldInfo, IdefixVtkFieldInfo


class IdefixGrid(AMRGridPatch):
    _id_offset = 0

    def __init__(self, id, index, level, dims):
        super().__init__(id, filename=index.index_filename, index=index)
        self.Parent = None
        self.Children = []
        self.Level = level
        self.ActiveDimensions = dims

    def __repr__(self):
        return "IdefixGrid_%04i (%s)" % (self.id, self.ActiveDimensions)


class IdefixHierarchy(GridIndex, ABC):
    grid = IdefixGrid

    def __init__(self, ds, dataset_type="idefix"):
        self.dataset_type = dataset_type
        self.dataset = weakref.proxy(ds)
        # for now, the index file is the dataset!
        self.index_filename = self.dataset.parameter_filename
        self.directory = os.path.dirname(self.index_filename)
        # float type for the simulation edges and must be float64 now
        self.float_type = np.float64
        super().__init__(ds, dataset_type)

    def _detect_output_fields(self):
        self.field_list = [
            (self.dataset_type, f) for f in self.dataset._detected_field_list
        ]

    def _count_grids(self):
        self.num_grids = 1

    def _parse_index(self):
        self.grid_left_edge[0][:] = self.ds.domain_left_edge[:]
        self.grid_right_edge[0][:] = self.ds.domain_right_edge[:]
        self.grid_dimensions[0][:] = self.ds.domain_dimensions[:]
        self.grid_particle_count[0][0] = 0
        self.grid_levels[0][0] = 1
        self.max_level = 1

        self._field_offsets = self._get_field_offset_index()

    def _populate_grid_objects(self):
        # the minimal form of this method is
        #
        # for g in self.grids:
        #     g._prepare_grid()
        #     g._setup_dx()
        #
        # This must also set:
        #   g.Children <= list of child grids
        #   g.Parent   <= parent grid
        # This is handled by the frontend because often the children must be identified.
        self.grids = np.empty(self.num_grids, dtype="object")
        for i in range(self.num_grids):
            g = self.grid(i, self, self.grid_levels.flat[i], self.grid_dimensions[i])
            g._prepare_grid()
            g._setup_dx()
            self.grids[i] = g

    @abstractmethod
    def _get_field_offset_index(self) -> dict[str, int]:
        HEADER_SIZE: int = 256
        with open(self.index_filename, "rb") as fh:
            fh.seek(HEADER_SIZE)
            field_index = ...
        return field_index  # type: ignore


class IdefixVtkHierarchy(IdefixHierarchy):
    def _get_field_offset_index(self) -> dict[str, int]:
        return self.ds._field_offset_index


class IdefixDmpHierarchy(IdefixHierarchy):
    def _get_field_offset_index(self) -> dict[str, int]:
        with open(self.index_filename, "rb") as fh:
            return dmp_io.get_field_offset_index(fh)


class IdefixDataset(Dataset, ABC):
    """A common abstraction for IdefixDmpDataset and IdefixVtkDataset."""

    _version_regexp = re.compile(r"v\d+\.\d+\.?\d*[-\w+]*")

    def __init__(
        self,
        filename,
        dataset_type=None,
        unit_system="cgs",
        units_override=None,
        *,
        geometry: str | None = None,
        inifile=None,
    ):
        dt = type(self)._dataset_type
        self.fluid_types += (dt,)

        self.geometry = geometry
        super().__init__(
            filename,
            dataset_type=dt,
            units_override=units_override,
            unit_system=unit_system,
        )
        self.inifile = inifile
        self._parse_inifile()

        self.storage_filename = None

        # idefix does not support grid refinement
        self.refine_by = 1

    def _parse_parameter_file(self):
        # base method, intended to be subclassed
        # parse the version hash
        self.parameters["code version"] = self._get_code_version()

        # idefix is never cosmological
        self.cosmological_simulation = 0
        self.current_redshift = 0.0
        self.omega_lambda = 0.0
        self.omega_matter = 0.0
        self.hubble_constant = 0.0

    def _parse_inifile(self) -> None:
        if self.inifile is None:
            return

        self.parameters.update(inifix.load(self.inifile))
        grid_ini = self.parameters["Grid"]

        msg_elems: list[str] = []
        for ax, vals in grid_ini.items():
            if vals[0] > 1:
                # more than one block is only relevant for mixing grid spacings,
                # but only "u" is supported
                msg_elems.append(f"found multiple blocks in direction {ax}; got {vals}")
            if any(_ != "u" for _ in vals[3::3]):
                msg_elems.append(f"found non-uniform block(s) in direction {ax}")
        if len(msg_elems) > 0:
            msg = (
                "Streched grid detected !\n"
                + "- "
                + "\n- ".join(msg_elems)
                + "\nThe grid will be treated as uniformly spaced in every direction. "
                "yt_idefix.load_stretched may be a better alternative loader for these data"
            )
            warnings.warn(msg)

    def _set_code_unit_attributes(self):
        # This is where quantities are created that represent the various
        # on-disk units.  These are the currently available quantities which
        # should be set, along with examples of how to set them to standard
        # values.
        #
        # self.length_unit = self.quan(1.0, "cm")
        # self.mass_unit = self.quan(1.0, "g")
        # self.time_unit = self.quan(1.0, "s")
        # self.time_unit = self.quan(1.0, "s")
        #
        # These can also be set:
        # self.velocity_unit = self.quan(1.0, "cm/s")
        # self.magnetic_unit = self.quan(1.0, "gauss")
        for key, unit in self.__class__.default_units.items():
            setdefaultattr(self, key, self.quan(1, unit))

    # The following methods are frontend-specific

    @abstractmethod
    def _get_header(self) -> str:
        pass

    def _get_code_version(self) -> str:
        # take the last line of the header
        # - in Idefix dumps there's only one line
        # - in Vtk files (Idefix or Pluto), there are two,
        #   the first of which isn't code specific
        header = self._get_header().splitlines()[-1]

        regexp = self.__class__._version_regexp

        match = re.search(regexp, header)
        version: str
        if match is None:
            warnings.warn(
                f"Could not determine code version from file header {header!r}"
            )
            return "unknown"

        return match.group()


class IdefixVtkDataset(IdefixDataset):
    _index_class = IdefixVtkHierarchy
    _field_info_class = IdefixVtkFieldInfo
    _dataset_type = "idefix-vtk"
    _required_header_keyword = "Idefix"

    @classmethod
    def _is_valid(cls, fn, *args, **kwargs) -> bool:
        try:
            header = vtk_io.read_header(fn)
        except Exception:
            return False
        else:
            return cls._required_header_keyword in header

    def _get_header(self) -> str:
        return vtk_io.read_header(self.parameter_filename)

    def _parse_parameter_file(self):
        super()._parse_parameter_file()

        # parse the grid
        with open(self.parameter_filename, "rb") as fh:
            # at this point, self.geometry represents the user input (possibly None)
            md = vtk_io.read_metadata(fh, geometry=self.geometry)
            coords = vtk_io.read_grid_coordinates(fh, md)
            self._field_offset_index = vtk_io.read_field_offset_index(
                fh, md["array_shape"]
            )

        self.parameters.update(md)

        self._detected_field_list = list(self._field_offset_index.keys())

        self.domain_dimensions = np.array(md["array_shape"])
        self.dimensionality = np.count_nonzero(self.domain_dimensions - 1)

        dle = np.array([arr.min() for arr in coords], dtype="float64")
        dre = np.array([arr.max() for arr in coords], dtype="float64")

        # temporary hack to prevent 0-width dimensions for 2D data
        dre = np.where(dre == dle, dle + 1, dre)
        self.domain_left_edge = dle
        self.domain_right_edge = dre

        # time wasn't stored in vtk files before Idefix 0.8
        self.current_time = md.get("time", -1)
        # periodicity was not stored in vtk files before Idefix 0.9
        self._periodicity = md.get("periodicity", (True, True, True))
        self.geometry = md["geometry"]


class IdefixDmpDataset(IdefixDataset):
    _index_class = IdefixDmpHierarchy
    _field_info_class = IdefixDmpFieldInfo
    _dataset_type = "idefix-dmp"

    @classmethod
    def _is_valid(cls, fn, *args, **kwargs):
        ok = bool(
            re.match(r"^(dump)\.\d{4}(\.dmp)$", os.path.basename(fn))
        )  # this is possibly too restrictive
        try:
            ok &= "Idefix" in dmp_io.read_header(fn)
        except Exception:
            ok = False
        return ok

    def _get_fields_metadata(self) -> tuple[IdefixFieldProperties, IdefixMetadata]:
        # read everything except large arrays
        return dmp_io.read_idefix_dmpfile(self.parameter_filename, skip_data=True)

    def _get_header(self) -> str:
        return dmp_io.read_header(self.parameter_filename)

    def _parse_parameter_file(self):
        super()._parse_parameter_file()

        fprops, fdata = self._get_fields_metadata()
        self._detected_field_list = [k for k in fprops if re.match(r"^V[sc]-", k)]

        # parse the grid
        axes = ("x1", "x2", "x3")
        self.domain_dimensions = np.concatenate([fprops[k][-1] for k in axes])
        self.dimensionality = np.count_nonzero(self.domain_dimensions - 1)

        # note that domain edges parsing is already implemented in a mutli-block
        # supporting fashion even though we specifically error out in case there's more
        # than one block.
        self.domain_left_edge = np.array(
            [fdata[f"xl{idir}"][0] for idir in "123"], dtype="float64"
        )
        self.domain_right_edge = np.array(
            [fdata[f"xr{idir}"][-1] for idir in "123"], dtype="float64"
        )

        self.current_time = fdata["time"]

        self._periodicity = tuple(bool(p) for p in fdata["periodicity"])
        self.geometry = fdata["geometry"]


class PlutoVtkDataset(IdefixVtkDataset):
    _version_regexp = re.compile(r"\d+\.\d+\.?\d*[-\w+]*")
    _required_header_keyword = "PLUTO"
