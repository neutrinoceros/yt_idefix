from __future__ import annotations

import logging
import os
import re
import warnings
import weakref
from abc import ABC, abstractmethod
from functools import cached_property
from typing import Literal, Sequence

import inifix
import numpy as np
from packaging.version import Version

import yt
from yt.data_objects.static_output import Dataset
from yt.funcs import setdefaultattr
from yt.geometry.grid_geometry_handler import GridIndex
from yt_idefix._typing import UnitLike

from ._io import C_io, dmp_io, vtk_io
from ._io.commons import IdefixFieldProperties, IdefixMetadata
from .definitions import _PlutoBaseUnits, pluto_def_constants
from .fields import BaseVtkFields, IdefixDmpFields, IdefixVtkFields, PlutoVtkFields

try:
    from yt.data_objects.index_subobjects.stretched_grid import StretchedGrid
except ImportError:
    from ._vendors.streched_grids import StretchedGrid  # type: ignore [no-redef]

try:
    from yt.utilities.lib.misc_utilities import _obtain_coords_and_widths
except ImportError:
    from ._vendors.streched_grids import (
        _obtain_coords_and_widths,  # type: ignore [no-redef]
    )

# import IO classes to ensure they are properly registered,
# even though we don't call them directly
from .io import IdefixDmpIO, IdefixVtkIO, PlutoVtkIO  # noqa

ytLogger = logging.getLogger("yt")

YT_VERSION = Version(yt.__version__)


class IdefixGrid(StretchedGrid):
    _id_offset = 0

    def __init__(self, id, cell_widths, filename, index, level, dims):
        super().__init__(id=id, filename=filename, index=index, cell_widths=cell_widths)
        self.Parent = None
        self.Children = []
        self.Level = level
        self.ActiveDimensions = dims

    def __repr__(self):
        if YT_VERSION >= Version("4.1"):
            # https://github.com/yt-project/yt/pull/3936
            return super().__repr__()
        else:
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

        assert self.num_grids == 1

        i = 0
        g = self.grid(
            id=i,
            index=self,
            filename=self.index_filename,
            cell_widths=self._cell_widths,
            level=self.grid_levels.flat[i],
            dims=self.grid_dimensions[i],
        )
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

    @abstractmethod
    @cached_property
    def _cell_widths(self):
        # must return a 3-tuple of 1D unyt_array
        # with unit "code_length" and dtype float64
        ...

    @abstractmethod
    @cached_property
    def _cell_centers(self):
        # must return a 3-tuple of 1D unyt_array
        # with unit "code_length" and dtype float64
        ...

    def _icoords_to_fcoords(self, icoords, ires, axes: Sequence[int] = (0, 1, 2)):
        # this is needed to support projections
        coords = np.empty(icoords.shape, dtype="f8")
        cell_widths = np.empty(icoords.shape, dtype="f8")
        for i, ax in enumerate(axes):
            coords[:, i], cell_widths[:, i] = _obtain_coords_and_widths(
                icoords[:, i],
                ires,
                self._cell_widths[ax],
                self.ds.domain_left_edge[ax].d,
            )
        return coords, cell_widths


class IdefixVtkHierarchy(IdefixHierarchy):
    def _get_field_offset_index(self) -> dict[str, int]:
        return self.ds._field_offset_index

    @cached_property
    def _cell_widths(self):
        with open(self.index_filename, "rb") as fh:
            cell_edges = vtk_io.read_grid_coordinates(fh, geometry=self.ds.geometry)

        cell_widths: list[np.ndarray] = []
        length_unit = self.ds.quan(1, "code_length")
        for idir, edges in enumerate(cell_edges[:3]):
            ncells = self.ds.domain_dimensions[idir]
            if ncells > 1:
                cell_widths.append(np.ediff1d(edges).astype("float64") * length_unit)
            else:
                cell_widths.append(np.array([self.ds.domain_width[idir]]) * length_unit)
        return tuple(cell_widths)

    @cached_property
    def _cell_centers(self):
        with open(self.index_filename, "rb") as fh:
            cell_edges = vtk_io.read_grid_coordinates(fh, geometry=self.ds.geometry)

        cell_centers: list[np.ndarray] = []
        length_unit = self.ds.quan(1, "code_length")
        for idir, edges in enumerate(cell_edges[:3]):
            ncells = self.ds.domain_dimensions[idir]
            if ncells > 1:
                e64 = edges.astype("float64")
                cell_centers.append(0.5 * (e64[1:] + e64[:-1]) * length_unit)
            else:
                cell_centers.append(np.array([edges[0]]) * length_unit)
        return tuple(cell_centers)


class IdefixDmpHierarchy(IdefixHierarchy):
    def _get_field_offset_index(self) -> dict[str, int]:
        with open(self.index_filename, "rb") as fh:
            return dmp_io.get_field_offset_index(fh)

    @cached_property
    def _cell_widths(self):
        _fprops, fdata = dmp_io.read_idefix_dmpfile(self.index_filename, skip_data=True)
        length_unit = self.ds.quan(1, "code_length")
        return tuple((fdata[f"xr{d}"] - fdata[f"xl{d}"]) * length_unit for d in "123")

    @cached_property
    def _cell_centers(self):
        _fprops, fdata = dmp_io.read_idefix_dmpfile(self.index_filename, skip_data=True)
        length_unit = self.ds.quan(1, "code_length")
        return tuple(fdata[f"x{d}"] * length_unit for d in "123")


class IdefixDataset(Dataset, ABC):
    """A common abstraction for IdefixDmpDataset and IdefixVtkDataset."""

    _version_regexp = re.compile(r"v\d+\.\d+\.?\d*[-\w+]*")
    _dataset_type: str  # defined in subclasses

    def __init__(
        self,
        filename,
        *,
        dataset_type: str | None = None,  # deleguated to child classes
        file_style: str | None = None,
        units_override: dict[str, UnitLike] | None = None,
        unit_system: Literal["cgs", "mks", "code"] = "cgs",
        default_species_fields: Literal["neutral", "ionized"] | None = None,
        # from here, frontend-specific arguments
        geometry: Literal["cartesian", "spherical", "cylindrical", "polar"]
        | None = None,
        inifile: str | os.PathLike[str] | None = None,
    ):
        self._geometry_from_user = geometry

        dt = type(self)._dataset_type
        self.fluid_types += (dt,)

        super().__init__(
            filename,
            dataset_type=dt,
            file_style=file_style,
            units_override=units_override,
            unit_system=unit_system,
            default_species_fields=default_species_fields,
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

        self._setup_geometry()

    def _setup_geometry(self) -> None:
        from_file = self.parameters.get("geometry")
        from_user = self._geometry_from_user
        if from_file is None and from_user is None:
            raise ValueError(
                "Geometry couldn't be parsed from file. "
                "The 'geometry' keyword argument must be specified."
            )
        elif from_user is not None:
            if from_file is not None and from_user != from_user:
                warnings.warn(
                    f"Got inconsistent geometry flags:\n"
                    f" - {from_file!r} (from file)\n"
                    f" - {from_user!r} (from user)\n"
                    "user-input prevails to allow working around hypothetical parsing bugs, "
                    "but it is very likely to result in an error in the general case."
                )
            self.geometry = from_user
        else:
            assert from_file is not None
            self.geometry = from_file

    def _parse_inifile(self) -> None:
        if self.inifile is None:
            return

        with open(self.inifile, "rb") as fh:
            self.parameters.update(inifix.load(fh))
        grid_ini = self.parameters["Grid"]

        msg_elems: list[str] = []
        for ax, vals in grid_ini.items():
            if vals[0] > 1:
                # more than one block is only relevant for mixing grid spacings,
                # but only "u" is supported
                msg_elems.append(f"found multiple blocks in direction {ax}; got {vals}")
            if any(_ != "u" for _ in vals[3::3]):
                msg_elems.append(f"found non-uniform block(s) in direction {ax}")

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
    _field_info_class: type[BaseVtkFields] = IdefixVtkFields
    _dataset_type = "idefix-vtk"
    _required_header_keyword = "Idefix"

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs) -> bool:
        try:
            header = vtk_io.read_header(filename)
        except Exception:
            return False
        else:
            return cls._required_header_keyword in header

    def _get_header(self) -> str:
        return vtk_io.read_header(self.parameter_filename)

    def _parse_parameter_file(self):

        # parse metadata
        with open(self.parameter_filename, "rb") as fh:
            md = vtk_io.read_metadata(fh)
        self.parameters.update(md)

        super()._parse_parameter_file()
        # from here self.geometry is assumed to be set

        # parse the grid
        with open(self.parameter_filename, "rb") as fh:
            coords = vtk_io.read_grid_coordinates(fh, geometry=self.geometry)
            self._field_offset_index = vtk_io.read_field_offset_index(
                fh, coords.array_shape
            )
        self._detected_field_list = list(self._field_offset_index.keys())

        self.domain_dimensions = np.array(coords.array_shape)
        self.dimensionality = np.count_nonzero(self.domain_dimensions - 1)

        dle = np.array([arr.min() for arr in coords.arrays], dtype="float64")
        dre = np.array([arr.max() for arr in coords.arrays], dtype="float64")

        # temporary hack to prevent 0-width dimensions for 2D data
        dre = np.where(dre == dle, dle + 1, dre)
        self.domain_left_edge = dle
        self.domain_right_edge = dre

        # time wasn't stored in vtk files before Idefix 0.8
        self.current_time = md.get("time", -1)

        # periodicity was not stored in vtk files before Idefix 0.9
        self._periodicity = md.get("periodicity", (True, True, True))


class IdefixDmpDataset(IdefixDataset):
    _index_class = IdefixDmpHierarchy
    _field_info_class = IdefixDmpFields
    _dataset_type = "idefix-dmp"

    @classmethod
    def _is_valid(cls, filename, *args, **kwargs) -> bool:
        try:
            header_string = dmp_io.read_header(filename)
            return re.match(r"Idefix .* Dump Data", header_string) is not None
        except Exception:
            return False

    def _get_fields_metadata(self) -> tuple[IdefixFieldProperties, IdefixMetadata]:
        # read everything except large arrays
        return dmp_io.read_idefix_dmpfile(self.parameter_filename, skip_data=True)

    def _get_header(self) -> str:
        return dmp_io.read_header(self.parameter_filename)

    def _parse_parameter_file(self):

        fprops, fdata = self._get_fields_metadata()
        self.parameters.update(fdata)

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

        super()._parse_parameter_file()


class PlutoVtkDataset(IdefixVtkDataset):
    _field_info_class = PlutoVtkFields
    _dataset_type = "pluto-vtk"
    _version_regexp = re.compile(r"\d+\.\d+\.?\d*[-\w+]*")
    _required_header_keyword = "PLUTO"

    def __init__(
        self,
        filename,
        *,
        dataset_type: str | None = None,  # deleguated to child classes
        file_style: str | None = None,  # unused
        units_override: dict[str, UnitLike] | None = None,
        unit_system: Literal["cgs", "mks", "code"] = "cgs",
        default_species_fields: Literal["neutral", "ionized"] | None = None,
        # from here, frontend-specific arguments
        geometry: Literal["cartesian", "spherical", "cylindrical", "polar"]
        | None = None,
        inifile: str | os.PathLike[str] | None = None,
        definitions_header: str | None = None,
    ):
        self._definitions_header: str | None
        if definitions_header is not None:
            self._definitions_header = os.fspath(definitions_header)
        else:
            self._definitions_header = None

        if YT_VERSION < Version("4.1.dev0"):
            # https://github.com/yt-project/yt/pull/3772
            filename = os.path.abspath(os.path.expanduser(filename))

        super().__init__(
            filename,
            dataset_type=dataset_type,
            file_style=file_style,
            units_override=units_override,
            unit_system=unit_system,
            geometry=geometry,
            inifile=inifile,
            default_species_fields=default_species_fields,
        )

    def _parse_parameter_file(self):
        self._parse_header_file()
        super()._parse_parameter_file()
        self._get_time()

    def _parse_header_file(self):
        """Read some metadata from header file 'definitions.h'."""
        geom_regexp = re.compile(r"^\s*#define\s+GEOMETRY\s+([A-Z]+)")
        unit_regexp = re.compile(r"^\s*#define\s+UNIT_(\w+)\s+(\S+)")
        constexpr = re.compile(r"CONST_\w+")

        # definitions.h is presumed to be along with data file
        if self._definitions_header is None:
            self._definitions_header = os.path.join(self.directory, "definitions.h")
        elif not os.path.isfile(self._definitions_header):
            raise FileNotFoundError(f"No such file {self._definitions_header!r}")

        if os.path.isfile(self._definitions_header):
            with open(self._definitions_header) as fh:
                body = fh.read()
            lines = C_io.strip_comments(body).split("\n")

            for line in lines:
                geom_match = re.fullmatch(geom_regexp, line)
                if geom_match is not None:
                    self.parameters["geometry"] = geom_match.group(1).lower()
                    continue

                unit_match = re.fullmatch(unit_regexp, line)
                if unit_match is not None:
                    unit = unit_match.group(1).lower() + "_unit"
                    expr = unit_match.group(2)
                    expr = re.sub(constexpr, self._get_constants, expr)
                    self.parameters[unit] = eval(expr)
        else:
            warnings.warn(
                f"Header file {self._definitions_header} couldn't be found. "
                "The code units are set to be 1.0 in cgs by default."
            )

    def _get_constants(self, match: re.Match) -> str:
        """Replace matched constant string with its value"""
        key = match.group()
        return str(pluto_def_constants[key])

    def _get_time(self):
        """Get current time from vtk.out."""
        log_file = os.path.join(self.directory, "vtk.out")
        match = re.search(r"\.(\d*)\.", self.parameter_filename)
        if match is None:
            raise RuntimeError(
                f"Failed to parse output number from file name {self.parameter_filename}"
            )
        index = int(match.group(1))

        # will be converted to actual unyt_quantity in _set_derived_attrs
        self.current_time = -1
        if os.path.isfile(log_file):
            log_regexp = re.compile(rf"^{index}\s(\S+)")
            with open(log_file) as fh:
                for line in fh.readlines():
                    log_match = re.search(log_regexp, line)
                    if log_match:
                        self.current_time = float(log_match.group(1))
                        break
                else:
                    ytLogger.warning(
                        "Failed to retrieve time from %s, setting current_time = -1",
                        log_file,
                    )
        else:
            ytLogger.warning("Missing log file %s, setting current_time = -1", log_file)

    def _set_code_unit_attributes(self):
        """Conversion between physical units and code units."""

        # Pluto's base units are length, velocity and density, but here we consider
        # length, mass and time as base units. Since it can make us easy to calculate
        # all units when self.units_override is not None.

        # Default values of Pluto's base units which are stored in self.parameters
        # if they can be read from definitions.h
        # Otherwise, they are set to unity in cgs.
        pluto_units = {
            "velocity_unit": self.quan(
                self.parameters.get("velocity_unit", 1.0), "cm/s"
            ),
            "density_unit": self.quan(
                self.parameters.get("density_unit", 1.0), "g/cm**3"
            ),
            "length_unit": self.quan(self.parameters.get("length_unit", 1.0), "cm"),
        }

        uo_size = len(self.units_override)
        if uo_size > 0 and uo_size < 3:
            ytLogger.info(
                "Less than 3 units were specified in units_override (got %s). "
                "Need to rely on PLUTO's internal units to derive other units",
                uo_size,
            )

        uo_cache = self.units_override.copy()
        while len(uo_cache) < 3:
            # If less than 3 units were passed into units_override,
            # the rest will be chosen from Pluto's units
            unit, value = pluto_units.popitem()
            # If any Pluto's base unit is specified in units_override, it'll be preserved
            if unit in uo_cache:
                continue
            uo_cache[unit] = value
            # Make sure the combination of units are able to derive base units
            # No need of validation and logging when no unit to be overrided
            if uo_size > 0:
                try:
                    self._validate_units_override_keys(uo_cache)
                except ValueError:
                    # It means the combination is invalid
                    del uo_cache[unit]
                else:
                    ytLogger.info("Relying on %s: %s.", unit, uo_cache[unit])

        bu = _PlutoBaseUnits(uo_cache)
        for unit, value in bu._data.items():
            setattr(self, unit, value)

        self.velocity_unit = self.length_unit / self.time_unit
        self.density_unit = self.mass_unit / self.length_unit**3
        self.magnetic_unit = (
            np.sqrt(4.0 * np.pi * self.density_unit) * self.velocity_unit
        )
        self.magnetic_unit.convert_to_units("gauss")
        self.temperature_unit = self.quan(1.0, "K")

    invalid_unit_combinations = [
        {"magnetic_unit", "velocity_unit", "density_unit"},
        {"velocity_unit", "time_unit", "length_unit"},
        {"density_unit", "length_unit", "mass_unit"},
    ]

    default_units = {
        "length_unit": "cm",
        "time_unit": "s",
        "mass_unit": "g",
        "velocity_unit": "cm/s",
        "magnetic_unit": "gauss",
        "temperature_unit": "K",
        # this is the one difference with Dataset.default_units:
        # we accept density_unit as a valid override
        "density_unit": "g/cm**3",
    }

    @classmethod
    def _validate_units_override_keys(cls, units_override):
        """Check that units in units_override are able to derive three base units:
        mass, length and time
        """

        # YT supports overriding other normalisations, this method ensures consistency
        # between supplied 'units_override' items and principles in PLUTO.

        # PLUTO's normalisations/units have 3 degrees of freedom. Therefore, any combinations
        # are valid other than the three cases listed explicitly.

        if "temperature_unit" in units_override:
            raise ValueError(
                "Temperature is not allowed in units_override, "
                "since it's always in Kelvin in PLUTO"
            )

        # Three units are enough for deriving others, more will likely cause conflict
        if len(units_override) > 3:
            raise ValueError(
                "More than 3 degrees of freedom were specified "
                f"in units_override ({len(units_override)} given)"
            )

        # check if provided overrides are allowed
        suo = set(units_override)
        if suo in cls.invalid_unit_combinations:
            raise ValueError(
                f"Combination {suo} passed to units_override "
                "cannot derive all units\n"
                f"Choose any other combinations, except for:\n {cls.invalid_unit_combinations}"
            )

        super(cls, cls)._validate_units_override_keys(units_override)
