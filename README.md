
# yt_idefix
[![PyPI](https://img.shields.io/pypi/v/yt-idefix.svg?logo=pypi&logoColor=white&label=PyPI)](https://pypi.org/project/yt_idefix/)
[![yt-project](https://img.shields.io/static/v1?label="works%20with"&message="yt"&color="blueviolet")](https://yt-project.org)

<!--- Tests and style --->
[![CI](https://github.com/neutrinoceros/yt_idefix/actions/workflows/ci.yml/badge.svg)](https://github.com/neutrinoceros/yt_idefix/actions/workflows/ci.yml)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/neutrinoceros/yt_idefix/main.svg)](https://results.pre-commit.ci/latest/github/neutrinoceros/yt_idefix/main)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/charliermarsh/ruff/main/assets/badge/v2.json)](https://github.com/charliermarsh/ruff)

A maturing yt frontend for Idefix and Pluto, packaged as an extension for yt.

## Installation

```shell
python -m pip install yt_idefix
```

## Supported formats

| Code   | format | supported since | additional dependencies |
|--------|:------:|:---------------:|-------------------------|
| Idefix | `.dmp` | v0.1.0          |                         |
| Idefix | `.vtk` | v0.3.0          |                         |
| Pluto  | `.vtk` | v0.9.0          |                         |
| Pluto  |  XDMF  | v1.1.0          | `h5py`                  |

## Usage

Integration with yt is seamless. *Installing* this plugin is all that's required to make yt
compatible with data formats supported by `yt_idefix` !

### Additional arguments to `yt.load`
The metadata are parsed from data file, definitions header file and inifile when loading dataset.

Definitions header file (`definitions.h` for Pluto, or `definitions.hpp` for Idefix) and inifile (`pluto.ini` and `idefix.ini` respectively) are discovered automatically if they match default names, are located along with data files, and unique. Otherwise, they can be specified explicitly as paths (either relative to data files or absolute paths) with parameters `definitions_header` and `inifile` respectively.

```python
ds = yt.load(
    "data.0010.vtk",
    definitions_header="../definitions.h",
    inifile="example.ini",
)
```

Geometry is parsed automatically whenever possible, but as a last resort, it can also be specified as a keyword argument (possible values are "cartesian", "spherical", "cylindrical" and "polar").

```python
ds = yt.load("data.0010.vtk", geometry="spherical")
```

The data are loaded as physical quantities with units. The default unit system is `cgs` in yt. Data is always interpreted as dimensionful.
For Pluto simulation, metadata is retrieved from `definitions.h` and `pluto.ini` to guess the proper on-disk units automatically.

Units may also be provided at runtime using the `units_override` argument
```python
ds = yt.load(
    "data.0010.vtk",
    units_override={
        "length_unit": (100.0, "au"),
        "mass_unit": yt.units.mass_sun,
    },
)
```
Note that other units will also be changed for consistency (Pluto).

Displayed units can also be controled using the `unit_system` argument.
Accepted values are `"cgs"` (default), `"mks"` and `"code"`.

```python
ds = yt.load("data.0010.vtk", unit_system="mks")
```

With Pluto data, units not specified with `units_override` will be derived consistently with given units, within the following rules:
1. Temperature unit cannot be overridden (always set to Kelvin)
2. No more than three units can be overridden at once (overconstrained systems are never validated for simplicity)
3. When given less than three overrides, base units in Pluto (ordered: velocity_unit, density_unit, length_unit) are assumed
4. The following combinations are not allowed

```python
{"magnetic_unit", "velocity_unit", "density_unit"},
{"velocity_unit", "time_unit", "length_unit"},
{"density_unit", "length_unit", "mass_unit"}
```

yt is able to provide some derived fields from existed fields, e.g., `"cell_volume"`. Fields related to element species can be created according to primordial abundances of H and He, through `default_species_fields` (`"neutral"` and `"ionized"`) parameters.

```python
ds = yt.load("data.0010.vtk", default_species_fields="ionized")
```

### Conventions on field names

Field names of on-disk fields for density, pressure, velocity and magnetic field components are always normalized to upper case, even if Pluto may use lowercase in some versions.

```python
>>> ds.field_list
[('pluto-vtk', 'PRS'),
 ('pluto-vtk', 'RHO'),
 ('pluto-vtk', 'VX1'),
 ('pluto-vtk', 'VX2'),
 ('pluto-vtk', 'VX3')]
```

This normalization is only applied to non-user-defined outputs and Pluto's ion
fraction outputs.
