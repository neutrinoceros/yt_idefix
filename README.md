
# yt_idefix
[![PyPI](https://img.shields.io/pypi/v/yt-idefix.svg?logo=pypi&logoColor=white&label=PyPI)](https://pypi.org/project/yt_idefix/)
[![PyPI](https://img.shields.io/badge/requires-Python%20≥%203.8-blue?logo=python&logoColor=white)](https://pypi.org/project/yt_idefix/)
[![yt-project](https://img.shields.io/static/v1?label="works%20with"&message="yt"&color="blueviolet")](https://yt-project.org)

<!--- Tests and style --->
[![CI](https://github.com/neutrinoceros/yt_idefix/actions/workflows/ci.yml/badge.svg)](https://github.com/neutrinoceros/yt_idefix/actions/workflows/ci.yml)
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/neutrinoceros/yt_idefix/main.svg)](https://results.pre-commit.ci/latest/github/neutrinoceros/yt_idefix/main)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Imports: isort](https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat&labelColor=ef8336)](https://pycqa.github.io/isort/)

A maturing yt frontend for Idefix and Pluto (vtk), packaged as an extension for yt.
This frontend is a candidate for integration in the core yt code base.

## Installation

```shell
pip install yt_idefix
```
## Usage

After importing `yt` itself, make sure to activate the extension
```python
import yt
import yt_idefix
```

Now `yt.load` will be able to read Pluto/Idefix `.vtk` output files, as well as Idefix `.dmp` dump files.


## Strecthed grids support

### yt_idefix 0.12.0 and newer

version 0.12.0 brings experimental *native* support for streched grids, which is under
active development upstream, in yt itself.

Slices should now work seamlessly even with older versions of yt, however
yt 4.1 will be required to perform projections correctly.

**update**: yt_idefix 0.13.3 is the last release allowing yt 4.0.x, yt_idefix 0.14.0 requires yt 4.1


### yt_idefix 0.11 and older (deprecated)

yt_idefix ships a specialized loader function for datasets with streched grids
`yt_idefix.load_stretched`. This function is only provided as a workaround 4.0.x limitations, but it's highly limited itself:
- no field magic (no aliasing, or dimensionalization, or automatic derived field generation)
- no lazy loading (all data has to reside in memory)
- projections are not supported
- only supports vtk outputs (not dumps)
