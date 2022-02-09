
# yt_idefix
[![PyPI](https://img.shields.io/pypi/v/yt-idefix.svg?logo=pypi&logoColor=white&label=PyPI)](https://pypi.org/project/yt_idefix/)
[![PyPI](https://img.shields.io/pypi/pyversions/yt-idefix/0.11.4?logo=python&logoColor=white&label=Python)](https://pypi.org/project/yt_idefix/)
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


### Strecthed grids support

yt_idefix ships a specialized loader function for datasets with streched grids
`yt_idefix.load_stretched`. This function is experimental. Here are its known
limitations.
- no field magic (no aliasing, or dimensionalization, or automatic derived field generation)
- no lazy loading (all data has to reside in memory)
- projections are not supported
- only supports vtk outputs (not dumps)
