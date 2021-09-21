
# yt_idefix
[![PyPI](https://img.shields.io/pypi/v/yt_idefix)](https://pypi.org/project/yt_idefix)
[![Supported Python Versions](https://img.shields.io/pypi/pyversions/yt_idefix/0.4.0)](https://pypi.org/project/yt_idefix/)
[![yt-project](https://img.shields.io/static/v1?label="works%20with"&message="yt"&color="blueviolet")](https://yt-project.org)

<!--- Tests and style --->
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/neutrinoceros/yt_idefix/main.svg)](https://results.pre-commit.ci/latest/github/neutrinoceros/yt_idefix/main)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Imports: isort](https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat&labelColor=ef8336)](https://pycqa.github.io/isort/)

A maturing yt frontend for Idefix, packaged as an extension for yt.
This frontend is a candidate for integration in the core yt code base.

## Installation

Make sure you have Python 3.7 or newer, then run
```shell
python3 -m pip install yt_idefix
```
If you don't already have yt it will be installed along.

## Usage

After importing `yt` itself, make sure to activate the extension
```python
import yt
import yt.extensions.idefix
```
Single dump files as well as time series can be loaded directly with `yt.load`, e.g.,
```python
ds = yt.load("dump.0054.dmp")
ts = yt.load("dump.00??.dmp")
```

But vtk files currently require a little additional work
```python
# load a single dataset
from yt_idefix.api import IdefixVtkDataset

ds = IdefixVtkDataset("data.0042.vtk")

# load time series
class IdefixVtkDatasetSeries(yt.DatasetSeries):
    _dataset_cls = IdefixVtkDataset


ts = IdefixVtkDatasetSeries("data.00??.vtk")
```
This is because of a bug that will be fixed in yt's next release, see after.

## Current limitations

As of version 0.5.0 of this project, I/O performances are yet to be optimized
for both dump and vtk files.

Vtk support is limited to cartesian geometries and currently slightly hacky. It
will received better care in upcoming releases.

As of yt 4.0.1:
- Non-uniform grids (using log spacing) are not supported, which makes this
  frontend of very limited use for Idefix.
- `yt.load()` is not suitable for vtk files because they are considered ambiguous
  since *all* vtk files are (erroneously) recognized as Athena.
  This bug is resolved on yt's dev branch (https://github.com/yt-project/yt/pull/3424)
