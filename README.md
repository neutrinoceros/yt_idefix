
# yt_idefix
[![PyPI](https://img.shields.io/pypi/v/yt_idefix)](https://pypi.org/project/yt_idefix)
[![Supported Python Versions](https://img.shields.io/pypi/pyversions/yt_idefix/0.1.0)](https://pypi.org/project/yt_idefix/)
[![yt-project](https://img.shields.io/static/v1?label="works%20with"&message="yt"&color="blueviolet")](https://yt-project.org)

<!--- Tests and style --->
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/neutrinoceros/yt_idefix/main.svg)](https://results.pre-commit.ci/latest/github/neutrinoceros/yt_idefix/main)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Imports: isort](https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat&labelColor=ef8336)](https://pycqa.github.io/isort/)

A maturing yt frontend for Idefix, packaged as an extension for yt.
This frontend is a candidate for integration in the core yt code base.

## Installation

Make sure you have Python 3.6 or newer, then run
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
Then you should be able to load Idefix data seamlessly with `yt.load`.


## Current limitations

As of version 0.1.0 of this project, the frontend allows one to read Idefix's
dumpfiles only, through the `IdefixDumpDataset` class. `IdefixVTKDataset` may be
implemented in the future, but won't be usable directly with `yt.load` before
the next yt bugfix release is available (see bellow).

The `IdefixDumpDataset` class is functional but far from optimized, it may take
much longer than strictly needed to perform queries. This will be adressed in
the future.

As of yt 4.0.1:
- Non-uniform grids (using log spacing) are not supported, which makes this
  frontend of very limited use for Idefix.
- `yt.load()` is not suitable for vtk files that are not produced by Athena.
  https://github.com/yt-project/yt/issues/3001 This bug is however resolved on
  the main branch (https://github.com/yt-project/yt/pull/3424)
