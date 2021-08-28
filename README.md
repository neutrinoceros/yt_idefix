
# yt_idefix
[![PyPI](https://img.shields.io/pypi/v/yt_idefix)](https://pypi.org/project/yt_idefix)
[![Supported Python Versions](https://img.shields.io/pypi/pyversions/yt_idefix/0.1.0)](https://pypi.org/project/yt_idefix/)
[![yt-project](https://img.shields.io/static/v1?label="works%20with"&message="yt"&color="blueviolet")](https://yt-project.org)

<!--- Tests and style --->
[![pre-commit.ci status](https://results.pre-commit.ci/badge/github/neturinoceros/yt_idefix/main.svg)](https://results.pre-commit.ci/latest/github/neturinoceros/yt_idefix/main)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Imports: isort](https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat&labelColor=ef8336)](https://pycqa.github.io/isort/)

A maturing yt frontend for Idefix, packaged as an extension for yt.
This frontend is a candidate for integration in the core yt code base.

## Installation

Make sure you have Python 3.6 or newer, then run
```shell
python3 -m pip install yt_idefix
```


## Usage

On top of import `yt` itself, make sure to activate the extension as
```python
import yt
import yt.extensions.idefix
```
Then you should be able to load Idefix snapshots seamlessly with `yt.load`.
