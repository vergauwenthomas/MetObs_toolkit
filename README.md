# MetObs-toolkit

[![PyPI version](https://badge.fury.io/py/metobs-toolkit.svg)](https://badge.fury.io/py/metobs-toolkit)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/metobs-toolkit.svg)](https://anaconda.org/conda-forge/metobs-toolkit)
[![Documentation Status](https://readthedocs.org/projects/metobs-toolkit/badge/?version=latest)](https://metobs-toolkit.readthedocs.io/en/latest/?badge=latest)
[![status](https://joss.theoj.org/papers/ffa3a79315bdf4c4793992a1de41193d/status.svg)](https://joss.theoj.org/papers/ffa3a79315bdf4c4793992a1de41193d)
[![Tests passing](https://github.com/vergauwenthomas/MetObs_toolkit/actions/workflows/main_workflow.yml/badge.svg?branch=master)](https://github.com/vergauwenthomas/MetObs_toolkit/actions/workflows/main_workflow.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10794417.svg)](https://doi.org/10.5281/zenodo.10794417)

[<img src="https://raw.githubusercontent.com/vergauwenthomas/MetObs_toolkit/master/docs/logo_small.jpeg" alt="drawing" style="width:200px;"/>](https://metobs-toolkit.readthedocs.io/en/latest/index.html)

The MetObs-toolkit provides a comprehensive framework for scientists to process, quality control, and analyze raw meteorological data. It is designed to be flexible, extensible, and user-friendly for a wide range of meteorological applications.

## Documentation

Full documentation, including installation instructions, usage examples, and API reference, is available at:

👉 [https://metobs-toolkit.readthedocs.io/en/latest/index.html](https://metobs-toolkit.readthedocs.io/en/latest/index.html)

Please ensure the documentation version matches your installed version of the toolkit.

## Installation

Install the latest release from PyPI:

```bash
pip install metobs-toolkit
```
Install using `conda`

```bash
conda install -c conda-forge metobs-toolkit 
```
or using `mamba`

```bash
mamba install metobs-toolkit
```

To install the latest development version from GitHub:

```bash
pip install git+https://github.com/vergauwenthomas/MetObs_toolkit.git@dev
```

## Usage

Import the package in Python:

```python
import metobs_toolkit

# Check your version
print(metobs_toolkit.__version__)
```

## Related Projects

* [MetObs_GUI](https://github.com/vergauwenthomas/MetObs_GUI): A graphical user interface for the MetObs-toolkit
* [JOSS publication on the MetObs-toolkit](https://joss.theoj.org/papers/10.21105/joss.05916#)
