# MetObs-toolkit

[![PyPI version](https://badge.fury.io/py/metobs-toolkit.svg)](https://badge.fury.io/py/metobs-toolkit)
[![Documentation Status](https://readthedocs.org/projects/metobs-toolkit/badge/?version=latest)](https://metobs-toolkit.readthedocs.io/en/latest/?badge=latest)
[![status](https://joss.theoj.org/papers/ffa3a79315bdf4c4793992a1de41193d/status.svg)](https://joss.theoj.org/papers/ffa3a79315bdf4c4793992a1de41193d)
[![Tests passing](https://github.com/vergauwenthomas/MetObs_toolkit/actions/workflows/main_workflow.yml/badge.svg?branch=master)](https://github.com/vergauwenthomas/MetObs_toolkit/actions/workflows/main_workflow.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10794417.svg)](https://doi.org/10.5281/zenodo.10794417)

<img src="https://raw.githubusercontent.com/vergauwenthomas/MetObs_toolkit/master/docs/logo_small.jpeg" alt="drawing" style="width:200px;"/>



The MetObs-toolkit provides a comprehensive framework for scientists to process raw meteorological data for analysis.

This repo contains all the software for the [metobs_toolkit](https://test.pypi.org/project/metobs-toolkit/).

## Documentation

Documentation can be found [here](https://metobs-toolkit.readthedocs.io/en/latest/index.html).

The Documentation contains examples, explanations, and the API documentation of
the toolkit. Make sure that the documentation version matches your version
of the toolkit.


## Installing the package
Install the package using pip:

```bash
pip install metobs-toolkit
 ```
To install the PyPi version of the toolkit. To install the github versions one can use these commands:

```bash
#main versions
pip install git+https://github.com/vergauwenthomas/MetObs_toolkit.git

#development version
pip install git+https://github.com/vergauwenthomas/MetObs_toolkit.git@dev

#specific release from github
pip install git+https://github.com/vergauwenthomas/MetObs_toolkit.git@v0.2.0
 ```
For some advanced quality control methods, the [Titanlib package](https://github.com/metno/titanlib) is used. Since the installation of titanlib requires a c++ compiler, we have chosen not to include it in the toolkit. If you want to use the Titanlib functionality you must install both the toolkit and Titanlib:

```bash
pip install metobs-toolkit titanlib
 ```

To use the package, import it in Python:


```python
import metobs_toolkit

#Check your version
metobs_toolkit.__version__
 ```
## Exercises and demos
In the context of a [FAIRNESS (COST action)](https://www.fairness-ca20108.eu/) summer school, a set of well-documented exercises and demos are made.

| Notebook  | Description  |       |
|:----------|:-------------|------:|
| [Introduction](https://github.com/vergauwenthomas/MetObs_toolkit/blob/master/fairness_demo_exercises/Introduction_01.ipynb) | Introduction to the toolkit | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/vergauwenthomas/MetObs_toolkit/blob/master/fairness_demo_exercises/Introduction_01.ipynb) |
| [Quality control](https://github.com/vergauwenthomas/MetObs_toolkit/blob/master/fairness_demo_exercises/Quality_control_excercise_02.ipynb) | Introduction to quality control methods | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/vergauwenthomas/MetObs_toolkit/blob/master/fairness_demo_exercises/Quality_control_excercise_02.ipynb)|
| [Filling gaps](https://github.com/vergauwenthomas/MetObs_toolkit/blob/master/fairness_demo_exercises/Gap_filling_excercise_03.ipynb) | Introduction to gap filling methods | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/vergauwenthomas/MetObs_toolkit/blob/master/fairness_demo_exercises/Gap_filling_excercise_03.ipynb)|
| [Analysis](https://github.com/vergauwenthomas/MetObs_toolkit/blob/master/fairness_demo_exercises/Urban_analysis_excercise_04.ipynb) | Introduction analysis methods | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/vergauwenthomas/MetObs_toolkit/blob/master/fairness_demo_exercises/Urban_analysis_excercise_04.ipynb)|

*NOTE: * these exercises are built on an early version of the toolkit! We recommend using the examples in the [documentation](https://metobs-toolkit.readthedocs.io/en/latest/index.html).



## Related
* A graphical user interface for the MetObs-Toolkit: [MetObs_GUI](https://github.com/vergauwenthomas/MetObs_GUI)
* A JOSS (Journal of Open Source Software) [publication on the MetObs-Toolkit](https://joss.theoj.org/papers/10.21105/joss.05916#)
