:html_theme.sidebar_secondary.remove:

MetObs Toolkit
============================
.. image:: logo_small.svg
    :alt: logo
    :width: 500

MetObs Toolkit is an open source project to make working with meteorological observations in python easier.


Description
-----------

The MetObs-toolkit is a Python package developed to facilitate the use of non-traditional meteorological observations.
The package provides automated quality control (QC) techniques to identify and
flag erroneous observations, and includes methods to fill data gaps.
Additionally, the package offers tools for analyzing the data, e.g. linkage
with popular land-use datasets is included, through the use of the Google Earth Engine, such
that microclimate effects can be investigated with the MetObs-toolkit.

The MetObs-toolkit provides a comprehensive framework for scientists to process
raw meteorological data for analysis by making intensive use of the `Pandas <https://pandas.pydata.org/>`_
and `GeoPandas <https://geopandas.org/en/stable/>`_ functionalities.

.. toctree::
   :hidden:
   :maxdepth: 5

   Home <self>
   Introduction <intro>
   Examples <examples/index>
   Documentation <reference/index>
   Referencing <paper/index>
   Extra topics <topics/index>
   Contributing <contributing_link>


How to install
=================

To use the package python 3.10 or higher is required.
To install the package one can use pip:

.. code-block:: console

   pip install metobs-toolkit

Or install the packge using `conda`:

.. code-block:: console

   conda install -c conda-forge metobs-toolkit


To install the latest development version: 

.. code-block:: console

   #development version
   pip install git+https://github.com/vergauwenthomas/MetObs_toolkit.git@dev


.. note::
   To install the package in a notebook, one has to add ! in front of the pip install command.


.. code-block:: python

   import metobs_toolkit

   #Check your version
   metobs_toolkit.__version__



Indices and tables
----------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
