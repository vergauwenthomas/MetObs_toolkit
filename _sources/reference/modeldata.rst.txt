.. _Modeldata api:

============
Modeldata
============

The Modeldata holds timeseries data which is the output of a model and methods
to interact with the Google Earth Engine to extract timeseries.

.. note::
   Keep in mind that the Modeldata does not hold gridded data, but timeseries of
   point extractions at the locations of the stations. Therefore a Modeldata
   instance is linked to a collection of stations defined by the `metadf` attribute.

Demo examples on the Modeldata class can be found here: :ref:`Extracting ERA5 timeseries`_ .


.. currentmodule:: metobs_toolkit

Constructor
-----------

.. autosummary::
   :toctree: api/

   Modeldata



Common methods
--------------------

.. autosummary::
   :toctree: api/

   Modeldata.add_obstype
   Modeldata.save_modeldata
   Modeldata.import_modeldata
   Modeldata.make_plot


GEE interactions
------------------
.. autosummary::
   :toctree: api/

   Modeldata.add_gee_dataset
   Modeldata.list_gee_datasets
   Modeldata.get_gee_dataset_data
   Modeldata.get_ERA5_data
   Modeldata.set_model_from_csv
