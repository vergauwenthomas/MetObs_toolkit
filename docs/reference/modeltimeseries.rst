.. _ModelTimeSeries api:

================
ModelTimeSeries
================
The ModelTimeSeries class stores timeseries data extracted from a dynamic GEE dataset (e.g., ERA5) for a specific observation type at a station. It is similar to the SensorData class and is typically used for comparison with observations, quality control, or gap filling.

ModelTimeSeries objects are usually accessed via the .modeldata attribute of Station or Dataset objects.

.. currentmodule:: metobs_toolkit.modeltimeseries

Constructor
-----------

.. autosummary::
   :toctree: api/

   ModelTimeSeries

Attributes
----------
A summary of all the attributes (and properties) of the ModelTimeSeries class.

.. autosummary::
   :toctree: api/

   ModelTimeSeries.df
   ModelTimeSeries.stationname
   ModelTimeSeries.tz
   ModelTimeSeries.start_datetime
   ModelTimeSeries.end_datetime
   ModelTimeSeries.freq

Methods
-------
A summary of all methods in the ModelTimeSeries class.

.. autosummary::
   :toctree: api/

   ModelTimeSeries.get_info
   ModelTimeSeries.make_plot
