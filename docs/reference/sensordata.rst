.. _SensorData api:

===========
SensorData
===========
The SensorData class stores timeseries data from sensors at a station. It is used for handling, processing, and analyzing observational data.

SensorData objects are usually accessed via the .sensordata attribute of Station or Dataset objects.

.. currentmodule:: metobs_toolkit.sensordata

Constructor
-----------

.. autosummary::
   :toctree: api/

   SensorData

Attributes
----------
A summary of all the attributes (and properties) of the SensorData class.

.. autosummary::
   :toctree: api/

   SensorData.df
   SensorData.stationname
   SensorData.tz
   SensorData.start_datetime
   SensorData.end_datetime
   SensorData.freq

Methods
-------
A summary of all methods in the SensorData class.

.. autosummary::
   :toctree: api/

   SensorData.get_info
   SensorData.convert_outliers_to_gaps
   SensorData.resample
   SensorData.get_qc_freq_statistics
   SensorData.fill_gap_with_modeldata
   SensorData.interpolate_gaps
   SensorData.convert_to_standard_units
   SensorData.gross_value_check
   SensorData.persistence_check
   SensorData.repetitions_check
   SensorData.step_check
   SensorData.window_variation_check
