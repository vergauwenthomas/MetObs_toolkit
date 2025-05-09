.. _Station api:

=========
Station
=========

The Station holds all the data related to a single station. This includeds:
 * observational data, stored as ``SensorData``,
 * metadata, stored as a ``Site``,
 * and timeseries of a model source, stored as ``ModelTimeSeries``.

.. currentmodule:: metobs_toolkit.station

Constructor
-----------

.. autosummary::
   :toctree: api/

   Station

Data attributes
----------------
A summary of all the attributes that hold or return data.

.. autosummary::
   :toctree: api/

   Station.sensordata
   Station.df
   Station.metadf
   Station.outliersdf
   Station.gapsdf
   Station.modeldatadf
   Station.present_observations


General methods and attributes
------------------------------

.. autosummary::
   :toctree: api/


   Station.get_info
   Station.resample



GEE related methods
------------------------------
.. autosummary::
   :toctree: api/

   Station.modeldatadf
   Station.get_static_gee_point_data
   Station.get_static_gee_buffer_fraction_data
   Station.get_LCZ
   Station.get_altitude
   Station.get_landcover_fractions
   Station.get_gee_timeseries_data


Gaps related methods
------------------------------
.. autosummary::
   :toctree: api/

   Station.gapsdf
   Station.convert_outliers_to_gaps
   Station.interpolate_gaps
   Station.fill_gaps_with_raw_modeldata
   Station.fill_gaps_with_debiased_modeldata
   Station.fill_gaps_with_diurnal_debiased_modeldata
   Station.fill_gaps_with_weighted_diurnal_debiased_modeldata



Quality control related methods
--------------------------------
.. autosummary::
   :toctree: api/

   Station.outliersdf
   Station.gross_value_check
   Station.persistence_check
   Station.repetitions_check
   Station.step_check
   Station.window_variation_check
   Station.get_qc_stats
   Station.convert_outliers_to_gaps



Visualisations
------------------
.. autosummary::
   :toctree: api/

   Station.make_plot_of_modeldata
   Station.make_plot
