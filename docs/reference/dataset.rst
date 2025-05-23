.. _Dataset api:

=========
Dataset
=========

The Dataset holds a collection of Stations. All methods applied on a dataset are applied on all present stations (and a target observationtype).

.. currentmodule:: metobs_toolkit.dataset

Constructor
-----------

.. autosummary::
   :toctree: api/

   Dataset

Data attributes
----------------
A summary of all the attributes that hold or return data.

.. autosummary::
   :toctree: api/

   Dataset.stations
   Dataset.obstypes
   Dataset.template
   Dataset.df
   Dataset.metadf
   Dataset.outliersdf
   Dataset.gapsdf
   Dataset.modeldatadf
   Dataset.present_observations


General methods and attributes
------------------------------

.. autosummary::
   :toctree: api/


   Dataset.get_info
   Dataset.subset_by_stations
   Dataset.get_station
   Dataset.rename_stations

   Dataset.sync_records
   Dataset.resample


   Dataset.add_new_observationtype
   Dataset.save_dataset_to_pkl
   Dataset.import_data_from_file


GEE related methods
------------------------------
.. autosummary::
   :toctree: api/

   Dataset.modeldatadf
   Dataset.import_gee_data_from_file
   Dataset.get_static_gee_point_data
   Dataset.get_static_gee_buffer_fraction_data
   Dataset.get_LCZ
   Dataset.get_altitude
   Dataset.get_landcover_fractions
   Dataset.get_gee_timeseries_data


Gaps related methods
------------------------------
.. autosummary::
   :toctree: api/

   Dataset.gapsdf
   Dataset.convert_outliers_to_gaps
   Dataset.interpolate_gaps
   Dataset.fill_gaps_with_raw_modeldata
   Dataset.fill_gaps_with_debiased_modeldata
   Dataset.fill_gaps_with_diurnal_debiased_modeldata
   Dataset.fill_gaps_with_weighted_diurnal_debiased_modeldata



Quality control related methods
--------------------------------
.. autosummary::
   :toctree: api/

   Dataset.outliersdf
   Dataset.gross_value_check
   Dataset.persistence_check
   Dataset.repetitions_check
   Dataset.step_check
   Dataset.window_variation_check
   Dataset.buddy_check
   Dataset.get_qc_stats
   Dataset.convert_outliers_to_gaps



Visualisations
------------------
.. autosummary::
   :toctree: api/

   Dataset.make_plot_of_modeldata
   Dataset.make_plot
   Dataset.make_gee_plot


=============================
Dataset related functions
=============================

.. currentmodule:: metobs_toolkit.dataset

.. autosummary::
   :toctree: api/

   import_dataset_from_pkl
