.. _Dataset api:

=========
Dataset
=========

The dataset is the heart of the MetObs Toolkit. It holds the observations and corresponding methods.

.. currentmodule:: metobs_toolkit

Constructor
-----------

.. autosummary::
   :toctree: api/

   Dataset

General methods and attributes
------------------------------

.. autosummary::
   :toctree: api/

   Dataset.show
   Dataset.get_info
   Dataset.show_settings
   Dataset.get_gaps_df
   Dataset.get_gaps_info
   Dataset.get_missing_obs_info
   Dataset.get_full_status_df

Common methods
--------------------

.. autosummary::
   :toctree: api/

   Dataset.get_station
   Dataset.update_gaps_and_missing_from_outliers
   Dataset.fill_gaps_linear
   Dataset.fill_gaps_era5
   Dataset.fill_gaps_automatic
   Dataset.fill_missing_obs_linear
   Dataset.get_analysis
   Dataset.apply_quality_control
   Dataset.apply_buddy_check
   Dataset.apply_titan_buddy_check
   Dataset.apply_titan_sct_resistant_check

Gap related methods
-----------------------

.. autosummary::
   :toctree: api/

   Dataset.interpolate_gaps

Extracting data
------------------
.. autosummary::
   :toctree: api/

   Dataset.get_modeldata
   Dataset.get_lcz
   Dataset.get_altitude
   Dataset.get_landcover



Updating Settings
------------------------

.. autosummary::
   :toctree: api/

   Dataset.update_settings
   Dataset.add_new_observationtype
   Dataset.add_new_unit
   Dataset.update_timezone
   Dataset.update_default_name
   Dataset.update_gap_and_missing_fill_settings
   Dataset.update_qc_settings
   Dataset.update_titan_qc_settings


Reading and writing related
------------------------------

.. autosummary::
   :toctree: api/

   Dataset.import_data_from_file
   Dataset.save_dataset
   Dataset.import_dataset
   Dataset.write_to_csv
   Dataset.sync_observations
   Dataset.coarsen_time_resolution


Plotting methods
------------------------------

.. autosummary::
   :toctree: api/

   Dataset.make_plot
   Dataset.make_interactive_plot
   Dataset.make_geo_plot
   Dataset.make_gee_plot
   Dataset.get_qc_stats
