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
   Dataset.import_data_from_file
   Dataset.save_dataset
   import_dataset
   Dataset.get_full_status_df
   Dataset.add_new_observationtype
   Dataset.add_new_unit
   Dataset.show_settings
   Dataset.get_station
   Dataset.write_to_csv
   Dataset.coarsen_time_resolution
   Dataset.sync_records


Updating Settings
-------------------

.. autosummary::
   :toctree: api/

   Dataset.update_file_paths
   Dataset.update_output_dir
   Dataset.update_default_name
   Dataset.update_qc_settings
   Dataset.update_titan_qc_settings



Gaps related methods
------------------------------
.. _Dataset_gaps api:

.. autosummary::
   :toctree: api/

   Dataset.get_gaps_fill_df
   Dataset.convert_outliers_to_gaps
   Dataset.find_gap
   Dataset.interpolate_gaps
   Dataset.fill_gaps_with_raw_modeldata
   Dataset.fill_gaps_with_debiased_modeldata
   Dataset.fill_gaps_with_diurnal_debiased_modeldata
   Dataset.fill_gaps_with_weighted_diurnal_debias_modeldata



Quality control related methods
--------------------------------

.. autosummary::
   :toctree: api/

   Dataset.get_qc_stats
   Dataset.apply_quality_control
   Dataset.apply_buddy_check
   Dataset.apply_titan_buddy_check



Visualisations
------------------

.. autosummary::
   :toctree: api/

   Dataset.make_plot
   Dataset.make_geo_plot
   Dataset.make_interactive_plot
   Dataset.make_gee_static_spatialplot
   Dataset.make_gee_dynamic_spatialplot

Extracting data
------------------

.. autosummary::
   :toctree: api/

   Dataset.get_modeldata
   Dataset.add_new_geemodeldata
   Dataset.get_lcz
   Dataset.get_altitude
   Dataset.get_landcover


Special methods
------------------

.. autosummary::
   :toctree: api/

   Dataset.__add__
