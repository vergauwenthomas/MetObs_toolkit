.. _Analysis api:

=========
Analysis
=========

The Analysis class holds observations (assumed to be correct) and methods for
analysing and aggregating the data. An Analysis is typical created from a
Dataset after the cleanup of the raw observations.

.. currentmodule:: metobs_toolkit

Constructor
-----------

.. autosummary::
   :toctree: api/

   Analysis


Data and filtering methods
----------------------------

.. autosummary::
   :toctree: api/


   Analysis.subset_period
   Analysis.get_possible_filter_keywords
   Analysis.apply_filter

Aggregation methods
---------------------

.. autosummary::
   :toctree: api/

   Analysis.get_analysis_records
   Analysis.aggregate_df
   Analysis.get_anual_statistics
   Analysis.get_diurnal_statistics
   Analysis.get_diurnal_statistics_with_reference
   Analysis.get_aggregated_cycle_statistics


Other methods
---------------

.. autosummary::
   :toctree: api/

   Analysis.get_lc_correlation_matrices
   Analysis.plot_correlation_heatmap
   Analysis.plot_correlation_variation
