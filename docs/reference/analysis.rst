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

   analysis.Analysis

.. note::
   It is most common to construct the Analysis directly from a Dataset using `metobs_toolkit.Dataset.get_analysis()` method.

General methods and attributes
------------------------------

.. autosummary::
   :toctree: api/

   analysis.Analysis.subset_period
   analysis.Analysis.apply_filter
   analysis.Analysis.aggregate_df
   analysis.Analysis.get_anual_statistics
   analysis.Analysis.get_diurnal_statistics
   analysis.Analysis.get_diurnal_statistics_with_reference
   analysis.Analysis.get_aggregated_cycle_statistics
   analysis.Analysis.get_lc_correlation_matrices



Plotting methods
------------------------------

.. autosummary::
   :toctree: api/

   analysis.Analysis.plot_correlation_heatmap
   analysis.Analysis.plot_correlation_variation
