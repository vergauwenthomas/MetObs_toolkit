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

   metobs_toolkit.analysis.Analysis

.. note::
   It is most common to construct the Analysis directly from a Dataset using `metobs_toolkit.Dataset.get_analysis()` method.

General methods and attributes
------------------------------

.. autosummary::
   :toctree: api/

   metobs_toolkit.analysis.Analysis.subset_period
   metobs_toolkit.analysis.Analysis.apply_filter
   metobs_toolkit.analysis.Analysis.aggregate_df
   metobs_toolkit.analysis.Analysis.get_anual_statistics
   metobs_toolkit.analysis.Analysis.get_diurnal_statistics
   metobs_toolkit.analysis.Analysis.get_diurnal_statistics_with_reference
   metobs_toolkit.analysis.Analysis.get_aggregated_cycle_statistics
   metobs_toolkit.analysis.Analysis.get_lc_correlation_matrices



Plotting methods
------------------------------

.. autosummary::
   :toctree: api/

   metobs_toolkit.analysis.Analysis.plot_correlation_heatmap
   metobs_toolkit.analysis.Analysis.plot_correlation_variation
