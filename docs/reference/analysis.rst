.. _Analysis api:

=========
Analysis
=========

The ``Analysis`` class holds 'good records', and some common methods to analyze the observations.

.. currentmodule:: metobs_toolkit.analysis

Constructor
-----------

.. autosummary::
   :toctree: api/

   Analysis


Data attributes
----------------
A summary of all the attributes that hold or return data.

.. autosummary::
   :toctree: api/

   Analysis.df


Methods
----------
A summary of all methods in the ``Analysis`` class.

.. autosummary::
   :toctree: api/

   Analysis.get_info
   Analysis.get_tz
   Analysis.apply_filter_on_metadata
   Analysis.apply_filter_on_records
   Analysis.subset_period
   Analysis.aggregate_df
   Analysis.plot_diurnal_cycle
   Analysis.plot_diurnal_cycle_with_reference_station
