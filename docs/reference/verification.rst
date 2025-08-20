.. _Verification api:

=============
Verification
=============

The Verification class holds all the data related to a single verification process. This includes:

 - verification dataframe where observations are compared against model values
 - metadata
 - Obstypes
 
.. currentmodule:: metobs_toolkit.verification

Constructor
-----------

.. autosummary::
   :toctree: api/

   Verification

Data attributes
----------------
A summary of all the attributes that hold or return data.

.. autosummary::
   :toctree: api/

   Verification.verifdf
   Verification.metadf
  

General methods and attributes
------------------------------

.. autosummary::
   :toctree: api/

   Verification.save_verification_to_pkl
   Verification.traditional_scores
   Verification.group_values
   Verification.diurnal_cycle_values

Visualisations
------------------

.. autosummary::
   :toctree: api/

   Verification.scatter_plot
   Verification.plot_diurnal_cycle
   Verification.timeseries_plot
   Verification.scores_plot


=============================
Verification related functions
=============================

.. currentmodule:: metobs_toolkit.verification

.. autosummary::
   :toctree: api/

   import_verification_from_pkl


