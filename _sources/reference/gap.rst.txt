.. _Gap api:

=========
Gap
=========

The class defenition of a gap.

.. note::
    For most usecases, there is no need to interact or initiate a `Gap`. For
    all relevant methods of `Gap` there is an equivalent method in `Dataset` but
    then applied on all gaps in the `Dataset`.

    See the :ref:`Dataset <Dataset_gaps api>` API.


.. currentmodule:: metobs_toolkit

API
------------------------

.. autosummary::
   :toctree: api/

   Gap
   Gap.get_info
   Gap.interpolate
   Gap.raw_model_gapfill
   Gap.debias_model_gapfill
   Gap.diurnal_debias_model_gapfill
   Gap.weighted_diurnal_debias_model_gapfill
