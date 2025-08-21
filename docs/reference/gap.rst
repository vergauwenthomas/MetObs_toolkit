.. _Gap api:

=========
Gap
=========
A gap is a part of the observational timeseries that is missing. It is thus related to a station and an observationtype.
The gaps are described by the ``Gap`` class, containing methods for manipulating the gap (i.g. filling a gap).


A regular user *should not directly interact with a ``Gap`` instance*. All methods
for filling gaps are accessible in the ``Station`` and ``Dataset`` classes.


.. currentmodule:: metobs_toolkit.gap

Constructor
-----------

.. autosummary::
   :toctree: api/

   Gap


Attributes
----------------
A summary of all the attributes (and properties) of the Gap class.

.. autosummary::
   :toctree: api/

   Gap.records
   Gap.obstype
   Gap.stationname
   Gap.fillsettings
   Gap.fillstatus
   Gap.start_datetime
   Gap.end_datetime
   Gap.df

Methods
----------
A summary of all methods in the Gap class.

.. autosummary::
   :toctree: api/

   Gap.get_info
   Gap.flag_can_be_filled
   Gap.flush_fill
   Gap.interpolate
   Gap.raw_model_gapfill
   Gap.debiased_model_gapfill
   Gap.diurnal_debiased_model_gapfill
   Gap.weighted_diurnal_debiased_model_gapfill
