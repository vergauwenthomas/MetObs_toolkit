.. _Obstype api:

=========
Obstype
=========

Parent class defenition for observation types.

.. currentmodule:: metobs_toolkit

Constructor for Obstype
------------------------

.. autosummary::
   :toctree: api/

   Obstype

General methods and attributes for Obstype
-------------------------------------------

.. autosummary::
   :toctree: api/

   Obstype.set_description
   Obstype.set_original_name
   Obstype.set_original_unit
   Obstype.add_unit
   Obstype.get_orig_name
   Obstype.get_description
   Obstype.get_all_units
   Obstype.get_standard_unit
   Obstype.get_info



Developers methods and attributes for Obstype
----------------------------------------------

.. autosummary::
   :toctree: api/

   Obstype.convert_to_standard_units
   Obstype.test_if_unit_is_known

.. _ModelObstype api:

==============================
Obstypes used for GEE datasets
===============================

A child of :ref:`Obstype api` that adds info on how this observationtype is represented in a GEE dataset.

There are two classes:
 * `ModelObstype` : Represent a scalar Obstype for which there exists a band in a GEE dataset.
 * `ModelObstype_Vectorfield` : Represent a vectorfield, for which the *u* and *v* components exists in bands of a GEE dataset.


ModelObstype
-----------------------------

All methods of `Obstype` are inhereted.

.. autosummary::
   :toctree: api/

   ModelObstype
   ModelObstype.get_modelunit
   ModelObstype.get_modelband


ModelObstype_Vectorfield
-----------------------------

All methods of `Obstype` are inhereted.

.. autosummary::
   :toctree: api/

   ModelObstype_Vectorfield
   ModelObstype_Vectorfield.get_modelunit
   ModelObstype_Vectorfield.get_modelband_u
   ModelObstype_Vectorfield.get_modelband_v
