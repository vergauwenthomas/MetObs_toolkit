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

   Obstype.get_plot_y_label
   Obstype.convert_to_standard_units
   Obstype.test_if_unit_is_known

.. _ModelObstype api:

==============
ModelObstype
==============

A child of :ref:`Obstype api` that adds info on how this observationtype is represented in modeloutput.
All methods of Obstype() are inhereted.

Constructor for ModelObstype
-----------------------------

.. autosummary::
   :toctree: api/

   ModelObstype

General methods and attributes for ModelObstype
------------------------------------------------

.. autosummary::
   :toctree: api/

   ModelObstype.get_mapped_datasets
   ModelObstype.get_bandname
   ModelObstype.get_modelunit
   ModelObstype.add_new_band

Developers methods and attributes for ModelObstype
---------------------------------------------------

.. autosummary::
   :toctree: api/

   ModelObstype.get_bandname_mapper
   ModelObstype.has_mapped_band


.. _ModelObstype_Vectorfield api:

==========================
ModelObstype_Vectorfield
==========================

A child of :ref:`Obstype api`, similar to :ref:`ModelObstype api`, that adds info on how this handle 2D vectorfields in the modeloutput.
A vectorfield in the modeloutput is defined by its components.

All methods of Obstype() are inhereted.

Constructor for ModelObstype_Vectorfield
-----------------------------------------

.. autosummary::
   :toctree: api/

   ModelObstype_Vectorfield

General methods and attributes for ModelObstype_Vectorfield
-------------------------------------------------------------

.. autosummary::
   :toctree: api/

   ModelObstype_Vectorfield.get_mapped_datasets
   ModelObstype_Vectorfield.get_modelunit
   ModelObstype_Vectorfield.add_new_band

Developers methods and attributes
------------------------------------

.. autosummary::
   :toctree: api/

   ModelObstype_Vectorfield.get_bandname_mapper
   ModelObstype_Vectorfield.has_mapped_band


Developers vectorfield conversion functions
---------------------------------------------
These functions are used by the :ref:`ModelObstype_Vectorfield api` to convert components to amplitudes and angles.

.. autosummary::
   :toctree: api/

   obstype_modeldata.compute_amplitude
   obstype_modeldata.compute_angle
