.. _Obstype api:

=========
Obstype
=========
The Obstype class represents an observation type, including its standard unit and description. It is used throughout the toolkit to define the nature and units of observational data.

A regular user typically interacts with Obstype instances via higher-level classes such as Dataset or Station, but direct use is possible for advanced customization.

.. currentmodule:: metobs_toolkit.obstypes

Constructor
-----------

.. autosummary::
   :toctree: api/

   Obstype

Attributes
----------------
A summary of all the attributes (and properties) of the Obstype class.

.. autosummary::
   :toctree: api/

   Obstype.name
   Obstype.std_unit
   Obstype.description
   Obstype.original_name
   Obstype.original_unit

Methods
----------
A summary of all methods in the Obstype class.

.. autosummary::
   :toctree: api/

   Obstype.get_compatible_units
   Obstype.is_compatible_with
   Obstype.get_info
   Obstype.convert_to_standard_units
