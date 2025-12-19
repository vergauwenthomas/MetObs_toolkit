.. _Settings api:

=========
Settings
=========

The Settings class is a singleton configuration class that stores global settings 
used across multiple functions and methods in the MetObs toolkit. Settings stored 
here affect behavior throughout the toolkit (e.g., timezone handling, label 
definitions, logging).


The Settings class uses the singleton pattern to ensure only one instance exists.
Settings are accessed and modified using class methods, so no instantiation is 
required by the user. That makes the `Settings` especially convenient for 
global configuration, and settings used at the deep modules without having to 
pass all the settings as arguments (e.g. default style settings for plots).

.. note::
   The Settings class is accessible directly from the ``metobs_toolkit`` module
   as ``metobs_toolkit.Settings``.

.. currentmodule:: metobs_toolkit.settings_collection.settings

Constructor
-----------

.. autosummary::
   :toctree: api/

   Settings


Getting and Setting Values
---------------------------

Methods for retrieving and modifying configuration values.

.. autosummary::
   :toctree: api/

   Settings.get
   Settings.set
   Settings.reset


Inspection Methods
-------------------

Methods for viewing and exporting the current configuration.

.. autosummary::
   :toctree: api/

   Settings.get_info
   Settings.to_dict


Usage Examples
--------------

Getting settings:

.. code-block:: python

   import metobs_toolkit

   # Get a top-level setting
   timezone = metobs_toolkit.Settings.get("store_tz")  # Returns 'UTC'

   # Get a nested setting using dot notation
   label = metobs_toolkit.Settings.get("label_def.goodrecord.label")  # Returns 'ok'

   # Get with a default fallback
   value = metobs_toolkit.Settings.get("nonexistent_key", "fallback")

Modifying settings:

.. code-block:: python

   # Set a top-level setting
   metobs_toolkit.Settings.set("store_tz", "Europe/Brussels")

   # Set a nested setting
   metobs_toolkit.Settings.set("log_level", "DEBUG")

Resetting settings:

.. code-block:: python

   # Reset a specific setting to its default
   metobs_toolkit.Settings.reset("store_tz")

   # Reset all settings to defaults
   metobs_toolkit.Settings.reset()

Viewing current settings:

.. code-block:: python

   # Print all settings
   metobs_toolkit.Settings.get_info()

   # Get settings as a dictionary
   config = metobs_toolkit.Settings.to_dict()




