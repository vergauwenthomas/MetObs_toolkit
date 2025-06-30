.. _Geedatasetmanagers api:

=====================
GEEDatasetManagers
=====================

A GEEDatasetManager is a class that handles the interaction between metobs_toolkit and a specific
dataset on Google Earth Engine. No modeldata is stored in the managers.

There are two GEEDatasetManagers:
   * ``GEEStaticDatasetManager``: Reference to a dataset on GEE that does not have a time-dimension.
   * ``GEEDynamicDatasetManager``: Reference to a dataset on GEE with a time dimension. ``ModelObstype``s are used to map a dataset band to an equivalent ``Obstype``.

.. currentmodule:: metobs_toolkit.geedatasetmanagers

Constructors
-------------

.. autosummary::
   :toctree: api/

   GEEStaticDatasetManager
   GEEDynamicDatasetManager


Methods for Static dataset managers
------------------------------------
A summary of all methods in the GEEStaticDatasetManager class.

.. autosummary::
   :toctree: api/

   GEEStaticDatasetManager.get_info
   GEEStaticDatasetManager.extract_static_point_data
   GEEStaticDatasetManager.extract_static_buffer_frac_data
   GEEStaticDatasetManager.make_gee_plot


Methods for Dynamic dataset managers
----------------------------------------
A summary of all methods in the GEEDynamicDatasetManager class.

.. autosummary::
   :toctree: api/

   GEEDynamicDatasetManager.get_info
   GEEDynamicDatasetManager.add_modelobstype
   GEEDynamicDatasetManager.extract_timeseries_data
   GEEDynamicDatasetManager.make_gee_plot


Special methods
------------------

The DatasetManager classes implements several Python special methods for convenience:

- ``__str__`` and ``__repr__``: String representations for printing and debugging.
