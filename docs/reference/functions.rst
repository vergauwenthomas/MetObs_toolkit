.. _functions api:

==================================
MetObs toolkit defaults
==================================

In the MetObs toolkit, default variables are defined.
These variables can be handy to use for some applications.


.. code-block:: python

   import metobs_toolkit

   #The default GEE datasetmanagers
   metobs_toolkit.default_GEE_datasets


==================================
MetObs toolkit Utility functions
==================================

A collection of utility functions is presented. These functions are not
part of a class, but are used in the toolkit.



Logging functions
------------------
Utility functions related to logging.

.. autosummary::
   :toctree: api/

   metobs_toolkit.add_FileHandler
   metobs_toolkit.add_StreamHandler

Authentication functions
-------------------------------
Utility functions related to authentication.

.. autosummary::
   :toctree: api/

   metobs_toolkit.connect_to_gee


Other utility functions
-------------------------------

.. autosummary::
   :toctree: api/

   metobs_toolkit.import_dataset_from_pkl
