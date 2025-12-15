.. _WhiteSet api:

=========
WhiteSet
=========

The WhiteSet class manages whitelisted records for quality control operations. It allows specific observations (identified by station name, observation type, and/or timestamp) to be excluded from outlier detection during QC checks. Whitelisted records participate in all QC calculations but are not flagged as outliers in the final results.

A regular user interacts with WhiteSet instances when calling QC methods on ``Dataset`` and ``Station`` classes, providing fine-grained control over which records should be protected from being flagged as outliers.

.. important::
   In practice, users only create WhiteSet instances to pass them as parameters to QC methods. A user does generally not need to interact with this class through its methods.

.. currentmodule:: metobs_toolkit.qc_collection.whitelist

Constructor
-----------

.. autosummary::
   :toctree: api/

   WhiteSet


Methods
----------
A summary of all methods in the WhiteSet class.

.. autosummary::
   :toctree: api/

   WhiteSet.get_info
   WhiteSet.create_sensorwhitelist



Index Structure
---------------

The ``white_records`` index can have one or more of the following levels:

* **name**: Station identifier (matches the station name)
* **obstype**: Observation type (e.g., 'temp', 'humidity')
* **datetime**: Specific timestamps to whitelist

If the 'datetime' level is absent, all timestamps for matching station/obstype combinations are whitelisted.

Timezone Handling
-----------------

.. important::
   When the ``white_records`` index contains a 'datetime' level, timestamps are automatically 
   handled during WhiteSet initialization:
   
   * **Timezone-aware timestamps** are converted to UTC
   * **Timezone-naive timestamps** are localized to UTC (with a warning)
   
   It is **strongly recommended** to provide timezone-aware timestamps to avoid ambiguity and 
   ensure correct matching during quality control operations.


Examples
--------

**Create a WhiteSet with datetime-only whitelisting:**

.. code-block:: python

   import pandas as pd
   import metobs_toolkit
   
   # Whitelist specific timestamps across all stations
   timestamps = pd.date_range('2022-09-01 00:00', periods=10, freq='1h', tz='UTC')
   whiteset = metobs_toolkit.WhiteSet(
       pd.Index(timestamps, name='datetime')
   )
   
   # Use in QC check
   dataset.gross_value_check(
       obstype='temp',
       lower_threshold=10.0,
       upper_threshold=25.0,
       whiteset=whiteset
   )


**Create a WhiteSet with station and datetime:**

.. code-block:: python

   # Whitelist specific timestamps for specific stations
   white_records = pd.MultiIndex.from_arrays([
       ['station1', 'station1', 'station2'],
       pd.to_datetime(['2022-09-01 12:00', '2022-09-01 13:00', '2022-09-01 14:00'])
   ], names=['name', 'datetime'])
   
   whiteset = metobs_toolkit.WhiteSet(white_records)
   
   dataset.persistence_check(
       obstype='temp',
       timewindow='2h',
       whiteset=whiteset
   )


**Create a WhiteSet with full specification:**

.. code-block:: python

   # Whitelist specific records with all levels
   white_records = pd.MultiIndex.from_arrays([
       ['station1', 'station1', 'station2'],
       ['temp', 'humidity', 'temp'],
       pd.to_datetime(['2022-09-01 12:00', '2022-09-01 13:00', '2022-09-01 14:00'])
   ], names=['name', 'obstype', 'datetime'])
   
   whiteset = metobs_toolkit.WhiteSet(white_records)



Related Classes
---------------

The WhiteSet uses the SensorWhiteSet class internally for station-specific and obstype-specific filtering:

.. currentmodule:: metobs_toolkit.qc_collection.whitelist

.. autosummary::
   :toctree: api/

   SensorWhiteSet


Notes
-----

* Whitelisted records participate in all QC calculations (e.g., mean, standard deviation in buddy checks) but are protected from being flagged as outliers.
* An empty WhiteSet (default) means no records are whitelisted.
* The WhiteSet is validated upon initialization to ensure it has the correct index structure.
* When using WhiteSet with Dataset-level QC methods, the whitelist is automatically filtered for each station and observation type combination.
