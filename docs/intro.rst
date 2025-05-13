**************************
Toolkit objects overview
**************************

This toolkit is a Python package based on object-oriented programming (OOP).
Here you can find a short description of the classes that are commenly used.
In the `introduction example` you can find a notebook example to illustrate these classes.


``Dataset``
------------

The ``Datset`` class is at the heart of the toolkit and it holds a set of ``Station`` s.
All methods applied on a ``Dataset`` are equally applied on all its ``Station`` s. Thus all methods that can be applied
on a ``Station`` can be applied on a ``Dataset``. See the API documentation for more details :ref:`Dataset API <Dataset api>` and :ref:`Station API <Station api>`.

.. code-block:: python

   your_dataset = metobs_toolkit.Dataset()

The dataset holds methods for
 - Importing raw data
 - Resampling/syncronizing timeseries
 - Extracting metadata
 - Visualizing
 - Quality control
 - Gap filling



``Station``
-----------

A ``Station`` holds:
- timeseries related to sensors (as `SensorData`)
- metadata related to the site of the station (as a `Site`)
- timeseries originating from models at the station site (as `ModelTimeSeries`).

All stations are assume to have a unique name. The data of a sensor is stored in `SensorData` objects, and all the metadata is stored by the `Site` object of the station.

You cannot create a new station, but you can extract a station from a `Dataset`.
A :ref:`Station <Station api>` is a class that has the same attributes and methods as a Dataset, but all the observations are limited to a specific station.

.. code-block:: python

   your_station = your_dataset.get_station(stationname = 'station_A')



``SensorData``
--------------
A ``SensorData`` object holds the timeseries data for a specific observation type (e.g., temperature, humidity, wind speed) at a station.
Each station can have multiple SensorData objects, one for each observation type.
SensorData manages the actual measurements, associated timestamps, and quality control labels for its observation type. If present, gaps are stored in the `SensorData`
SensorData objects are not created directly by users; they are managed by the toolkit when importing or processing data.

In pracktiche you do not need to interact directly with this class. You can inspect the observations by using the `df` attribute on a `Station` or `Dataset`.

.. code-block:: python

   records_dataframe = your_dataset.df

See the API documentation :ref:`SensorData API <SensorData api>` for more details.



``Gaps``
--------
A *gap* is a period in the timeseries where observations are missing. Gaps are automatically detected when importing raw data or after resampling/synchronizing timeseries. Each gap is represented by a `Gap` object, which is stored by the corresponding `SensorData`.

Regular users do not interact directly with `Gap` objects. Instead, gaps can be inspected and filled via methods on the `Dataset` and `Station` classes. For example, you can view gaps using the `gapsdf` attribute or visualize them with plotting methods.

.. code-block:: python

   # Inspect gaps for a dataset
   print(your_dataset.gapsdf)


See the API documentation :ref:`Gap API <Gap api>` for more details.


``Analysis``
-----------
The :ref:`Analysis <Analysis api>` class is created from a Dataset and holds the observations that are assumed to be correct. In contrast to the Dataset, the Analysis methods do not change the observations but the focus is on filtering and aggregation.
The Analysis methods are focussed on  aggregating the observations to get insight into diurnal/seasonal patterns and landcover effects.


See the `Analysis example` for more details.

.. code-block:: python

   your_dataset_analysis = metobs_toolkit.Analysis(Dataholder=dataset)

.. note::

   Creating an Analysis of a Station is not recommended, since there is not much scientific value in it.


``Geedatasetmanagers``
----------------------
A ``Geedatasetmanager`` is a class that manages the interaction between the toolkit and a specific dataset on Google Earth Engine (GEE).
These managers do not store modeldata themselves (that is done in the `ModelTimeSeries`), but provide the interface to extract and interpret data from GEE.

There are two types of Geedatasetmanagers:

* ``GEEStaticDatasetManager``: Handles GEE datasets without a time dimension (static). Used to extract static properties (e.g., land cover, altitude, LCZ) at station locations or within buffers.
* ``GEEDynamicDatasetManager``: Handles GEE datasets with a time dimension (dynamic). Used to extract timeseries data (e.g., ERA5 temperature) at station locations. This manager uses ``ModelObstype`` definitions to map GEE dataset bands to observation types and handle unit conversions.

Default managers for common datasets are provided and accessible via the `metobs_toolkit.default_GEE_datasets`. You can also define your own for custom GEE datasets.

See the API documentation :ref:`Geedatasetmanagers API <Geedatasetmanagers api>` and the `Gee example` for more details.


``ModelTimeSeries``
-------------------
A ``ModelTimeSeries`` object stores timeseries data extracted from a dynamic GEE dataset (e.g., ERA5) for a specific observation type at a station. It is similar to the `SensorData` class.
These timeseries represent modelled or reanalysis data, and are typically used for comparison with observations, quality control, or gap filling.

ModelTimeSeries are stored in the `Station` objects.
Regular users do not interact directly with ``ModelTimeSeries`` objects. Instead, modeldata can be inspected via the `.modeldatadf` attribute on the `Dataset` and `Station` classes.

.. code-block:: python

   # Access modelled temperature timeseries for a station
   temp_modeldata = your_station.modeldata['temp']

   # View the timeseries DataFrame
   print(temp_modeldata.df)

   # Plot the modelled data
   temp_modeldata.plot()

See the API documentation :ref:`ModelTimeSeries API <ModelTimeSeries api>` and the `Gee example` for more details.


``Obstype and ModelObstype``
---------------------------
An ``Obstype`` defines an observation type, such as temperature, humidity, or wind speed.
It specifies the standard name, standard unit, and a description for the observation type.
Obstypes are used throughout the toolkit to ensure consistency in data handling, unit conversion, and quality control.

A ``ModelObstype`` extends the concept of an Obstype to model or reanalysis data (e.g., from GEE datasets). In addition to the standard attributes, a ModelObstype defines the corresponding band name and unit in the model dataset. This allows the toolkit to map model data bands to observation types and handle unit conversions automatically.

You typically do not need to create these objects directly; common obstypes and modelobstypes are predefined and used internally by the toolkit and GEE dataset managers.

See the API documentation :ref:`Obstype API <Obstype api>` and :ref:`ModelObstype API <ModelObstype api>` for more details.
