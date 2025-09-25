.. _Site api:

======
Site
======

The Site holds all spatial metadata related to a station. This includes:

* geographic coordinates (latitude, longitude, altitude),
* Local Climate Zone (LCZ) information,
* Google Earth Engine (GEE) extracted data,
* and additional metadata.

The Site is a component of the ``Station`` class and provides spatial context for meteorological observations.

.. currentmodule:: metobs_toolkit.site

Constructor
-----------

.. autosummary::
   :toctree: api/

   Site

Data attributes
----------------
A summary of all the attributes that hold or return spatial and metadata information.

.. autosummary::
   :toctree: api/

   Site.stationname
   Site.lat
   Site.lon
   Site.altitude
   Site.LCZ
   Site.buffered_fractions
   Site.extradata
   Site.metadf

Setters
-------
Methods for updating site metadata and spatial information.

.. autosummary::
   :toctree: api/

   Site.set_altitude
   Site.set_LCZ
   Site.set_geedata
   Site.set_gee_buffered_frac_data
   Site.add_metadata

Status check methods
--------------------
Methods to check the availability of specific metadata.

.. autosummary::
   :toctree: api/

   Site.flag_has_altitude
   Site.flag_has_LCZ


General methods
---------------

.. autosummary::
   :toctree: api/

   Site.get_info


Special methods
---------------

The `Site` class implements several Python special methods for convenience:

- ``__add__``: Merge two Site objects, combining metadata and preferring non-NaN values from the first site.

- ``__eq__``: Test equality between two Site objects based on location and metadata.

- ``__repr__``: String representation showing site name and key spatial information.



Notes
-----

The Site class integrates with Google Earth Engine to extract environmental data such as:

- Static point data (elevation, land cover)
- Buffered fraction data (land cover fractions within specified radii)
- Local Climate Zone classifications

All spatial coordinates are stored in WGS84 coordinate system, and the ``metadf`` property returns a GeoDataFrame for spatial operations.