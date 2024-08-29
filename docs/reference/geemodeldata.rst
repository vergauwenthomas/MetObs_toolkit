.. _Geemodeldata api:

===============
GEE Modeldata
===============

The Gee (Google Earth Engine) Modeldata is the bridge between the Metobs-toolkit
and the GEE services. There are two classes:

 * `GeeStaticModelData`: This class handles GEE Datasets that do not have a time dimension (static). This class is used to extract GEE dataset values at the location of the station (or buffers arround them).
 * `GeeDynamicModelData`: This class handles GEE Dataset that have a time dimension. This class is used to extract timeseries of GEE dataset values at the stations locations.

Both classes can hold metadata (=Coordinates of the stations), and the ´GeeDynamicModelData` class can hold timeseries data.

.. note::
    Extracting data from a GEE dataset, can be done directly from a `GeeStaticModelData`
    or a `GeeDynamicModelData`. In addition, one can call extractions also directly from
    a `Dataset`, see the *Extracting data* section in the `Dataset` API documentation.


.. note::
   Keep in mind that the Modeldata does not hold gridded data, but timeseries of
   point extractions at the locations of the stations. Therefore a Modeldata
   instance is linked to a collection of stations defined by the `metadf` attribute.

Demo examples on the Modeldata class can be found here: :ref:`Extracting ERA5 timeseries`_ .


GeeStaticModelData
--------------------

.. currentmodule:: metobs_toolkit

.. autosummary::
   :toctree: api/

   GeeStaticModelData
   GeeStaticModelData.get_info
   GeeStaticModelData.extract_static_point_data
   GeeStaticModelData.extract_static_buffer_frac_data
   GeeStaticModelData.make_gee_plot


GeeDynamicModelData
--------------------

.. autosummary::
   :toctree: api/

   GeeDynamicModelData
   GeeDynamicModelData.get_info
   GeeDynamicModelData.add_modelobstype
   GeeDynamicModelData.extract_timeseries_data
   GeeDynamicModelData.save_modeldata
   GeeDynamicModelData.set_modeldata_from_csv
   GeeDynamicModelData.import_modeldata_from_pkl
   GeeDynamicModelData.make_plot
   GeeDynamicModelData.make_gee_plot
