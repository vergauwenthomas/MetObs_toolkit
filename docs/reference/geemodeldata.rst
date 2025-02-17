.. _Geemodeldata api:

===============
GEE Modeldata
===============

The Gee (Google Earth Engine) Modeldata is the bridge between the Metobs-toolkit
and the GEE services. There are two classes:

 * `GeeStaticDataset`: This class handles GEE Datasets that do not have a time dimension (static). This class is used to extract GEE dataset values at the location of the station (or buffers arround them).
 * `GeeDynamicDataset`: This class handles GEE Dataset that have a time dimension. This class is used to extract timeseries of GEE dataset values at the stations locations.

Both classes can hold metadata (=Coordinates of the stations), and the ´GeeDynamicDataset` class can hold timeseries data.

.. note::
    Extracting data from a GEE dataset, can be done directly from a `GeeStaticDataset`
    or a `GeeDynamicDataset`. In addition, one can call extractions also directly from
    a `Dataset`, see the *Extracting data* section in the `Dataset` API documentation.


.. note::
   Keep in mind that the Modeldata does not hold gridded data, but timeseries of
   point extractions at the locations of the stations. Therefore a Modeldata
   instance is linked to a collection of stations defined by the `metadf` attribute.

Demo examples on the Modeldata class can be found here: :ref:`Extracting ERA5 timeseries`_ .


GeeStaticDataset
--------------------

.. currentmodule:: metobs_toolkit

.. autosummary::
   :toctree: api/

   GeeStaticDataset
   GeeStaticDataset.get_info
   GeeStaticDataset.extract_static_point_data
   GeeStaticDataset.extract_static_buffer_frac_data
   GeeStaticDataset.make_gee_plot


GeeDynamicDataset
--------------------

.. autosummary::
   :toctree: api/

   GeeDynamicDataset
   GeeDynamicDataset.get_info
   GeeDynamicDataset.add_modelobstype
   GeeDynamicDataset.extract_timeseries_data
   GeeDynamicDataset.save_modeldata
   GeeDynamicDataset.set_modeldata_from_csv
   import_modeldata_from_pkl
   GeeDynamicDataset.make_plot
   GeeDynamicDataset.make_gee_plot
