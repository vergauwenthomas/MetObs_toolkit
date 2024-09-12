#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 17:14:10 2024

@author: thoverga
"""


import logging
import pandas as pd
import pytz

# from metobs_toolkit.landcover_functions import (
#     connect_to_gee,
#     lcz_extractor,
#     height_extractor,
#     lc_fractions_extractor,
# )

from metobs_toolkit.modeldata import (
    GeeStaticModelData,
    GeeDynamicModelData,
)

from metobs_toolkit.gee_api import connect_to_gee


logger = logging.getLogger(__name__)


class DatasetModelData:
    """Extension on the metobs_toolkit.Dataset class with GEE modeldata related methods"""

    def get_modeldata(
        self,
        Model,
        obstypes=["temp"],
        stations=None,
        startdt=None,
        enddt=None,
        get_all_bands=False,
        drive_filename=None,
        drive_folder="gee_timeseries_data",
        force_direct_transfer=False,
        force_to_drive=False,
    ):
        """Extract Timeseries data from a Gee dataset at your stations.

        The link with a Gee dataset is done by specifying a GeeDyanmicModelData.
        The extracted data is returned in the form of a DataFrame and is stored
        in the Model.


        Parameters
        ----------
        Model : str or GeeDynamicModelData
            The Gee dataset to download timeseries from. If Model is a str, it
            is looked for by name in the `Dataset.gee_dataset.keys()`.
        obstypes : str or list of strings, optional
            The name of the ModelObstypes (or ModelObstype_Vectorfield) to extract timeseries for. The requested
            ModelObstypes must be known by the Model. The default is ['temp'].
        stations : string or list of strings, optional
            Stationnames to subset the modeldata to. If None, all stations will be used. The default is None.
        startdt : datetime.datetime, optional
            Start datetime of the model timeseries. If None, the start datetime of the dataset is used.
            A startdt must be specified when the dataset does not have data records (metadata-only). The default is None.
        enddt : datetime.datetime, optional
            End datetime of the model timeseries. If None, the last datetime of the dataset is used.
            A enddt must be specified when the dataset does not have data records (metadata-only). The default is None.
        get_all_bands : bool, optional
            If True, all values (over all bands) are extracted. If the band is
            linked to a ModelObstye, then the name of the modelObstype is used
            instead of the bandname. If True, the obstypes argument is ignored.
            The default is False.
        drive_filename : str or None, optional
            If given, the data will be saved as this filename on your Google Drive.
            This argument will only take effect when the data is written to
            Google Drive. If None, a custom filename is created. The default is
            None.
        drive_folder: str
            The name of the folder on your Google Drive to save the drive_filename
            in. If the folder, does not exist it will be created instead. This
            argument will only take effect when the data is written to Google
            Drive.
        force_direct_transfer: bool, optional
            If True, the data is demanded as a direct transfer (no file written
            to Google Drive). If the request is tol large, a GEE error is raised.
            The default is False.
        force_to_drive: bool, optional
            If True, The gee data is written to a file on your drive. Direct
            transfer of data is prohibited. The default is False.


        Returns
        -------
        Model : metobs_toolkit.GeeDynamicModelData
            The Gee Modeldata with the extracted timeseries stored in. The
            timeseries are stored in the `Model.modeldf`.

        Note
        --------
        If a timezone unaware datetime is given as an argument, it is interpreted
        as if it has the same timezone as the observations.

        Note
        ------
        When extracting large amounts of data, the timeseries data will be
        written to a file and saved on your Google Drive. In this case, you need
        to provide the Modeldata with the data using the .set_model_from_csv()
        method.


        See Also
        -----------
        GeeDynamicModelData: Gee Modeldata dataset for time-evolving data.
        metobs_toolkit.ModelObstype: A Obstype for Modeldata (linked to a band).
        metobs_toolkit.ModelObstype_Vectorfield: A Vectorfield version of ModelObstype.
        Dataset.add_new_geemodeldata(): Add a new Gee Modeldata to your dataset.

        Examples
        --------

        Create a ``Dataset`` and fill it with data (and metadata).

        >>> import metobs_toolkit
        >>>
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> dataset.import_data_from_file(
        ...                         input_data_file=metobs_toolkit.demo_datafile,
        ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                         template_file=metobs_toolkit.demo_template,
        ...                         )

        We will now extract modeldata, directly through the use of the GEE (
        Google Earht Engine) API. The Modeldata will extract the timeseries,
        of the stations present in the Dataset.

        By default, each `Dataset` is equipped with default Gee datasets.

        >>> dataset.gee_datasets
        {'lcz': GeeStaticModelData instance of lcz  (no metadata has been set) , 'altitude': GeeStaticModelData instance of altitude  (no metadata has been set) , 'worldcover': GeeStaticModelData instance of worldcover  (no metadata has been set) , 'ERA5-land': Empty GeeDynamicModelData instance of ERA5-land }

        As you can see, is "ERA5-land" a default GeeDynamicModelData. We will
        use it for this example.


        If the data transfer is to big, a file .csv file is written in your
        Google Drive. You must download that file, and import it using the
        ``Modeldata.set_model_from_csv()`. To limit the transfer of data,
        we will download timeseries for a single station, and a specific timeperiod.

        >>> import datetime
        >>>
        >>> tstart = datetime.datetime(2022, 9, 5)
        >>> tend = datetime.datetime(2022, 9, 6)
        >>>
        >>> sta = dataset.get_station('vlinder02')

        We specify the Gee dataset to use

        >>> ERA5_model = sta.gee_datasets['ERA5-land']
        >>> ERA5_model.get_info()
        Empty GeeDynamicModelData instance of ERA5-land
        ------ Details ---------
        <BLANKLINE>
         * name: ERA5-land
         * location: ECMWF/ERA5_LAND/HOURLY
         * value_type: numeric
         * scale: 2500
         * is_static: False
         * is_image: False
         * is_mosaic: False
         * credentials:
         * time res: 1h
        <BLANKLINE>
         -- Known Modelobstypes --
        <BLANKLINE>
         * temp : ModelObstype instance of temp (linked to band: temperature_2m)
            (conversion: Kelvin --> Celsius)
         * pressure : ModelObstype instance of pressure (linked to band: surface_pressure)
            (conversion: pa --> pa)
         * wind : ModelObstype_Vectorfield instance of wind (linked to bands: u_component_of_wind_10m and v_component_of_wind_10m)
            vectorfield that will be converted to:
              * wind_speed
              * wind_direction
            (conversion: m/s --> m/s)
        <BLANKLINE>
         -- Metadata --
        <BLANKLINE>
        No metadf is set.
        <BLANKLINE>
         -- Modeldata --
        <BLANKLINE>
        No model data is set.

        Now we download temperature timeseries of ERA5 data at the location
        of "vlinder02" for the period of interest.

        >>> # Collect ERA5 2mT timeseries at your station
        >>> sta.get_modeldata(
        ...                     Model= ERA5_model,
        ...                     obstypes=["temp"],
        ...                     startdt=tstart,
        ...                     enddt=tend,
        ...                     force_direct_transfer=True)
        GeeDynamicModelData instance of ERA5-land with modeldata

        The timeseries are stored in the Model itself.

        >>> ERA5_model.modeldf.head()
                                              temp
        name      datetime
        vlinder02 2022-09-05 00:00:00+00:00  20.27
                  2022-09-05 01:00:00+00:00  20.11
                  2022-09-05 02:00:00+00:00  19.74
                  2022-09-05 03:00:00+00:00  19.54
                  2022-09-05 04:00:00+00:00  19.00

        ERA5_data contains 1 timeseries of temperature data, automatically
        converted to the toolkit standard unit (Celcius).

        """
        # Model conversion and check
        if isinstance(Model, str):
            if Model not in self.gee_datasets.keys():
                raise MetobsDatasetGeeModelDataHandlingError(
                    f"{Model} is not a known GeeDynamicModelData of {self}."
                )
            Model = self.gee_datasets[str(Model)]
        elif isinstance(Model, GeeDynamicModelData):
            pass
        else:
            raise MetobsDatasetGeeModelDataHandlingError(
                f"{Model} is not a (known) GeeDynamicModelData (of {self})."
            )

        # Check obstypes
        if isinstance(obstypes, str):
            obstypes = [obstypes]  # convert to list
        for obstype in obstypes:
            if obstype not in Model.modelobstypes.keys():
                raise MetobsDatasetGeeModelDataHandlingError(
                    f"{obstype} is not a known Modelobstype of {Model}."
                )

        # Filters
        if startdt is None:
            if self.df.empty:
                raise MetobsDatasetGeeModelDataHandlingError(
                    'Specifying a startdt is required for a "metadata-only" dataset.'
                )
            startdt = self.df.index.get_level_values("datetime").min()
        else:
            startdt = self._datetime_arg_check(startdt)

        if enddt is None:
            if self.df.empty:
                raise MetobsDatasetGeeModelDataHandlingError(
                    'Specifying a enddt is required for a "metadata-only" dataset.'
                )
            enddt = self.df.index.get_level_values("datetime").max()
        else:
            enddt = self._datetime_arg_check(enddt)

        # # make sure bounds include required range
        startdt = startdt.floor(Model.time_res)
        enddt = enddt.ceil(Model.time_res)

        if stations is not None:
            if isinstance(stations, str):
                metadf = self.metadf.loc[[stations]]
            if isinstance(stations, list):
                metadf = self.metadf.iloc[self.metadf.index.isin(stations)]
        else:
            metadf = self.metadf

        Model.set_metadf(metadf)

        Model.extract_timeseries_data(
            startdt_utc=startdt,
            enddt_utc=enddt,
            obstypes=obstypes,
            get_all_bands=get_all_bands,
            drive_filename=drive_filename,
            drive_folder=drive_folder,
            force_direct_transfer=force_direct_transfer,
            force_to_drive=force_to_drive,
        )

        return Model

    def get_lcz(self):
        """Extract local climate zones for all stations.

        Function to extract the Local CLimate zones (LCZ) from the
        wudapt global LCZ map on the Google engine for all stations.

        A 'LCZ' column will be added to the metadf, and the series is returned.

        Returns
        -------
        lcz_series : pandas.Series()
            A series with the stationnames as index and the LCZ as values.

        Warning
        ---------
        This method makes use of GEE API. Make sure that you have accesss and
        user rights to use the GEE API.

        See Also
        --------
        connect_to_gee: Set up a new connection/credentials to the GEE service.
        get_altitude: Extract altitudes for all stations
        get_landcover: Extract landcover fractions for all stations


        Examples
        --------
        Create a ``Dataset`` and fill it with data (and metadata).

        >>> import metobs_toolkit
        >>>
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> dataset.import_data_from_file(
        ...                         input_data_file=metobs_toolkit.demo_datafile,
        ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                         template_file=metobs_toolkit.demo_template,
        ...                         )

        At this point, there is no LCZ information present in the metadata

        >>> dataset.metadf.head()
                     lat   lon        school                  geometry dataset_resolution                  dt_start                    dt_end
        name
        vlinder01  50.98  3.82         UGent  POINT (3.81576 50.98044)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder02  51.02  3.71         UGent   POINT (3.7097 51.02238)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder03  51.32  4.95   Heilig Graf  POINT (4.95211 51.32458)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder04  51.34  4.93   Heilig Graf  POINT (4.93473 51.33552)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder05  51.05  3.68  Sint-Barbara  POINT (3.67518 51.05266)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00

        We use the GEE API to extract the LCZ for all the stations present in
        the metadf (if coordinates are present).

        >>> lcz_series = dataset.get_lcz()
        >>> lcz_series.head()
                                  lcz
        name
        vlinder01  Low plants (LCZ D)
        vlinder02       Large lowrise
        vlinder03        Open midrise
        vlinder04      Sparsely built
        vlinder05       Water (LCZ G)

        The LCZ are automatically added to the metadf as 'lcz' column.

        >>> dataset.metadf['lcz']
         name
        vlinder01     Low plants (LCZ D)
        vlinder02          Large lowrise
        vlinder03           Open midrise
        vlinder04         Sparsely built
        vlinder05          Water (LCZ G)
                            ...
        vlinder24    Dense Trees (LCZ A)
        vlinder25          Water (LCZ G)
        vlinder26           Open midrise
        vlinder27        Compact midrise
        vlinder28           Open lowrise
        Name: lcz, Length: 28, dtype: object

        """

        # Get lcz modeldata class
        if "lcz" not in self.gee_datasets.keys():
            raise MetobsDatasetGeeModelDataHandlingError(
                f"No lcz modeldataset is defined in the .gee_datasets: {self.gee_datasets}"
            )
        modl = self.gee_datasets["lcz"]

        # Set the metadf of the model
        modl.set_metadf(self.metadf)

        # connect to gee
        connect_to_gee()

        # Extract lcz
        lczdf = modl.extract_static_point_data()

        # drop column if it was already present
        if "lcz" in self.metadf:
            self.metadf = self.metadf.drop(columns=["lcz"])

        # update metadata
        self.metadf = self.metadf.merge(
            lczdf,
            how="left",
            left_index=True,
            right_index=True,
        )
        return lczdf

    def get_altitude(self):
        """Extract Altitudes for all stations.

        Function to extract the Altitude from the SRTM Digital Elevation Data
        global map on the Google engine for all stations.

        An 'altitude' column will be added to the metadf, and the series is returned.

        Returns
        -------
        altitude_series : pandas.Series()
            A series with the stationnames as index and the altitudes as values.

        Warning
        ---------
        This method makes use of GEE API. Make sure that you have access and
        user rights to use the GEE API.

        See Also
        --------
        connect_to_gee: Set up a new connection/credentials to the GEE service.
        get_lcz: Extract LCZ for all stations
        get_landcover: Extract landcover fractions for all stations

        Examples
        --------
        Create a ``Dataset`` and fill it with data (and metadata).

        >>> import metobs_toolkit
        >>>
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> dataset.import_data_from_file(
        ...                         input_data_file=metobs_toolkit.demo_datafile,
        ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                         template_file=metobs_toolkit.demo_template,
        ...                         )

        At this point, there is no altitude information present in the metadata

        >>> dataset.metadf.head()
                     lat   lon        school                  geometry dataset_resolution                  dt_start                    dt_end
        name
        vlinder01  50.98  3.82         UGent  POINT (3.81576 50.98044)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder02  51.02  3.71         UGent   POINT (3.7097 51.02238)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder03  51.32  4.95   Heilig Graf  POINT (4.95211 51.32458)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder04  51.34  4.93   Heilig Graf  POINT (4.93473 51.33552)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder05  51.05  3.68  Sint-Barbara  POINT (3.67518 51.05266)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00

        We use the GEE API to extract the LCZ for all the stations present in
        the metadf (if coordinates are present).

        >>> alt_series = dataset.get_altitude()
        >>> alt_series.head()
                   altitude
        name
        vlinder01        12
        vlinder02         7
        vlinder03        30
        vlinder04        25
        vlinder05         0

        The altitudes are automatically added to the metadf as 'altitude' column.

        >>> dataset.metadf['altitude']
        name
        vlinder01    12
        vlinder02     7
        vlinder03    30
        vlinder04    25
        vlinder05     0
                     ..
        vlinder24    12
        vlinder25    12
        vlinder26    24
        vlinder27    12
        vlinder28     7
        Name: altitude, Length: 28, dtype: int64

        """
        # Get lcz modeldata class
        if "altitude" not in self.gee_datasets.keys():
            raise MetobsDatasetGeeModelDataHandlingError(
                f"No altitude (DEM) modeldataset is defined in the .gee_datasets: {self.gee_datasets}"
            )
        modl = self.gee_datasets["altitude"]

        # Set the metadf of the model
        modl.set_metadf(self.metadf)

        # connect to gee
        connect_to_gee()

        # Extract altitude
        demdf = modl.extract_static_point_data()

        # drop column if it was already present
        if "altitude" in self.metadf:
            self.metadf = self.metadf.drop(columns=["altitude"])

        # update metadata
        self.metadf = self.metadf.merge(
            demdf,
            how="left",
            left_index=True,
            right_index=True,
        )
        return demdf

    def get_landcover(
        self,
        buffers=[100],
        aggregate=False,
        overwrite=True,
        gee_map="worldcover",
    ):
        """Extract landcover for all stations.

        Extract the landcover fractions in a buffer with a specific radius for
        all stations. If an aggregation scheme is defined, one can choose to
        aggregate the landcover classes.

        The landcover fractions will be added to the Dataset.metadf if overwrite
        is True. Presented as separate columns where each column represents the
        landcovertype and corresponding buffer.

        Parameters
        ----------
        buffers : list of numerics, optional
            The list of buffer radia in dataset units (meters for ESA worldcover) . The default is 100.
        aggregate : bool, optional
            If True, the classes will be aggregated with the corresponding
            aggregation scheme. The default is False.
        overwrite : bool, optional
            If True, the Datset.metadf will be updated with the generated
            landcoverfractions. The default is True.
        gee_map : str, optional
            The name of the dataset to use. This name should be present in the
            settings.gee['gee_dataset_info']. If aggregat is True, an aggregation
            scheme should included as well. The default is 'worldcover'

        Returns
        -------
        frac_df : pandas.DataFrame
            A Dataframe with index: name, buffer_radius and the columns are the
            fractions.

        Warning
        ---------
        This method makes use of GEE API. Make sure that you have access and
        user rights to use the GEE API.

        Warning
        ---------
        It can happen that for stations located on small islands, or close to
        the coast, the sea-mask is not used as a landcover fraction.

        See Also
        --------
        connect_to_gee: Set up a new connection/credentials to the GEE service.
        get_altitude: Extract altitudes for all stations
        get_lcz: Extract lcz for all stations
        make_gee_plot: Make an interactive plot of a GEE dataset

        Examples
        --------
        Create a ``Dataset`` and fill it with data (and metadata).

        >>> import metobs_toolkit
        >>>
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> dataset.import_data_from_file(
        ...                         input_data_file=metobs_toolkit.demo_datafile,
        ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                         template_file=metobs_toolkit.demo_template,
        ...                         )

        At this point, there is no landcover information present in the metadata.

        >>> dataset.metadf.head()
                     lat   lon        school                  geometry dataset_resolution                  dt_start                    dt_end
        name
        vlinder01  50.98  3.82         UGent  POINT (3.81576 50.98044)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder02  51.02  3.71         UGent   POINT (3.7097 51.02238)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder03  51.32  4.95   Heilig Graf  POINT (4.95211 51.32458)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder04  51.34  4.93   Heilig Graf  POINT (4.93473 51.33552)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder05  51.05  3.68  Sint-Barbara  POINT (3.67518 51.05266)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00

        We use the GEE API to extract the landcover fractions for all the
        stations present in the metadf (if coordinates are present).

        >>> frac_df = dataset.get_landcover(buffers=[50, 100, 250, 500],
        ...                                 aggregate=False)
        >>> frac_df
                                 Grassland  Cropland  Tree cover  Built-up  Permanent water bodies  Herbaceous wetland  Bare / sparse vegetation  Shrubland
        name      buffer_radius
        vlinder01 50              5.46e-01  4.54e-01    0.00e+00      0.00                0.00e+00                 0.0                  0.00e+00        0.0
                  100             3.46e-01  6.36e-01    0.00e+00      0.02                0.00e+00                 0.0                  0.00e+00        0.0
                  250             3.19e-01  6.41e-01    4.21e-03      0.04                0.00e+00                 0.0                  0.00e+00        0.0
                  500             3.90e-01  4.66e-01    4.91e-02      0.09                2.43e-04                 0.0                  0.00e+00        0.0
        vlinder02 50              6.29e-01  0.00e+00    1.43e-02      0.36                0.00e+00                 0.0                  0.00e+00        0.0
        ...                            ...       ...         ...       ...                     ...                 ...                       ...        ...
        vlinder27 500             4.31e-03  0.00e+00    6.11e-02      0.92                1.08e-02                 0.0                  3.23e-04        0.0
        vlinder28 50              4.99e-02  0.00e+00    1.60e-01      0.79                0.00e+00                 0.0                  0.00e+00        0.0
                  100             1.88e-01  0.00e+00    3.02e-01      0.51                0.00e+00                 0.0                  0.00e+00        0.0
                  250             1.28e-01  0.00e+00    5.94e-01      0.28                0.00e+00                 0.0                  0.00e+00        0.0
                  500             1.16e-01  1.62e-04    5.49e-01      0.34                0.00e+00                 0.0                  0.00e+00        0.0
        <BLANKLINE>
        [112 rows x 8 columns]


        The landcover fractions are automatically added to the metadf.


        >>> dataset.metadf.columns
        Index(['lat', 'lon', 'school', 'geometry', 'dataset_resolution', 'dt_start',
           'dt_end', 'Grassland_100m', 'Grassland_250m', 'Grassland_50m',
           'Grassland_500m', 'Cropland_100m', 'Cropland_250m', 'Cropland_50m',
           'Cropland_500m', 'Tree cover_100m', 'Tree cover_250m', 'Tree cover_50m',
           'Tree cover_500m', 'Built-up_100m', 'Built-up_250m', 'Built-up_50m',
           'Built-up_500m', 'Permanent water bodies_100m',
           'Permanent water bodies_250m', 'Permanent water bodies_50m',
           'Permanent water bodies_500m', 'Herbaceous wetland_100m',
           'Herbaceous wetland_250m', 'Herbaceous wetland_50m',
           'Herbaceous wetland_500m', 'Bare / sparse vegetation_100m',
           'Bare / sparse vegetation_250m', 'Bare / sparse vegetation_50m',
           'Bare / sparse vegetation_500m', 'Shrubland_100m', 'Shrubland_250m',
           'Shrubland_50m', 'Shrubland_500m'],
          dtype='object')


        """

        # Get lcz modeldata class
        if gee_map not in self.gee_datasets.keys():
            raise MetobsDatasetGeeModelDataHandlingError(
                f"No {gee_map} modeldataset is defined in the .gee_datasets: {self.gee_datasets}"
            )
        modl = self.gee_datasets[gee_map]

        # Set the metadf of the model
        modl.set_metadf(self.metadf)

        # connect to gee
        connect_to_gee()

        # Extract altitude
        # demdf = modl.extract_static_data()
        # return modl

        # connect to gee
        connect_to_gee()

        df_list = []
        for buffer in buffers:
            logger.info(
                f"Extracting landcover from {gee_map} with buffer radius = {buffer}"
            )
            # Extract landcover fractions for all stations
            freqsdf = modl.extract_static_buffer_frac_data(
                bufferradius=buffer, agg_bool=aggregate
            )
            df_list.append(freqsdf)

        lcdf = pd.concat(df_list).fillna(0)
        lcdf = lcdf.sort_index()

        # create a dataframe with buffer reference in the column name and name as index
        flatted_lcdf = lcdf.reset_index()
        flatted_lcdf["buffer_radius"] = (
            flatted_lcdf["buffer_radius"].astype(int).astype(str)
        )
        flatted_lcdf = flatted_lcdf.set_index(["name", "buffer_radius"])
        flatted_lcdf = flatted_lcdf.unstack()
        flatted_lcdf.columns = flatted_lcdf.columns.map("_".join)
        flatted_lcdf.columns = [col + "m" for col in flatted_lcdf.columns]

        # Update self.metadf
        self.metadf.drop(columns=flatted_lcdf.columns, errors="ignore", inplace=True)
        self.metadf = self.metadf.merge(
            flatted_lcdf, how="left", left_index=True, right_index=True
        )

        return lcdf

    def add_new_geemodeldata(self, Modeldata, overwrite=False):
        """Add a new GeeModeldata to the Dataset.

        The GeeModelData is in the form of a GeeStaticModelData or GeeDynamicModelData.


        Parameters
        ----------
        Modeldata : GeeStaticModelData or GeeDynamicModelData
            The GeeModeldata to add.
        overwrite : bool, optional
            If there is already a GeeModelData with the same name present
            in the Dataset, and if overwrite is True, then it will be
            overwritten. The default is False.

        Returns
        -------
        None.

        See Also
        -----------
        GeeStaticModelData: Gee Modeldata dataset without time dimension.
        GeeDynamicModelData: Gee Modeldata dataset for time-evolving data.

        Examples
        ----------
        As an example, we will add a new Gee Modeldata (precipitation satellite
        product), to the Dataset.

        >>> import metobs_toolkit
        >>>
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> dataset.gee_datasets
        {'lcz': GeeStaticModelData instance of lcz  (no metadata has been set) , 'altitude': GeeStaticModelData instance of altitude  (no metadata has been set) , 'worldcover': GeeStaticModelData instance of worldcover  (no metadata has been set) , 'ERA5-land': Empty GeeDynamicModelData instance of ERA5-land }

        As we can see, there are default GeeModelData's present in the Dataset.

        Now we will create a new GeeDynamicModelData linking to this GEE dataset:
        https://developers.google.com/earth-engine/datasets/catalog/JAXA_GPM_L3_GSMaP_v8_operational.

        We create a ModelObstype reflecting the precipitation, and
        linking it with the corresponding band of the GEE dataset.

        >>> dataset.obstypes #See which obstypes are already present
        {'temp': Obstype instance of temp, 'humidity': Obstype instance of humidity, 'radiation_temp': Obstype instance of radiation_temp, 'pressure': Obstype instance of pressure, 'pressure_at_sea_level': Obstype instance of pressure_at_sea_level, 'precip': Obstype instance of precip, 'precip_sum': Obstype instance of precip_sum, 'wind_speed': Obstype instance of wind_speed, 'wind_gust': Obstype instance of wind_gust, 'wind_direction': Obstype instance of wind_direction}

        >>> hourly_precip = dataset.obstypes['precip']

        Now we have a regular obstype for precipitation. To use in a Gee Modeldata,
        we must make a ModelObstype from it (=add the link with the band).

        >>> hourly_precip_Geedataset = metobs_toolkit.ModelObstype(
        ...               obstype=hourly_precip,
        ...               model_unit='mm', #See GEE
        ...               model_band="hourlyPrecipRateGC")

        (Add a new unit to the regular Obstype if the model_unit is not known.)

        Now we can create the Gee Modeldata.

        >>> precip_satelite = metobs_toolkit.GeeDynamicModelData(
        ...        name='precip_GSMaP',
        ...        location="JAXA/GPM_L3/GSMaP/v8/operational", #See GEE
        ...        value_type='numeric', #See GEE
        ...        scale=50, #See GEE
        ...        time_res='1h', #See GEE
        ...        modelobstypes=[hourly_precip_Geedataset],
        ...        is_image=False, #See GEE
        ...        is_mosaic=False)

        Now we can add it to the Dataset.

        >>> dataset.add_new_geemodeldata(Modeldata=precip_satelite)
        >>> dataset.gee_datasets
        {'lcz': GeeStaticModelData instance of lcz  (no metadata has been set) , 'altitude': GeeStaticModelData instance of altitude  (no metadata has been set) , 'worldcover': GeeStaticModelData instance of worldcover  (no metadata has been set) , 'ERA5-land': Empty GeeDynamicModelData instance of ERA5-land , 'precip_GSMaP': Empty GeeDynamicModelData instance of precip_GSMaP }

        Now we can use the Gee Modeldata with our dataset (gap filling, plotting,
        extracting timeseries, ...). As an example we will download the timeseries
        of precipitation (form the satelite product) at the locations of the stations.

        First, we need to import the metadata.

        >>> dataset.import_data_from_file(
        ...         input_data_file=metobs_toolkit.demo_datafile,
        ...         input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...         template_file=metobs_toolkit.demo_template)
        >>> dataset.metadf.head()
                     lat   lon        school                  geometry dataset_resolution                  dt_start                    dt_end
        name
        vlinder01  50.98  3.82         UGent  POINT (3.81576 50.98044)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder02  51.02  3.71         UGent   POINT (3.7097 51.02238)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder03  51.32  4.95   Heilig Graf  POINT (4.95211 51.32458)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder04  51.34  4.93   Heilig Graf  POINT (4.93473 51.33552)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder05  51.05  3.68  Sint-Barbara  POINT (3.67518 51.05266)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00

        Now we can extract timeseries at the location of the stations.

        >>> import datetime
        >>> tstart = datetime.datetime(2022,9,4)
        >>> tend = datetime.datetime(2022,9,5)
        >>> dataset.get_modeldata(
        ...                  Model='precip_GSMaP',
        ...                  obstypes=['precip'],
        ...                  startdt=tstart,
        ...                  enddt=tend)
        GeeDynamicModelData instance of precip_GSMaP with modeldata

        >>> dataset.gee_datasets['precip_GSMaP'].modeldf.head()
                                                 precip
        name      datetime
        vlinder01 2022-09-04 00:00:00+00:00       0
                  2022-09-04 01:00:00+00:00       0
                  2022-09-04 02:00:00+00:00       0
                  2022-09-04 03:00:00+00:00       0
                  2022-09-04 04:00:00+00:00       0

        """

        # Check instance
        if not (
            (isinstance(Modeldata, GeeStaticModelData))
            | (isinstance(Modeldata, GeeDynamicModelData))
        ):
            raise MetobsDatasetGeeModelDataHandlingError(
                f"{Modeldata} is not an instance of GeeStaticModelData or GeeDynamicModelData."
            )
        # Check name is unique
        if Modeldata.name in self.gee_datasets.keys():
            if overwrite:
                logger.info(f"Overwriting the {Modeldata.name} known Modeldata.")
            else:
                raise MetobsDatasetGeeModelDataHandlingError(
                    f"{Modeldata.name} is already a known name of a Modeldata."
                )

        self.gee_datasets[Modeldata.name] = Modeldata


# =============================================================================
# Errors
# =============================================================================


class MetobsDatasetGeeModelDataHandlingError(Exception):
    """Exception raised for errors in the Dataset - GEE modeldata interactions."""

    pass


# =============================================================================
# Docstring test
# =============================================================================
if __name__ == "__main__":
    from metobs_toolkit.doctest_fmt import setup_and_run_doctest

    setup_and_run_doctest()
