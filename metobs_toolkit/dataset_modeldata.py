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
    connect_to_gee,
)


logger = logging.getLogger(__name__)


class DatasetModelData:
    """Extension on the metobs_toolkit.Dataset class with GEE modeldata related methods"""

    def get_modeldata(
        self,
        Model,
        obstypes="temp",
        stations=None,
        startdt=None,
        enddt=None,
        get_all_bands=False,
        drive_filename=None,
        drive_folder="gee_timeseries_data",
    ):
        # TODO: update docstring
        """Make Modeldata for the Dataset.

        Make a metobs_toolkit.Modeldata object with modeldata at the locations
        of the stations present in the dataset. This Modeldata stores timeseries
        of model data for each station.

        Parameters
        ----------
        modelname : str, optional
            Which dataset to download timeseries from. This is only used when
            no modeldata is provided. The default is 'ERA5_hourly'.
        modeldata : metobs_toolkit.Modeldata, optional
            Use the modelname attribute and the gee information stored in the
            modeldata instance to extract timeseries.
        obstype : String, optional
            Name of the observationtype you want to apply gap filling on. The
            modeldata must contain this observation type as well. The
            default is 'temp'.
        stations : string or list of strings, optional
            Stationnames to subset the modeldata to. If None, all stations will be used. The default is None.
        startdt : datetime.datetime, optional
            Start datetime of the model timeseries. If None, the start datetime of the dataset is used. The default is None.
        enddt : datetime.datetime, optional
            End datetime of the model timeseries. If None, the last datetime of the dataset is used. The default is None.

        Returns
        -------
        Modl : metobs_toolkit.Modeldata
            The extracted modeldata for period and a set of stations.

        Note
        --------
        If a timezone unaware datetime is given as an argument, it is interpreted
        as if it has the same timezone as the observations.

        Note
        ------
        When extracting large amounts of data, the timeseries data will be
        written to a file and saved on your google drive. In this case, you need
        to provide the Modeldata with the data using the .set_model_from_csv()
        method.

        See Also
        --------

        Modeldata: Modeldata class.
        Modeldata.add_obstype: add a new obstype and band to the Modeldata.
        Modeldata.add_gee_dataset: add a new Google earth engine Modeldata dataset.
        Modeldata.set_model_from_csv: Read GEE modeldata from a csv file.
        fill_gaps_with_raw_modeldata: Raw modeldata gapfill method.


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

        We will now extract modeldata, directly trough the use of the GEE (
        google earht engine) API. The Modeldata will extract timeseries,
        of the stations present in the Dataset.

        If the data transfer is to big, a file .csv file is writen in your
        google Drive. You must download that file, and import it using the
        ``Modeldata.set_model_from_csv()`. To limit the transfer of data,
        we will dowload timeseries for a single station, and a specific timeperiod.

        >>> import datetime
        >>>
        >>> tstart = datetime.datetime(2022, 9, 5)
        >>> tend = datetime.datetime(2022, 9, 6)
        >>>
        >>> sta = dataset.get_station('vlinder02')

        Now we download temperature timeseries of ERA5 data at the location
        of "vlinder02" for the period of interest.

        >>> # Collect ERA5 2mT timeseries at your station
        >>> ERA5_data = sta.get_modeldata(
        ...                     modelname="ERA5_hourly",
        ...                     modeldata=None,
        ...                     obstype="temp",
        ...                     stations=None,
        ...                     startdt=tstart,
        ...                     enddt=tend)
        (When using the .set_model_from_csv() method, make sure the modelname of your Modeldata is ERA5_hourly)

        ERA5_data contains 1 timeseries of temperature data, automatically
        converted to the toolkit standard unit (Celcius).

        >>> print(ERA5_data)
        Modeldata instance containing:
            * Modelname: ERA5_hourly
            * 1 timeseries
            * The following obstypes are available: ['temp']
            * Data has these units: ['Celsius']
            * From 2022-09-05 00:00:00+00:00 --> 2022-09-06 00:00:00+00:00 (with tz=UTC)
        <BLANKLINE>
        (Data is stored in the .df attribute)

        """
        # Model conversion and check
        if isinstance(Model, str):
            if Model not in self.gee_datasets.keys():
                raise MetobsDatasetGeeModelDataHandlingError(
                    f"{Model} is not a known GeeDynamicModelData of {self}."
                )
            Model = self.gee_datasets[str(Model)]
        elif isinstance(Model, type(GeeDynamicModelData)):
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
                    f"{obstype} is not a knonw Modelobstype of {Model}."
                )

        # Filters
        if startdt is None:
            startdt = self.df.index.get_level_values("datetime").min()
        else:
            startdt = self._datetime_arg_check(startdt)

        if enddt is None:
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
            drive_folder="gee_timeseries_data",
        )

        return Model

    def get_lcz(self):
        """Extract local climate zones for all stations.

        Function to extract the Local CLimate zones (LCZ) from the
        wudapt global LCZ map on the Google engine for all stations.

        A 'LCZ' column will be added to the metadf, and series is returned.

        Returns
        -------
        lcz_series : pandas.Series()
            A series with the stationnames as index and the LCZ as values.

        Warning
        ---------
        This methods makes use of GEE API. Make sure that you have acces and
        user rights to use the GEE API.

        See Also
        --------
        connect_to_gee: Setup a new connection/credentials to the GEE service.
        get_altitude: Extract altitudes for all stations
        get_landcover: Extract landcoverfractions for all stations


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

        At this point there is no LCZ information present in the metadata

        >>> dataset.metadf.head()
                         lat       lon        school                  geometry dataset_resolution                  dt_start                    dt_end
        name
        vlinder01  50.980438  3.815763         UGent  POINT (3.81576 50.98044)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder02  51.022379  3.709695         UGent   POINT (3.7097 51.02238)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder03  51.324583  4.952109   Heilig Graf  POINT (4.95211 51.32458)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder04  51.335522  4.934732   Heilig Graf  POINT (4.93473 51.33552)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder05  51.052655  3.675183  Sint-Barbara  POINT (3.67518 51.05266)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00

        We use the GEE API to extract the LCZ for all the stations present in
        the metadf (if coordinates are present).

        >>> lcz_series = dataset.get_lcz()
        >>> lcz_series.head()
            name
        vlinder01    Low plants (LCZ D)
        vlinder02         Large lowrise
        vlinder03          Open midrise
        vlinder04        Sparsely built
        vlinder05         Water (LCZ G)
        Name: lcz, dtype: object

        The LCZ are automatically added to the metadf as 'lcz' column.

        >>> dataset.metadf['lcz']
            name
        vlinder01         Low plants (LCZ D)
        vlinder02              Large lowrise
        vlinder03               Open midrise
        vlinder04             Sparsely built
        vlinder05              Water (LCZ G)
        vlinder06    Scattered Trees (LCZ B)
        vlinder07            Compact midrise
        vlinder08            Compact midrise
        vlinder09    Scattered Trees (LCZ B)
        vlinder10            Compact midrise
        vlinder11               Open lowrise
        vlinder12              Open highrise
        vlinder13            Compact midrise
        vlinder14         Low plants (LCZ D)
        vlinder15         Low plants (LCZ D)
        vlinder16              Water (LCZ G)
        vlinder17    Scattered Trees (LCZ B)
        vlinder18         Low plants (LCZ D)
        vlinder19            Compact midrise
        vlinder20            Compact midrise
        vlinder21             Sparsely built
        vlinder22         Low plants (LCZ D)
        vlinder23         Low plants (LCZ D)
        vlinder24        Dense Trees (LCZ A)
        vlinder25              Water (LCZ G)
        vlinder26               Open midrise
        vlinder27            Compact midrise
        vlinder28               Open lowrise
        Name: lcz, dtype: object
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

        A 'altitude' column will be added to the metadf, and series is returned.

        Returns
        -------
        altitude_series : pandas.Series()
            A series with the stationnames as index and the altitudes as values.

        Warning
        ---------
        This methods makes use of GEE API. Make sure that you have acces and
        user rights to use the GEE API.

        See Also
        --------
        connect_to_gee: Setup a new connection/credentials to the GEE service.
        get_lcz: Extract LCZ for all stations
        get_landcover: Extract landcoverfractions for all stations

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

        At this point there is no altitude information present in the metadata

        >>> dataset.metadf.head()
                         lat       lon        school                  geometry dataset_resolution                  dt_start                    dt_end
        name
        vlinder01  50.980438  3.815763         UGent  POINT (3.81576 50.98044)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder02  51.022379  3.709695         UGent   POINT (3.7097 51.02238)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder03  51.324583  4.952109   Heilig Graf  POINT (4.95211 51.32458)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder04  51.335522  4.934732   Heilig Graf  POINT (4.93473 51.33552)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder05  51.052655  3.675183  Sint-Barbara  POINT (3.67518 51.05266)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00

        We use the GEE API to extract the LCZ for all the stations present in
        the metadf (if coordinates are present).

        >>> alt_series = dataset.get_altitude()
        >>> alt_series.head()
        name
        vlinder01    12
        vlinder02     7
        vlinder03    30
        vlinder04    25
        vlinder05     0
        Name: altitude, dtype: int64

        The altitudes are automatically added to the metadf as 'altitude' column.

        >>> dataset.metadf['altitude']
        name
        vlinder01    12
        vlinder02     7
        vlinder03    30
        vlinder04    25
        vlinder05     0
        vlinder06     0
        vlinder07     7
        vlinder08     7
        vlinder09    19
        vlinder10    14
        vlinder11     6
        vlinder12     9
        vlinder13    10
        vlinder14     4
        vlinder15    41
        vlinder16     4
        vlinder17    83
        vlinder18    35
        vlinder19    75
        vlinder20    44
        vlinder21    19
        vlinder22     3
        vlinder23     1
        vlinder24    12
        vlinder25    12
        vlinder26    24
        vlinder27    12
        vlinder28     7
        Name: altitude, dtype: int64

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
        all stations. If an aggregation scheme is define, one can choose to
        aggregate the landcoverclasses.

        The landcover fractions will be added to the Dataset.metadf if overwrite
        is True. Presented as seperate columns where each column represent the
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
        This methods makes use of GEE API. Make sure that you have acces and
        user rights to use the GEE API.

        Warning
        ---------
        It can happen that for stations located on small islands, or close to
        the coast, the sea-mask is not used as a landcover fraction.

        See Also
        --------
        connect_to_gee: Setup a new connection/credentials to the GEE service.
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

        At this point there is no landcover information present in the metadata.

        >>> dataset.metadf.head()
                         lat       lon        school                  geometry dataset_resolution                  dt_start                    dt_end
        name
        vlinder01  50.980438  3.815763         UGent  POINT (3.81576 50.98044)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder02  51.022379  3.709695         UGent   POINT (3.7097 51.02238)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder03  51.324583  4.952109   Heilig Graf  POINT (4.95211 51.32458)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder04  51.335522  4.934732   Heilig Graf  POINT (4.93473 51.33552)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder05  51.052655  3.675183  Sint-Barbara  POINT (3.67518 51.05266)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00

        We use the GEE API to extract the landcover fractions for all the
        stations present in the metadf (if coordinates are present).

        >>> frac_df = dataset.get_landcover(buffers=[50, 100, 250, 500],
        ...                                 aggregate=False)
        >>> frac_df
                                   Grassland  Cropland  Tree cover  Built-up  Permanent water bodies  Herbaceous wetland  Bare / sparse vegetation  Shrubland
        name      buffer_radius
        vlinder01 50              0.545691  0.454309    0.000000  0.000000                0.000000                 0.0                       NaN        NaN
                  100             0.345583  0.636198    0.000000  0.018219                0.000000                 0.0                       NaN        NaN
                  250             0.318707  0.640718    0.004210  0.036365                0.000000                 0.0                  0.000000        NaN
                  500             0.390234  0.466117    0.049143  0.094263                0.000243                 0.0                  0.000000        0.0
        vlinder02 50              0.629257  0.000000    0.014283  0.356460                0.000000                 0.0                       NaN        NaN
        ...                            ...       ...         ...       ...                     ...                 ...                       ...        ...
        vlinder27 500             0.004314  0.000000    0.061079  0.923446                0.010837                 0.0                  0.000323        0.0
        vlinder28 50              0.049876  0.000000    0.159655  0.790469                0.000000                 0.0                       NaN        NaN
                  100             0.187910  0.000000    0.302041  0.510049                0.000000                 0.0                       NaN        NaN
                  250             0.128338  0.000000    0.593612  0.278050                0.000000                 0.0                  0.000000        NaN
                  500             0.115984  0.000162    0.548787  0.335067                0.000000                 0.0                  0.000000        0.0
        <BLANKLINE>
        [112 rows x 8 columns]


        The landcover fractions are automatically added to the metadf.


        >>> dataset.metadf.columns
        Index(['lat', 'lon', 'school', 'geometry', 'dataset_resolution', 'dt_start',
           'dt_end', 'Grassland_50m', 'Cropland_50m', 'Tree cover_50m',
           'Built-up_50m', 'Permanent water bodies_50m', 'Herbaceous wetland_50m',
           'Bare / sparse vegetation_50m', 'Shrubland_50m', 'Grassland_100m',
           'Cropland_100m', 'Tree cover_100m', 'Built-up_100m',
           'Permanent water bodies_100m', 'Herbaceous wetland_100m',
           'Bare / sparse vegetation_100m', 'Shrubland_100m', 'Grassland_250m',
           'Cropland_250m', 'Tree cover_250m', 'Built-up_250m',
           'Permanent water bodies_250m', 'Herbaceous wetland_250m',
           'Bare / sparse vegetation_250m', 'Shrubland_250m', 'Grassland_500m',
           'Cropland_500m', 'Tree cover_500m', 'Built-up_500m',
           'Permanent water bodies_500m', 'Herbaceous wetland_500m',
           'Bare / sparse vegetation_500m', 'Shrubland_500m'],
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


class MetobsDatasetGeeModelDataHandlingError(Exception):
    """Exception raised for errors in the Dataset - GEE modeldata interactions."""

    pass
