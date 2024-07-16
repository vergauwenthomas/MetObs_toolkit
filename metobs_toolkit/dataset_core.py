#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module contains the Dataset class and all its methods.

A Dataset holds all observations and is at the center of the
MetObs-toolkit.
"""

import os
import sys
import copy
from datetime import timedelta
import pytz
import logging
import pandas as pd
import numpy as np
import pickle

from metobs_toolkit.settings import Settings
from metobs_toolkit.data_import import import_data_from_csv, import_metadata_from_csv
from metobs_toolkit.template import Template

from metobs_toolkit.printing import print_dataset_info
from metobs_toolkit.landcover_functions import (
    connect_to_gee,
    lcz_extractor,
    height_extractor,
    lc_fractions_extractor,
    _validate_metadf,
)

# from metobs_toolkit.plotting_functions import (
#     geospatial_plot,
#     timeseries_plot,
#     # qc_stats_pie,
#     folium_plot,
#     add_stations_to_folium_map,
#     make_folium_html_plot,
# )

from metobs_toolkit.qc_checks import (
    # gross_value_check,
    # persistance_check,
    # repetitions_check,
    duplicate_timestamp_check,
    #     step_check,
    #     window_variation_check,
    invalid_input_check,
    #     toolkit_buddy_check,
    #     titan_buddy_check,
    #     titan_sct_resistant_check,
)


# from metobs_toolkit.qc_statistics import get_freq_statistics
from metobs_toolkit.writing_files import write_dataset_to_csv

# from metobs_toolkit.missingobs import Missingob_collection

from metobs_toolkit.gap import (
    # Gap,
    remove_gaps_from_obs,
    # remove_gaps_from_outliers,
    missing_timestamp_and_gap_check,
    get_gaps_indx_in_obs_space,
    get_station_gaps,
    # apply_interpolate_gaps,
    # make_gapfill_df,
    # apply_debias_era5_gapfill,
    # gaps_to_df,
)


from metobs_toolkit.df_helpers import (
    # multiindexdf_datetime_subsetting,
    fmt_datetime_argument,
    # init_multiindex,
    init_multiindexdf,
    init_triple_multiindexdf,
    metadf_to_gdf,
    conv_applied_qc_to_df,
    get_freqency_series,
    value_labeled_doubleidxdf_to_triple_idxdf,
    xs_save,
    concat_save,
)

from metobs_toolkit.obstypes import tlk_obstypes
from metobs_toolkit.obstypes import Obstype as Obstype_class


from metobs_toolkit.analysis import Analysis
from metobs_toolkit.modeldata import Modeldata
from metobs_toolkit.datasetbase import _DatasetBase


logger = logging.getLogger(__name__)


# =============================================================================
# Dataset class
# =============================================================================


class Dataset(_DatasetBase):
    """Objects holding observations and methods on observations."""

    def __init__(self):
        """Construct all the necessary attributes for Dataset object."""
        logger.info("Initialise dataset")

        _DatasetBase.__init__(self)  # holds df, metadf, obstypes and settings

        # Dataset with outlier observations
        self.outliersdf = init_triple_multiindexdf()

        self.missing_obs = None  # becomes a Missingob_collection after import
        self.gaps = None  # becomes a list of gaps

        self.gapfilldf = init_multiindexdf()
        self.missing_fill_df = init_multiindexdf()

        # Template for mapping data and metadata
        self.template = Template()

        self._istype = "Dataset"
        self._freqs = pd.Series(dtype=object)

        self._applied_qc = pd.DataFrame(columns=["obstype", "checkname"])
        self._qc_checked_obstypes = []  # list with qc-checked obstypes

    def __str__(self):
        """Represent as text."""
        if self.df.empty:
            if self._istype == "Dataset":
                return "Empty instance of a Dataset."
            elif self._istype == "Station":
                return "Empty instance of a Station."
            else:
                return "Empty instance of a Analysis."

        add_info = ""
        n_stations = self.df.index.get_level_values("name").unique().shape[0]
        n_obs_tot = self.df.shape[0]
        n_outl = self.outliersdf.shape[0]
        startdt = self.df.index.get_level_values("datetime").min()
        enddt = self.df.index.get_level_values("datetime").max()

        if (not self.metadf["lat"].isnull().all()) & (
            not self.metadf["lon"].isnull().all()
        ):
            add_info += "    *Coordinates are available for all stations."

        return (
            f"{self._istype} instance containing: \n \
    *{n_stations} stations \n \
    *{self.df.columns.to_list()} observation types \n \
    *{n_obs_tot} observation records \n \
    *{n_outl} records labeled as outliers \n \
    *{len(self.gaps)} gaps \n \
    *{self.missing_obs.series.shape[0]} missing observations \n \
    *records range: {startdt} --> {enddt} (total duration:  {enddt - startdt}) \n \
    *time zone of the records: {self.settings.time_settings['timezone']} \n "
            + add_info
        )

    def __repr__(self):
        """Info representation."""
        return self.__str__()

    def __add__(self, other, gapsize=None):
        """Addition of two Datasets."""
        # important !!!!!

        # the toolkit makes a new dataframe, and assumes the df from self and other
        # to be the input data.
        # This means that missing obs, gaps, invalid and duplicated records are
        # being looked for in the concatenation of both dataset, using their current
        # resolution !

        new = Dataset()
        self_obstypes = self.df.columns.to_list().copy()
        #  ---- df ----

        # check if observation of self are also in other
        assert all([(obs in other.df.columns) for obs in self_obstypes])
        # subset obstype of other to self
        other.df = other.df[self.df.columns.to_list()]

        # remove duplicate rows
        common_indexes = self.df.index.intersection(other.df.index)
        other.df = other.df.drop(common_indexes)

        # set new df
        new.df = concat_save([self.df, other.df])
        new.df = new.df.sort_index()

        #  ----- outliers df ---------

        other_outliers = other.outliersdf.reset_index()
        other_outliers = other_outliers[other_outliers["obstype"].isin(self_obstypes)]
        other_outliers = other_outliers.set_index(["name", "datetime", "obstype"])
        new.outliersdf = concat_save([self.outliersdf, other_outliers])
        new.outliersdf = new.outliersdf.sort_index()

        #  ------- Gaps -------------
        # Gaps have to be recaluculated using a frequency assumtion from the
        # combination of self.df and other.df, thus NOT the native frequency if
        # their is a coarsening allied on either of them.
        new.gaps = []

        # ---------- missing ---------
        # Missing observations have to be recaluculated using a frequency assumtion from the
        # combination of self.df and other.df, thus NOT the native frequency if
        # their is a coarsening allied on either of them.
        new.missing_obs = None

        # ---------- metadf -----------
        # Use the metadf from self and add new rows if they are present in other
        new.metadf = concat_save([self.metadf, other.metadf])
        new.metadf = new.metadf.drop_duplicates(keep="first")
        new.metadf = new.metadf.sort_index()

        # ------- specific attributes ----------

        # Template (units and descritpions) are taken from self
        new.template = self.template

        # Inherit Settings from self
        new.settings = copy.deepcopy(self.settings)

        # Applied qc:
        # TODO:  is this oke to do?
        new._applied_qc = pd.DataFrame(columns=["obstype", "checkname"])
        new._qc_checked_obstypes = []  # list with qc-checked obstypes

        # set init_dataframe to empty
        # NOTE: this is not necesarry but users will use this method when they
        # have a datafile that is to big. So storing and overloading a copy of
        # the very big datafile is invalid for these cases.
        new.input_df = pd.DataFrame()

        # ----- Apply IO QC ---------
        # Apply only checks that are relevant on records in between self and other
        # OR
        # that are dependand on the frequency (since the freq of the .df is used,
        # which is not the naitive frequency if coarsening is applied on either. )

        # missing and gap check
        if gapsize is None:
            gapsize = new.settings.gap["gaps_settings"]["gaps_finder"]["gapsize_n"]

        # note gapsize is now defined on the frequency of self
        new.missing_obs, new.gaps = missing_timestamp_and_gap_check(
            df=new.df,
            gapsize_n=self.settings.gap["gaps_settings"]["gaps_finder"]["gapsize_n"],
        )

        # duplicate check
        new.df, dup_outl_df = duplicate_timestamp_check(
            df=new.df,
            checks_info=new.settings.qc["qc_checks_info"],
            checks_settings=new.settings.qc["qc_check_settings"],
        )

        if not dup_outl_df.empty:
            new.update_outliersdf(add_to_outliersdf=dup_outl_df)

        # update the order and which qc is applied on which obstype
        checked_obstypes = list(self.obstypes.keys())

        checknames = ["duplicated_timestamp"]  # KEEP order

        new._applied_qc = concat_save(
            [
                new._applied_qc,
                conv_applied_qc_to_df(
                    obstypes=checked_obstypes, ordered_checknames=checknames
                ),
            ],
            ignore_index=True,
        )

        return new

    def show(self, show_all_settings=False, max_disp_n_gaps=5):
        """Show detailed information of the Dataset.

        A function to print out a detailed overview information about the Dataset.

        Parameters
        ----------
        show_all_settings : bool, optional
            If True all the settings are printed out. The default is False.
        max_disp_n_gaps: int, optional
            The maximum number of gaps to display detailed information of.
        Returns
        -------
        None.

        Examples
        --------
        .. code-block:: python

            >>> import metobs_toolkit
            >>>
            >>> #Create your Dataset
            >>> dataset = metobs_toolkit.Dataset() #empty Dataset

            >>>
            >>> #Add observations to the Dataset
            >>> dataset.update_settings(
            ...                         input_data_file=metobs_toolkit.demo_datafile,
            ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
            ...                         template_file=metobs_toolkit.demo_template,
            ...                         )
            >>> dataset.import_data_from_file()
            >>>
            >>> # Print out details
            >>> dataset.show()
            --------  General ---------
            ...

        """
        logger.info("Show basic info of dataset.")

        print_dataset_info(self, show_all_settings)

    def get_info(self, show_all_settings=False, max_disp_n_gaps=5):
        """Alias of show().

        A function to print out a detailed overview information about the Dataset.

        Parameters
        ----------
        show_all_settings : bool, optional
            If True all the settings are printed out. The default is False.
        max_disp_n_gaps: int, optional
            The maximum number of gaps to display detailed information of.

        Returns
        -------
        None.

        Examples
        --------
        .. code-block:: python

            >>> import metobs_toolkit
            >>>
            >>> #Create your Dataset
            >>> dataset = metobs_toolkit.Dataset() #empty Dataset

            >>>
            >>> #Add observations to the Dataset
            >>> dataset.update_settings(
            ...                         input_data_file=metobs_toolkit.demo_datafile,
            ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
            ...                         template_file=metobs_toolkit.demo_template,
            ...                         )
            >>> dataset.import_data_from_file()
            >>>
            >>> # Print out details
            >>> dataset.get_info()
            --------  General ---------
            ...

        """
        self.show(show_all_settings, max_disp_n_gaps)

    def save_dataset(
        self, outputfolder=None, filename="saved_dataset.pkl", overwrite=False
    ):
        """Save a Dataset instance to a (pickle) file.

        Parameters
        ----------
        outputfolder : str or None, optional
            The path to the folder to save the file. If None, the outputfolder
            from the Settings is used. The default is None.
        filename : str, optional
            The name of the output file. The default is 'saved_dataset.pkl'.
        overwrite : bool, optional
            If True, the target file will be overwritten if it exist. The
            default is False.

        Returns
        -------
        None.

        Examples
        --------
        .. code-block:: python

            >>> import metobs_toolkit
            >>> import os
            >>>
            >>> # Import data into a Dataset
            >>> dataset = metobs_toolkit.Dataset()
            >>> dataset.update_settings(
            ...            input_data_file=metobs_toolkit.demo_datafile,
            ...            input_metadata_file=metobs_toolkit.demo_metadatafile,
            ...            template_file=metobs_toolkit.demo_template)
            >>>
            >>> dataset.import_data_from_file()
            >>>
            >>> # Save dataset to a .pkl file
            >>> dataset.save_dataset(outputfolder=os.getcwd(),
            ...                     filename='your_saved_dataset.pkl',
            ...                     overwrite=True)
            Dataset saved in ...

        """
        # check if outputfolder is known and exists
        if outputfolder is None:
            outputfolder = self.settings.IO["output_folder"]
            assert (
                outputfolder is not None
            ), "No outputfolder is given, and no outputfolder is found in the settings."

        assert os.path.isdir(outputfolder), f"{outputfolder} is not a directory!"

        # check file extension in the filename:
        if filename[-4:] != ".pkl":
            filename += ".pkl"

        full_path = os.path.join(outputfolder, filename)
        if (os.path.isfile(full_path)) & overwrite:
            logger.info(f"The file {full_path} will be overwritten!")
            os.remove(full_path)

        # check if file exists
        assert not os.path.isfile(full_path), f"{full_path} is already a file!"

        with open(full_path, "wb") as outp:
            pickle.dump(self, outp, pickle.HIGHEST_PROTOCOL)

        print(f"Dataset saved in {full_path}")
        logger.info(f"Dataset saved in {full_path}")

    def import_dataset(self, folder_path=None, filename="saved_dataset.pkl"):
        """Import a Dataset instance from a (pickle) file.

        Parameters
        ----------
        folder_path : str or None, optional
            The path to the folder to save the file. If None, the outputfolder
            from the Settings is used. The default is None.
        filename : str, optional
            The name of the output file. The default is 'saved_dataset.pkl'.

        Returns
        -------
        metobs_toolkit.Dataset
            The Dataset instance.

        Examples
        --------
        .. code-block:: python

            import metobs_toolkit
            import os

            # Initialize an empty Dataset
            empty_dataset = metobs_toolkit.Dataset()

            # Import the dataset
            dataset=empty_dataset.import_dataset(folder_path=os.getcwd(),
                                                 filename='your_saved_dataset.pkl')

        """
        # check if folder_path is known and exists
        if folder_path is None:
            folder_path = self.settings.IO["output_folder"]
            assert (
                folder_path is not None
            ), "No folder_path is given, and no outputfolder is found in the settings."

        assert os.path.isdir(folder_path), f"{folder_path} is not a directory!"

        full_path = os.path.join(folder_path, filename)

        # check if file exists
        assert os.path.isfile(full_path), f"{full_path} does not exist."

        with open(full_path, "rb") as inp:
            dataset = pickle.load(inp)

        # convert metadf to a geodataframe (if coordinates are available)
        dataset.metadf = metadf_to_gdf(dataset.metadf)

        return dataset

    def add_new_observationtype(self, Obstype):
        """Add a new observation type to the known observation types.

        The observation can only be added if it is not already present in the
        knonw observation types. If that is the case that you probably need to
        use use the Dataset.add_new_unit() method.

        Parameters
        ----------

            The new Obstype to add.
        Returns
        -------
        None.

        Examples
        --------
        .. code-block:: python

            >>> import metobs_toolkit
            >>> co2_concentration = metobs_toolkit.Obstype(obsname='co2',
            ...                                            std_unit='ppm')
            >>> #add other units to it (if needed)
            >>> co2_concentration.add_unit(unit_name='ppb',
            ...                            conversion=['x / 1000'], #1 ppb = 0.001 ppm
            ...                           )
            >>> #Set a description
            >>> co2_concentration.set_description(desc='The CO2 concentration measured at 2m above surface')
            >>> #Add it to a Dataset
            >>> dataset = metobs_toolkit.Dataset()
            >>> dataset.add_new_observationtype(co2_concentration)
            >>> dataset.obstypes
            {'temp': Obstype instance of temp, 'humidity': Obstype instance of humidity, 'radiation_temp': Obstype instance of radiation_temp, 'pressure': Obstype instance of pressure, 'pressure_at_sea_level': Obstype instance of pressure_at_sea_level, 'precip': Obstype instance of precip, 'precip_sum': Obstype instance of precip_sum, 'wind_speed': Obstype instance of wind_speed, 'wind_gust': Obstype instance of wind_gust, 'wind_direction': Obstype instance of wind_direction, 'co2': Obstype instance of co2}
        """
        # Test if the obstype is of the correct class.
        if not isinstance(Obstype, Obstype_class):
            sys.exit(
                f"{Obstype} is not an instance of metobs_toolkit.obstypes.Obstype."
            )

        # Test if the obsname is already in use
        if Obstype.name in self.obstypes.keys():
            logger.warning(
                f"{Obstype.name} is already a known observation type: {self.obstypes[Obstype.name]}"
            )
            return

        # Update the known obstypes
        logger.info(f"Adding {Obstype} to the list of knonw observation types.")
        self.obstypes[Obstype.name] = Obstype

    def add_new_unit(self, obstype, new_unit, conversion_expression=[]):
        """Add a new unit to a known observation type.

        Parameters
        ----------
        obstype : str
            The observation type to add the new unit to.
        new_unit : str
            The new unit name.
        conversion_expression : list or str, optional
            The conversion expression to the standard unit of the observation
            type. The expression is a (list of) strings with simple algebraic
            operations, where x represent the value in the new unit, and the
            result is the value in the standard unit. Two examples for
            temperature (with a standard unit in Celsius):

                ["x - 273.15"] #if the new_unit is Kelvin
                ["x-32.0", "x/1.8"] #if the new unit is Farenheit

            The default is [].

        Returns
        -------
        None.

        Examples
        --------
        .. code-block:: python

            >>> import metobs_toolkit
            >>>
            >>> #Create your Dataset
            >>> dataset = metobs_toolkit.Dataset() #empty Dataset
            >>>
            >>> #Add new unit to a known obstype
            >>> dataset.add_new_unit(obstype = 'temp',
            ...                           new_unit= 'your_new_unit',
            ...                           conversion_expression = ['x+3', 'x * 2'])
            >>> # The conversion means: 1 [your_new_unit] = (1 + 3) * 2 [°C]
            >>> dataset.obstypes['temp'].get_info()
            temp observation with:
                 * standard unit: Celsius
                 * data column as None in None
                 * known units and aliases: {'Celsius': ['celsius', '°C', '°c', 'celcius', 'Celcius'], 'Kelvin': ['K', 'kelvin'], 'Farenheit': ['farenheit'], 'your_new_unit': []}
                 * description: 2m - temperature
                 * conversions to known units: {'Kelvin': ['x - 273.15'], 'Farenheit': ['x-32.0', 'x/1.8'], 'your_new_unit': ['x+3', 'x * 2']}
                 * originates from data column: None with None as native unit.
        """
        # test if observation is present
        if not obstype in self.obstypes.keys():
            logger.warning(f"{obstype} is not a known obstype! No unit can be added.")
            return

        # check if the unit is already present
        is_present = self.obstypes[obstype].test_if_unit_is_known(new_unit)
        if is_present:
            logger.info(
                f"{new_unit} is already a known unit of {self.obstypes[obstype]}"
            )
            return

        self.obstypes[obstype].add_unit(
            unit_name=new_unit, conversion=conversion_expression
        )

    def show_settings(self):
        """Show detailed information of the stored Settings.

        A function that prints out all the settings, structured per thematic.

        Returns
        -------
        None.

        Examples
        --------
        .. code-block:: python

            >>> import metobs_toolkit
            >>>
            >>> #Create your Dataset
            >>> dataset = metobs_toolkit.Dataset() #empty Dataset
            >>>
            >>> #Show default settings
            >>> dataset.show_settings()
            All settings:...

        """
        self.settings.show()

    def get_station(self, stationname):
        """Filter out one station of the Dataset.

        Extract a metobs_toolkit.Station object from the dataset by name.

        Parameters
        ----------
        stationname : string
            The name of the station.

        Returns
        -------
        metobs_toolkit.Station
            The station object.

        Examples
        --------
        .. code-block:: python

            >>> import metobs_toolkit
            >>>
            >>> #Create your Dataset
            >>> dataset = metobs_toolkit.Dataset() #empty Dataset

            >>>
            >>> #Add observations to the Dataset
            >>> dataset.update_settings(
            ...                         input_data_file=metobs_toolkit.demo_datafile,
            ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
            ...                         template_file=metobs_toolkit.demo_template,
            ...                         )
            >>> dataset.import_data_from_file()
            >>>
            >>> dataset.get_station('vlinder12')
            Station instance containing:
                 *1 stations
                 *['temp', 'humidity', 'wind_speed', 'wind_direction'] observation types
                 *4320 observation records
                 *0 records labeled as outliers
                 *0 gaps
                 *0 missing observations
                 *records range: 2022-09-01 00:00:00+00:00 --> 2022-09-15 23:55:00+00:00 (total duration:  14 days 23:55:00)
                 *time zone of the records: UTC
                 *Coordinates are available for all stations.

        """
        from metobs_toolkit.station import Station

        logger.info(f"Extract {stationname} from dataset.")

        # important: make sure all station attributes are of the same time as dataset.
        # so that all methods can be inherited.

        try:
            sta_df = self.df.xs(stationname, level="name", drop_level=False)
            sta_metadf = self.metadf.loc[stationname].to_frame().transpose()
            sta_metadf.index.name = "name"
        except KeyError:
            logger.warning(f"{stationname} not found in the dataset.")
            return None

        try:
            sta_outliers = self.outliersdf.xs(
                stationname, level="name", drop_level=False
            )
        except KeyError:
            sta_outliers = init_triple_multiindexdf()

        sta_gaps = get_station_gaps(self.gaps, stationname)
        sta_missingobs = self.missing_obs.get_station_missingobs(stationname)

        try:
            sta_gapfill = self.gapfilldf.xs(stationname, level="name", drop_level=False)
        except KeyError:
            sta_gapfill = init_multiindexdf()

        try:
            sta_missingfill = self.missing_fill_df.xs(
                stationname, level="name", drop_level=False
            )
        except KeyError:
            sta_missingfill = init_multiindexdf()

        return Station(
            name=stationname,
            df=sta_df,
            outliersdf=sta_outliers,
            gaps=sta_gaps,
            missing_obs=sta_missingobs,
            gapfilldf=sta_gapfill,
            missing_fill_df=sta_missingfill,
            metadf=sta_metadf,
            obstypes=self.obstypes,
            template=self.template,
            settings=self.settings,
            _qc_checked_obstypes=self._qc_checked_obstypes,
            _applied_qc=self._applied_qc,
        )

    # =============================================================================
    #   Gap Filling
    # =============================================================================
    def get_modeldata(
        self,
        modelname="ERA5_hourly",
        modeldata=None,
        obstype="temp",
        stations=None,
        startdt=None,
        enddt=None,
    ):
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

        Examples
        --------
        .. code-block:: python

            import metobs_toolkit

            # Import data into a Dataset
            dataset = metobs_toolkit.Dataset()
            dataset.update_settings(
                        input_data_file=metobs_toolkit.demo_datafile,
                        input_metadata_file=metobs_toolkit.demo_metadatafile,
                        template_file=metobs_toolkit.demo_template,
                        )
            dataset.import_data_from_file()

            # To limit data transfer, we define a short period
            import datetime

            tstart = datetime.datetime(2022, 9, 5)
            tend = datetime.datetime(2022, 9, 6)


            # Collect ERA5 2mT timeseries at your stations
            ERA5_data = dataset.get_modeldata(
                                    modelname="ERA5_hourly",
                                    modeldata=None,
                                    obstype="temp",
                                    stations=None,
                                    startdt=tstart,
                                    enddt=tend)

        """
        if modeldata is None:
            Modl = Modeldata(modelname)

        else:
            Modl = modeldata
            modelname = Modl.modelname

        # Filters

        if startdt is None:
            startdt = self.df.index.get_level_values("datetime").min()
        else:
            startdt = fmt_datetime_argument(
                startdt, self.settings.time_settings["timezone"]
            )

        if enddt is None:
            enddt = self.df.index.get_level_values("datetime").max()
        else:
            enddt = fmt_datetime_argument(
                enddt, self.settings.time_settings["timezone"]
            )

        # make shure bounds include required range
        Model_time_res = Modl.mapinfo[Modl.modelname]["time_res"]
        startdt = startdt.floor(Model_time_res)
        enddt = enddt.ceil(Model_time_res)

        if stations is not None:
            if isinstance(stations, str):
                metadf = self.metadf.loc[[stations]]
            if isinstance(stations, list):
                metadf = self.metadf.iloc[self.metadf.index.isin(stations)]
        else:
            metadf = self.metadf

        # Convert to UTC

        startdt_utc = startdt.astimezone(pytz.utc)
        enddt_utc = enddt.astimezone(pytz.utc)

        # fill modell with data
        if modelname == "ERA5_hourly":
            Modl.get_ERA5_data(
                metadf=metadf,
                startdt_utc=startdt_utc,
                enddt_utc=enddt_utc,
                obstypes=obstype,
            )

        else:
            Modl.get_gee_dataset_data(
                mapname=modelname,
                metadf=metadf,
                startdt_utc=startdt_utc,
                enddt_utc=enddt_utc,
                obstypes=obstype,
            )

        print(
            f"(When using the .set_model_from_csv() method, make sure the modelname of your Modeldata is {modelname})"
        )
        logger.info(
            f"(When using the .set_model_from_csv() method, make sure the modelname of your Modeldata is {modelname})"
        )
        return Modl

    def get_analysis(self, add_gapfilled_values=False):
        """Create an Analysis instance from the Dataframe.

        Parameters
        ----------
        add_gapfilled_values : bool, optional
            If True, all filled values (from gapfill and missing observation fill),
            are added to the analysis records aswell. The default is False.

        Returns
        -------
        metobs_toolkit.Analysis
            The Analysis instance of the Dataset.

        Examples
        --------

        .. code-block:: python

            >>> import metobs_toolkit
            >>>
            >>> # Import data into a Dataset
            >>> dataset = metobs_toolkit.Dataset()
            >>> dataset.update_settings(
            ...                         input_data_file=metobs_toolkit.demo_datafile,
            ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
            ...                         template_file=metobs_toolkit.demo_template,
            ...                         )
            >>> dataset.import_data_from_file()
            >>> dataset.coarsen_time_resolution(freq='1h')
            >>>
            >>> # Create an Analysis from the dataset
            >>> analysis = dataset.get_analysis()
            >>> analysis
            Analysis instance containing:
                 *28 stations
                 *['temp', 'humidity', 'wind_speed', 'wind_direction'] observation types
                 *10080 observation records
                 *Coordinates are available for all stations.
            <BLANKLINE>
                 *records range: 2022-09-01 00:00:00+00:00 --> 2022-09-15 23:00:00+00:00 (total duration:  14 days 23:00:00)     *Coordinates are available for all stations.
            <BLANKLINE>

        """
        # combine all to obsspace and include gapfill
        if add_gapfilled_values:
            mergedf = self.combine_all_to_obsspace()

            # gapsfilled labels
            gapfill_settings = self.settings.gap["gaps_fill_info"]
            gapfilllabels = [val for val in gapfill_settings["label"].values()]

            # missingfilled labels
            missingfill_settings = self.settings.missing_obs["missing_obs_fill_info"]
            missingfilllabels = [val for val in missingfill_settings["label"].values()]

            # get all labels
            fill_labels = gapfilllabels.copy()
            fill_labels.extend(missingfilllabels)
            fill_labels.append("ok")

            df = mergedf[mergedf["label"].isin(fill_labels)]
            df = df[["value"]]
            df = df.unstack(level="obstype")
            df = df.droplevel(level=0, axis=1)
        else:
            df = self.df

        return Analysis(
            obsdf=df,
            metadf=self.metadf,
            settings=self.settings,
            obstypes=self.obstypes,
        )

    def write_to_csv(
        self,
        obstype=None,
        filename=None,
        include_outliers=True,
        include_fill_values=True,
        add_final_labels=True,
        use_tlk_obsnames=True,
        overwrite_outliers_by_gaps_and_missing=True,
        seperate_metadata_file=True,
    ):
        """Write Dataset to a csv file.

        Write the dataset to a file where the observations, metadata and
        (if available) the quality labels per observation type are merged
        together.

        A final qualty control label for each
        quality-controlled-observation type can be added in the outputfile.

        The file will be written to the outputfolder specified in the settings.

        Parameters
        ----------
        obstype : string, optional
            Specify an observation type to subset all observations to. If None,
            all available observation types are written to file. The default is
            None.
        filename : string, optional
            The name of the output csv file. If none, a standard-filename
            is generated based on the period of data. The default is None.
        include_outliers : bool, optional
            If True, the outliers will be present in the csv file. The default is True.
        include_fill_values : bool, optional
            If True, the filled gap and missing observation values will be
            present in the csv file. The default is True.
        add_final_labels : bool, optional
            If True, a column is added containing the final label of an observation. The default is True.
        use_tlk_obsnames : bool, optional
            If True, the standard naming of the metobs_toolkit is used, else
            the original names for obstypes is used. The default is True.
        overwrite_outliers_by_gaps_and_missing : bool, optional
            If the gaps and missing observations are updated using outliers,
            interpret these records as gaps/missing outliers if True. Else these
            will be interpreted as outliers. The default is True.
        seperate_metadata_file : bool, optional
            If true, the metadat is written to a seperate file, else the metadata
            is merged to the observation in one file. The default is True.
        Returns
        -------
        None.

        Examples
        --------
        .. code-block:: python

            >>> import metobs_toolkit
            >>> import os
            >>>
            >>> # Import data into a Dataset
            >>> dataset = metobs_toolkit.Dataset()
            >>> dataset.update_settings(
            ...            input_data_file=metobs_toolkit.demo_datafile,
            ...            input_metadata_file=metobs_toolkit.demo_metadatafile,
            ...            template_file=metobs_toolkit.demo_template)
            >>>
            >>> dataset.import_data_from_file()
            >>>
            >>> # Save dataset to a .csv file
            >>> dataset.update_settings(output_folder = os.getcwd())
            >>> dataset.write_to_csv(filename='your_saved_table.csv')
            write metadata to file: ...
            write dataset to file: ...

        """
        logger.info("Writing the dataset to a csv file")

        assert (
            not self.settings.IO["output_folder"] is None
        ), "Specify Settings.output_folder in order to export a csv."

        assert os.path.isdir(
            self.settings.IO["output_folder"]
        ), f'The outputfolder: \
            {self.settings.IO["output_folder"]} is not found. '

        # combine all dataframes
        mergedf = self.combine_all_to_obsspace(
            overwrite_outliers_by_gaps_and_missing=overwrite_outliers_by_gaps_and_missing
        )  # with outliers
        # Unstack mergedf
        # remove duplicates
        mergedf = mergedf[~mergedf.index.duplicated(keep="first")]

        # drop outliers if required
        if not include_outliers:
            outlier_labels = [
                var["outlier_flag"] for var in self.settings.qc["qc_checks_info"]
            ]
            mergedf = mergedf[~mergedf["label"].isin(outlier_labels)]

        # drop fill values if required
        if not include_fill_values:
            fill_labels = [
                "gap fill",
                "missing observation fill",
            ]  # toolkit representation labels
            mergedf = mergedf[~mergedf["toolkit_representation"].isin(fill_labels)]

        if obstype is not None:
            mergedf = xs_save(mergedf, obstype, level="obstype", drop_level=False)

        # Map obstypes columns
        if not use_tlk_obsnames:
            mapper = {
                col: self.obstypes[col].get_orig_name() for col in self.obstypes.keys()
            }
            mergedf = mergedf.reset_index()
            mergedf["new_names"] = mergedf["obstype"].map(mapper)
            mergedf = mergedf.drop(columns=["obstype"])
            mergedf = mergedf.rename(columns={"new_names": "obstype"})
            mergedf = mergedf.set_index(["name", "datetime", "obstype"])

        mergedf = mergedf.unstack("obstype")

        # to one level for the columns
        mergedf.columns = [" : ".join(col).strip() for col in mergedf.columns.values]

        # columns to write
        write_dataset_to_csv(
            df=mergedf,
            metadf=self.metadf,
            filename=filename,
            outputfolder=self.settings.IO["output_folder"],
            seperate_metadata_file=seperate_metadata_file,
        )

    # =============================================================================
    #     Quality control
    # =============================================================================

    def combine_all_to_obsspace(
        self,
        repr_outl_as_nan=False,
        overwrite_outliers_by_gaps_and_missing=True,
    ):
        """Make one dataframe with all observations and their labels.

        Combine all observations, outliers, missing observations and gaps into
        one Dataframe. All observation types are combined an a label is added
        in a serperate column.

        When gaps and missing records are updated from outliers one has to choice
        to represent these records as outliers or gaps. There can not be duplicates
        in the return dataframe.

        By default the observation values of the outliers are saved, one can
        choice to use these values or NaN's.
        following checks!

        Parameters
        ----------
        repr_outl_as_nan : bool, optional
            If True, Nan's are use for the values of the outliers. The
            default is False.
        overwrite_outliers_by_gaps_and_missing : Bool, optional
            If True, records that are labeld as gap/missing and outlier are
            labeled as gaps/missing. This has only effect when the gaps/missing
            observations are updated from the outliers. The default is True.

         Returns
         ---------
         combdf : pandas.DataFrame()
            A dataframe containing a continious time resolution of records, where each
            record is labeld.

        Examples
        --------
        .. code-block:: python

            >>> import metobs_toolkit
            >>>
            >>> # Import data into a Dataset
            >>> dataset = metobs_toolkit.Dataset()
            >>> dataset.update_settings(
            ...                         input_data_file=metobs_toolkit.demo_datafile,
            ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
            ...                         template_file=metobs_toolkit.demo_template,
            ...                         )
            >>> dataset.import_data_from_file()
            >>> dataset.coarsen_time_resolution(freq='1h')
            >>>
            >>> # Apply quality control on the temperature observations
            >>> dataset.apply_quality_control(obstype='temp') #Using the default QC settings
            >>> dataset
            Dataset instance containing:
                 *28 stations
                 *['temp', 'humidity', 'wind_speed', 'wind_direction'] observation types
                 *10080 observation records
                 *1676 records labeled as outliers
                 *0 gaps
                 *3 missing observations
                 *records range: 2022-09-01 00:00:00+00:00 --> 2022-09-15 23:00:00+00:00 (total duration:  14 days 23:00:00)
                 *time zone of the records: UTC
                 *Coordinates are available for all stations.
            >>>
            >>> # Combine all records to one dataframe in Observation-resolution
            >>> overview_df = dataset.combine_all_to_obsspace()
            >>> overview_df.head(12)
                                                                    value  ... toolkit_representation
            name      datetime                  obstype                    ...
            vlinder01 2022-09-01 00:00:00+00:00 humidity        65.000000  ...            observation
                                                temp            18.800000  ...            observation
                                                wind_direction  65.000000  ...            observation
                                                wind_speed       1.555556  ...            observation
                      2022-09-01 01:00:00+00:00 humidity        65.000000  ...            observation
                                                temp            18.400000  ...            observation
                                                wind_direction  55.000000  ...            observation
                                                wind_speed       1.416667  ...            observation
                      2022-09-01 02:00:00+00:00 humidity        68.000000  ...            observation
                                                temp            17.100000  ...            observation
                                                wind_direction  45.000000  ...            observation
                                                wind_speed       1.583333  ...            observation
            <BLANKLINE>
            [12 rows x 3 columns]

        """
        # TODO: label values from settings not hardcoding

        # TODO: use the repr_outl_as_nan argumenten here
        # =============================================================================
        # Stack observations and outliers
        # =============================================================================
        df = self.df
        # better save than sorry
        present_obstypes = list(self.obstypes.keys())
        df = df[present_obstypes]

        # to tripple index
        df = (
            df.stack(future_stack=True)
            .reset_index()
            .rename(columns={"level_2": "obstype", 0: "value"})
            .set_index(["name", "datetime", "obstype"])
        )

        df["label"] = "ok"
        df["toolkit_representation"] = "observation"

        # outliers
        outliersdf = self.outliersdf.copy()
        outliersdf["toolkit_representation"] = "outlier"

        # Careful! Some outliers exist on inport frequency (duplicated, invalid)
        # So only use the outliers for which station-datetime-obstype are present in the
        # dataset.df
        outliersdf = outliersdf[outliersdf.index.isin(df.index)]

        # remove outliers from the observations
        df = df[~df.index.isin(outliersdf.index)]

        # =============================================================================
        # Stack gaps
        # =============================================================================
        # add gapfill and remove the filled records from gaps
        gapsfilldf = self.gapfilldf.copy()

        # to triple index
        gapsfilldf = value_labeled_doubleidxdf_to_triple_idxdf(
            gapsfilldf, known_obstypes=list(self.obstypes.keys())
        )
        gapsfilldf["toolkit_representation"] = "gap fill"

        gapsidx = get_gaps_indx_in_obs_space(
            gapslist=self.gaps,
            obsdf=self.df,
            outliersdf=self.outliersdf,
            resolutionseries=self.metadf["dataset_resolution"],
        )

        gapsdf = pd.DataFrame(index=gapsidx, columns=present_obstypes)
        gapsdf = (
            gapsdf.stack(future_stack=True)
            .reset_index()
            .rename(columns={"level_2": "obstype", 0: "value"})
            .set_index(["name", "datetime", "obstype"])
        )

        gapsdf["label"] = self.settings.gap["gaps_info"]["gap"]["outlier_flag"]
        gapsdf["toolkit_representation"] = "gap"

        # Remove gaps from df
        df = df[~df.index.isin(gapsdf.index)]

        if overwrite_outliers_by_gaps_and_missing:
            outliersdf = outliersdf.drop(index=gapsdf.index, errors="ignore")

        # Remove gapfill values records from the gaps
        gapsdf = gapsdf.drop(index=gapsfilldf.index)

        # =============================================================================
        # Stack missing
        # =============================================================================
        missingfilldf = self.missing_fill_df.copy()
        missingfilldf = value_labeled_doubleidxdf_to_triple_idxdf(
            missingfilldf, known_obstypes=list(self.obstypes.keys())
        )
        missingfilldf["toolkit_representation"] = "missing observation fill"

        # add missing observations if they occure in observation space
        missingidx = self.missing_obs.get_missing_indx_in_obs_space(
            self.df, self.metadf["dataset_resolution"]
        )

        missingdf = pd.DataFrame(index=missingidx, columns=present_obstypes)

        missingdf = (
            missingdf.stack(future_stack=True)
            .reset_index()
            .rename(columns={"level_2": "obstype", 0: "value"})
            .set_index(["name", "datetime", "obstype"])
        )

        missingdf["label"] = self.settings.gap["gaps_info"]["missing_timestamp"][
            "outlier_flag"
        ]
        missingdf["toolkit_representation"] = "missing observation"

        # Remove missing from df
        df = df[~df.index.isin(missingdf.index)]

        if overwrite_outliers_by_gaps_and_missing:
            outliersdf = outliersdf.drop(index=missingdf.index, errors="ignore")

        # Remove missingfill values records from the missing
        missingdf = missingdf.drop(index=missingfilldf.index)

        # =============================================================================
        # combine all
        # =============================================================================

        combdf = concat_save(
            [df, outliersdf, gapsdf, gapsfilldf, missingdf, missingfilldf]
        ).sort_index()
        combdf.index.names = ["name", "datetime", "obstype"]
        # To be shure?
        combdf = combdf[~combdf.index.duplicated(keep="first")]
        return combdf

    def update_outliersdf(self, add_to_outliersdf):
        """Update the outliersdf attribute."""
        self.outliersdf = concat_save([self.outliersdf, add_to_outliersdf])

    def coarsen_time_resolution(
        self, origin=None, origin_tz=None, freq=None, method=None, limit=None
    ):
        """Resample the observations to coarser timeresolution.

        The assumed dataset resolution (stored in the metadf attribute) will be
        updated.

        Parameters
        ----------
        origin : datetime.datetime, optional
            Define the origin (first timestamp) for the obervations. The origin
            is timezone naive, and is assumed to have the same timezone as the
            obervations. If None, the earliest occuring timestamp is used as
            origin. The default is None.
        origin_tz : str, optional
            Timezone string of the input observations. Element of
            pytz.all_timezones. If None, the timezone from the settings is
            used. The default is None.
        freq : DateOffset, Timedelta or str, optional
            The offset string or object representing target conversion.
            Ex: '15min' is 15 minutes, '1h', is one hour. If None, the target time
            resolution of the dataset.settings is used. The default is None.
        method : 'nearest' or 'bfill', optional
            Method to apply for the resampling. If None, the resample method of
            the dataset.settings is used. The default is None.
        limit : int, optional
            Limit of how many values to fill with one original observations. If
            None, the target limit of the dataset.settings is used. The default
            is None.

        Returns
        -------
        None.

        Examples
        --------
        .. code-block:: python

            >>> import metobs_toolkit
            >>>
            >>> # Import data into a Dataset
            >>> dataset = metobs_toolkit.Dataset()
            >>> dataset.update_settings(
            ...                         input_data_file=metobs_toolkit.demo_datafile,
            ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
            ...                         template_file=metobs_toolkit.demo_template,
            ...                         )
            >>> dataset.import_data_from_file()
            >>> dataset.coarsen_time_resolution(freq='15min') #to 15 minutes resolution
            >>> dataset.df[['temp', 'humidity']].head()
                                                 temp  humidity
            name      datetime
            vlinder01 2022-09-01 00:00:00+00:00  18.8      65
                      2022-09-01 00:15:00+00:00  18.7      65
                      2022-09-01 00:30:00+00:00  18.7      65
                      2022-09-01 00:45:00+00:00  18.6      65
                      2022-09-01 01:00:00+00:00  18.4      65

        """
        if freq is None:
            freq = self.settings.time_settings["target_time_res"]
        if method is None:
            method = self.settings.time_settings["resample_method"]
        if limit is None:
            limit = int(self.settings.time_settings["resample_limit"])
        if origin_tz is None:
            origin_tz = self.settings.time_settings["timezone"]

        logger.info(
            f"Coarsening the timeresolution to {freq} using \
                    the {method}-method (with limit={limit})."
        )

        # test if coarsening the resolution is valid for the dataset
        # 1. If resolution-dep-qc is applied --> coarsening is not valid and will result in a broken dataset

        if (
            self._applied_qc[
                ~self._applied_qc["checkname"].isin(
                    ["duplicated_timestamp", "invalid_input"]
                )
            ].shape[0]
            > 0
        ):
            logger.warning(
                "Coarsening time resolution is not possible because quality control checks that are resolution depening are already performed on the Dataset."
            )
            logger.info(
                "(Apply coarsening_time_resolution BEFORE applying quality control.)"
            )
            return

        # TODO: implement buffer method
        df = self.df.reset_index()

        if origin is None:
            # find earlyest timestamp, if it is on the hour, use it else use the following hour
            tstart = df["datetime"].min()

            if tstart.minute != 0 or tstart.second != 0 or tstart.microsecond != 0:
                # Round up to nearest hour
                tstart = tstart.ceil(freq=freq)
        else:
            origin_tz_aware = pytz.timezone(origin_tz).localize(origin)
            tstart = origin_tz_aware.astimezone(
                pytz.timezone(self.settings.time_settings["timezone"])
            )

        # Coarsen timeresolution

        if method == "nearest":
            df = (
                df.set_index("datetime")
                .groupby("name")
                .resample(freq, origin=tstart)
                .nearest(limit=limit)
            )

        elif method == "bfill":
            df = (
                df.set_index("datetime")
                .groupby("name")
                .resample(freq, origin=tstart)
                .bfill(limit=limit)
            )

        else:
            logger.warning(f"The coarsening method: {method}, is not implemented yet.")
            df = df.set_index(["name", "datetime"])

        if "name" in df.columns:
            df = df.drop(columns=["name"])

        # Update resolution info in metadf
        self.metadf["dataset_resolution"] = pd.to_timedelta(freq)
        # update df
        self.df = df

        # Remove gaps and missing from the observatios
        # most gaps and missing are already removed but when increasing timeres,
        # some records should be removed as well.
        self.df = remove_gaps_from_obs(gaplist=self.gaps, obsdf=self.df)
        self.df = self.missing_obs.remove_missing_from_obs(obsdf=self.df)

    def sync_observations(
        self,
        tolerance,
        verbose=True,
        _force_resolution_minutes=None,
        _drop_target_nan_dt=False,
    ):
        """Simplify and syncronize the observation timestamps.

        To simplify the resolution (per station), a tolerance is use to shift timestamps. The tolerance indicates the
        maximum translation in time that can be applied to an observation.

        The sycronisation tries to group stations that have an equal simplified resolution, and syncronize them. The origin
        of the sycronized timestamps will be set to round hours, round 10-minutes or round-5 minutes if possible given the tolerance.

        The observations present in the input file are used.

        After syncronization, the IO outliers, missing observations and gaps are recomputed.

        Parameters
        ----------
        tolerance :  Timedelta or str
            The tolerance string or object representing the maximum translation in time.
            Ex: '5min' is 5 minutes, '1h', is one hour.
        verbose : bool, optional
            If True, a dataframe illustrating the mapping from original datetimes to simplified and syncronized is returned. The default is True.
        _drop_target_nan_dt : bool, optional
            If record has no target datetime, the datetimes will be listed as Nat. To remove them,
            set this to True. Default is False.
        _force_resolution_minutes : bool, optional
            force the resolution estimate to this frequency in minutes. If None, the frequency is estimated. The default is None.
        Note
        --------
        Keep in mind that this method will overwrite the df, outliersdf, missing timestamps and gaps.

        Note
        --------
        Because the used observations are from the input file, previously coarsend timeresolutions are ignored.


        Returns
        -------
        pandas.DataFrame (if verbose is True)
            A dataframe containing the original observations with original timestamps and the corresponding target timestamps.

        """
        # get columns pressent in metadf, because the input df can have columns
        # that does not have to be mapped to the toolkit

        assert (
            not self.input_df.empty
        ), "To syncronize a dataset, the (pure) input dataframe cannot be empty."

        init_meta_cols = self.metadf.columns.copy()
        df = self.input_df

        self.df = init_multiindexdf()
        self.outliersdf = init_triple_multiindexdf()
        self.gapfilldf = init_multiindexdf()
        self.missing_obs = None
        self.gaps = None

        # find simplified resolution
        if _force_resolution_minutes is None:
            simplified_resolution = get_freqency_series(
                df=df,
                method="median",
                simplify=True,
                max_simplify_error=tolerance,
            )
        else:
            if isinstance(_force_resolution_minutes, list):
                # TODO
                print(
                    "foce resolution minutes as a list is not implemented yet, sorry."
                )
            else:
                stations = self.metadf.index
                freq_series = pd.Series(
                    index=stations,
                    data=[timedelta(minutes=float(_force_resolution_minutes))]
                    * len(stations),
                )
                simplified_resolution = freq_series

        logger.debug(f"Syncronizing to these resolutions: {simplified_resolution}")

        occuring_resolutions = simplified_resolution.unique()

        df = df.reset_index()

        def find_simple_origin(tstart, tolerance):
            if tstart.minute == 0 and tstart.second == 0 and tstart.microsecond == 0:
                return tstart  # already a round hour

            # try converting to a round hour
            tstart_round_hour = tstart.round("60min")
            if abs(tstart - tstart_round_hour) <= pd.to_timedelta(tolerance):
                return tstart_round_hour

            # try converting to a tenfold in minutes
            tstart_round_tenfold = tstart.round("10min")
            if abs(tstart - tstart_round_tenfold) <= pd.to_timedelta(tolerance):
                return tstart_round_tenfold

            # try converting to a fivefold in minutes
            tstart_round_fivefold = tstart.round("5min")

            if abs(tstart - tstart_round_fivefold) <= pd.to_timedelta(tolerance):
                return tstart_round_fivefold

            # no suitable conversion found
            return tstart

        merged_df = pd.DataFrame()
        _total_verbose_df = pd.DataFrame()
        for occur_res in occuring_resolutions:
            group_stations = simplified_resolution[
                simplified_resolution == occur_res
            ].index.to_list()
            logger.info(
                f" Grouping stations with simplified resolution of {pd.to_timedelta(occur_res)}: {group_stations}"
            )
            groupdf = df[df["name"].isin(group_stations)]

            tstart = groupdf["datetime"].min()
            tend = groupdf["datetime"].max()

            # find a good origin point
            origin = find_simple_origin(tstart=tstart, tolerance=tolerance)

            # Create records index
            target_records = pd.date_range(
                start=origin, end=tend, freq=pd.Timedelta(occur_res)
            ).to_series()

            target_records.name = "target_datetime"
            # convert records to new target records, station per station

            for sta in group_stations:
                stadf = groupdf[groupdf["name"] == sta]
                # Drop all nan values! these will be added later from the outliersdf
                stadf = stadf.set_index(["name", "datetime"])

                # drop all records per statiotion for which there are no obsecvations
                present_obs = list(self.obstypes.keys())

                stadf = stadf.loc[stadf[present_obs].dropna(axis=0, how="all").index]

                stadf = stadf.reset_index()

                mergedstadf = pd.merge_asof(
                    left=stadf.sort_values("datetime"),
                    right=target_records.to_frame(),
                    right_on="target_datetime",
                    left_on="datetime",
                    direction="nearest",
                    tolerance=pd.Timedelta(tolerance),
                )
                if _drop_target_nan_dt:
                    mergedstadf = mergedstadf.dropna(subset="target_datetime")
                # possibility 1: record is mapped crrectly
                correct_mapped = mergedstadf[~mergedstadf["target_datetime"].isnull()]

                # possibility2: records that ar not mapped to target
                # not_mapped_records =mergedstadf[mergedstadf['target_datetime'].isnull()]

                # possibilyt 3 : no suitable candidates found for the target
                # these will be cached by the missing and gap check
                # no_record_candidates = target_records[~target_records.isin(mergedstadf['target_datetime'])].values

                merged_df = concat_save([merged_df, correct_mapped])

                if verbose:
                    _total_verbose_df = concat_save([_total_verbose_df, mergedstadf])

        # overwrite the df with the synced observations
        merged_df = (
            merged_df.rename(
                columns={
                    "datetime": "original_datetime",
                    "target_datetime": "datetime",
                }
            )
            .set_index(["name", "datetime"])
            .drop(["original_datetime"], errors="ignore", axis=1)
            .sort_index()
        )
        # self.df = merged_df

        # Recompute the dataset attributes, apply qc, gap and missing searches, etc.
        self._construct_dataset(
            df=merged_df,
            freq_estimation_method="highest",
            freq_estimation_simplify=False,
            freq_estimation_simplify_error=None,
            fixed_freq_series=simplified_resolution,
            update_full_metadf=False,
        )  # Do not overwrite full metadf, only the frequencies

        self.metadf = self.metadf[
            [col for col in self.metadf.columns if col in init_meta_cols]
        ]

        if verbose:
            _total_verbose_df = _total_verbose_df.rename(
                columns={
                    "datetime": "original_datetime",
                    "target_datetime": "datetime",
                }
            ).set_index(["name", "datetime"])
            return _total_verbose_df

    def import_data_from_file(
        self,
        input_data_file=None,
        input_metadata_file=None,
        template_file=None,
        freq_estimation_method=None,
        freq_estimation_simplify=None,
        freq_estimation_simplify_error=None,
        kwargs_data_read={},
        kwargs_metadata_read={},
    ):
        """Read observations from a csv file.

        The paths (data, metadata and template) are stored in the settings if
        Dataset.update_settings() is called on this object. These paths can be
        updated by adding them as argument to this method.

        The input data (and metadata) are interpreted by using a template
        (json file).

        An estimation of the observational frequency is made per station. This is used
        to find missing observations and gaps.

        The Dataset attributes are set and the following checks are executed:
                * Duplicate check
                * Invalid input check
                * Find missing observations
                * Find gaps


        Parameters
        ----------
        input_data_file : string, optional
            Path to the input data file with observations. If None, the input
            data path in the settings is used.
        input_metadata_file : string, optional
            Path to the input metadata file. If None, the input metadata path
            in the settings is used.
        template_file : string, optional
            Path to the template (json) file to be used on the observations
            and metadata. If None, the template path in the settings is used.
        freq_estimation_method : 'highest' or 'median', optional
            Select wich method to use for the frequency estimation. If
            'highest', the highest apearing frequency is used. If 'median', the
            median of the apearing frequencies is used. If None, the method
            stored in the
            Dataset.settings.time_settings['freq_estimation_method'] is used.
            The default is None.
        freq_estimation_simplify : bool, optional
            If True, the likely frequency is converted to round hours, or round minutes.
            The "freq_estimation_simplify_error' is used as a constrain. If the constrain is not met,
            the simplification is not performed. If None, the method
            stored in the
            Dataset.settings.time_settings['freq_estimation_simplify'] is used.
            The default is None.
        freq_estimation_simplify_error : Timedelta or str, optional
            The tolerance string or object representing the maximum translation in time to form a simplified frequency estimation.
            Ex: '5min' is 5 minutes, '1h', is one hour. If None, the method
            stored in the
            Dataset.settings.time_settings['freq_estimation_simplify_error'] is
            used. The default is None.
        kwargs_data_read : dict, optional
            Keyword arguments collected in a dictionary to pass to the
            pandas.read_csv() function on the data file. The default is {}.
        kwargs_metadata_read : dict, optional
            Keyword arguments collected in a dictionary to pass to the
            pandas.read_csv() function on the metadata file. The default is {}.

        Note
        --------
        In pracktice, the default arguments will be sufficient for most applications.

        Note
        --------
        If options are present in the template, these will have priority over the arguments of this function.

        Warning
        ---------
        All CSV data files must be in *UTF-8 encoding*. For most CSV files,
        this condition is already met. To make sure, in Microsoft Excel (or
        similar), you can specify to export as **`CSV UTF-8`**. If you
        encounter an error, mentioning a `"/ueff..."` tag in a CSV file, it is
        often solved by converting the CSV to UTF-8.

        Returns
        -------
        None.

        Examples
        --------
        .. code-block:: python

            >>> import metobs_toolkit
            >>>
            >>> # Import data into a Dataset
            >>> dataset = metobs_toolkit.Dataset()
            >>> dataset.update_settings(
            ...                         input_data_file=metobs_toolkit.demo_datafile,
            ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
            ...                         template_file=metobs_toolkit.demo_template,
            ...                         )
            >>> dataset.import_data_from_file()

        """
        logger.info(f'Importing data from file: {self.settings.IO["input_data_file"]}')

        # Update paths to the input files, if given.
        if input_data_file is not None:
            self.update_settings(input_data_file=input_data_file)
        if input_metadata_file is not None:
            self.update_settings(input_metadata_file=input_metadata_file)
        if template_file is not None:
            self.update_settings(template_file=template_file)

        if freq_estimation_method is None:
            freq_estimation_method = self.settings.time_settings[
                "freq_estimation_method"
            ]
        if freq_estimation_simplify is None:
            freq_estimation_simplify = self.settings.time_settings[
                "freq_estimation_simplify"
            ]
        if freq_estimation_simplify_error is None:
            freq_estimation_simplify_error = self.settings.time_settings[
                "freq_estimation_simplify_error"
            ]

        assert self.settings.templatefile is not None, "No templatefile is specified."

        # Read template
        self.template.read_template_from_file(jsonpath=self.settings.templatefile)

        # overload the timezone from template to the settings
        if not self.template._get_tz() is None:
            self.update_timezone(self.template.get_tz())
            logger.info(
                f"Set timezone = {self.template.get_tz()} from options in template."
            )

        # Read observations into pandas dataframe
        df = import_data_from_csv(
            input_file=self.settings.IO["input_data_file"],
            template=self.template,
            known_obstypes=list(self.obstypes.keys()),
            kwargs_data_read=kwargs_data_read,
        )
        # Set timezone information
        df.index = df.index.tz_localize(
            tz=self.settings.time_settings["timezone"],
            ambiguous="infer",
            nonexistent="shift_forward",
        )

        # drop Nat datetimes if present
        df = df.loc[pd.notnull(df.index)]

        logger.debug(
            f'Data from {self.settings.IO["input_data_file"]} \
                     imported to dataframe {df.head()}.'
        )

        if self.settings.IO["input_metadata_file"] is None:
            logger.warning(
                "No metadata file is defined,\
                    no metadata attributes can be set!"
            )

        else:
            logger.info(
                f'Importing metadata from file: {self.settings.IO["input_metadata_file"]}'
            )
            meta_df = import_metadata_from_csv(
                input_file=self.settings.IO["input_metadata_file"],
                template=self.template,
                kwargs_metadata_read=kwargs_metadata_read,
            )
            # in dataset of one station
            if self.template._is_data_single_station():
                # logger.warning("No station names find in the observations!")

                # If there is ONE name in the metadf, than we use that name for
                # the df, else we use the default name
                if ("name" in meta_df.columns) & (meta_df.shape[0] == 1):
                    name = meta_df["name"].iloc[0]
                    df["name"] = name
                    logger.warning(
                        f"One stationname found in the metadata: {name}, this name is used for the data."
                    )
                else:
                    df["name"] = str(self.settings.app["default_name"])
                    # for later merging, we add the name column with the default
                    # also in the metadf
                    meta_df["name"] = str(self.settings.app["default_name"])
                    logger.warning(
                        f'Assume the dataset is for ONE station with the \
                        default name: {self.settings.app["default_name"]}.'
                    )

            # merge additional metadata to observations
            logger.debug(f"Head of data file, before merge: {df.head()}")
            logger.debug(f"Head of metadata file, before merge: {meta_df.head()}")

            meta_cols = [
                colname for colname in meta_df.columns if not colname.startswith("_")
            ]
            additional_meta_cols = list(set(meta_cols).difference(df.columns))

            if bool(additional_meta_cols):
                logger.debug(
                    f"Merging metadata ({additional_meta_cols}) to dataset data by name."
                )
                additional_meta_cols.append("name")  # merging on name
                # merge deletes datetime index somehow? so add it back.
                df_index = df.index
                df = df.merge(
                    right=meta_df[additional_meta_cols], how="left", on="name"
                )
                df.index = df_index

        # update dataset object

        # Remove stations whith only one observation (no freq estimation)
        station_counts = df["name"].value_counts()
        issue_station = station_counts[station_counts < 2].index.to_list()
        logger.warning(
            f"These stations will be removed because of only having one record: {issue_station}"
        )
        df = df[~df["name"].isin(issue_station)]

        # convert dataframe to multiindex (datetime - name)
        df = df.set_index(["name", df.index])

        # Sort by name and then by datetime (to avoid negative freq)
        df = df.sort_index(level=["name", "datetime"])

        # dataframe with all data of input file
        self.input_df = df.sort_index(level=["name", "datetime"])
        # Construct all attributes of the Dataset
        self._construct_dataset(
            df=df,
            freq_estimation_method=freq_estimation_method,
            freq_estimation_simplify=freq_estimation_simplify,
            freq_estimation_simplify_error=freq_estimation_simplify_error,
        )

    def _construct_dataset(
        self,
        df,
        freq_estimation_method,
        freq_estimation_simplify,
        freq_estimation_simplify_error,
        fixed_freq_series=None,
        update_full_metadf=True,
    ):
        """Construct the Dataset class from a IO dataframe.

        The df, metadf, outliersdf, gaps, missing timestamps and observationtypes attributes are set.


        The observations are converted to the toolkit standard units if possible.

        Qc on IO is applied (duplicated check and invalid check) + gaps and missing
        values are defined by assuming a frequency per station.

        Parameters
        ----------
        df : pandas.dataframe
            The dataframe containing the input observations and metadata.
        freq_estimation_method : 'highest' or 'median'
            Select wich method to use for the frequency estimation. If
            'highest', the highest apearing frequency is used. If 'median', the
            median of the apearing frequencies is used.
        freq_estimation_simplify : bool
            If True, the likely frequency is converted to round hours, or round minutes.
            The "freq_estimation_simplify_error' is used as a constrain. If the constrain is not met,
            the simplification is not performed.
        freq_estimation_simplify_error : Timedelta or str, optional
            The tolerance string or object representing the maximum translation in time to form a simplified frequency estimation.
            Ex: '5min' is 5 minutes, '1h', is one hour.
        fixed_freq_series : pandas.series or None, optional
            If you do not want the frequencies to be recalculated, one can pass the
            frequency series to update the metadf["dataset_resolution"]. If None, the frequencies will be estimated. The default is None.
        update_full_metadf : bool, optional
            If True, the full Dataset.metadf will be updated. If False, only the frequency columns in the Dataset.metadf will be updated. The default is True.


        Returns
        -------
        None.

        """
        # Convert dataframe to dataset attributes
        self._initiate_df_attribute(dataframe=df, update_metadf=update_full_metadf)

        # Check observation types and convert units if needed.
        self._setup_of_obstypes_and_units()

        # Apply quality control on Import resolution
        self._apply_qc_on_import()

        if fixed_freq_series is None:
            freq_series = get_freqency_series(
                df=self.df,
                method=freq_estimation_method,
                simplify=freq_estimation_simplify,
                max_simplify_error=freq_estimation_simplify_error,
            )

            freq_series_import = freq_series

        else:
            if "assumed_import_frequency" in self.metadf.columns:
                freq_series_import = self.metadf[
                    "assumed_import_frequency"
                ]  # No update
            else:
                freq_series_import = fixed_freq_series
            freq_series = fixed_freq_series

        # add import frequencies to metadf (after import qc!)
        self.metadf["assumed_import_frequency"] = freq_series_import
        self.metadf["dataset_resolution"] = freq_series

        # Remove gaps and missing from the observations AFTER timecoarsening
        self.df = remove_gaps_from_obs(gaplist=self.gaps, obsdf=self.df)
        self.df = self.missing_obs.remove_missing_from_obs(obsdf=self.df)

    def _initiate_df_attribute(self, dataframe, update_metadf=True):
        """Initialize dataframe attributes."""
        logger.info(f"Updating dataset by dataframe with shape: {dataframe.shape}.")

        # Create dataframe with fixed order of observational columns
        obs_col_order = [
            col for col in list(self.obstypes.keys()) if col in dataframe.columns
        ]

        self.df = dataframe[obs_col_order].sort_index()

        if update_metadf:
            metadf = dataframe

            # create metadataframe
            metacolumns = list(self.template._get_metadata_column_map().values())
            metacolumns.remove("name")  # This is the index
            metadf = metadf.reindex(columns=metacolumns)
            metadf.index = metadf.index.droplevel("datetime")  # drop datetimeindex

            # drop dubplicates due to datetime
            metadf = metadf[~metadf.index.duplicated(keep="first")]

            # "lat' and 'lon' column are required, so add them as empty if not present
            if "lat" not in metadf.columns:
                metadf["lat"] = np.nan
            if "lon" not in metadf.columns:
                metadf["lon"] = np.nan

            self.metadf = metadf_to_gdf(metadf)

    def _apply_qc_on_import(self):
        # if the name is Nan, remove these records from df, and metadf (before)
        # they end up in the gaps and missing obs
        if np.nan in self.df.index.get_level_values("name"):
            logger.warning(
                f'Following observations are not linked to a station name and will be removed: {xs_save(self.df, np.nan, "name")}'
            )
            self.df = self.df[~self.df.index.get_level_values("name").isna()]
        if np.nan in self.metadf.index:
            logger.warning(
                f"Following station will be removed from the Dataset {self.metadf[self.metadf.index.isna()]}"
            )
            self.metadf = self.metadf[~self.metadf.index.isna()]

        # find missing obs and gaps, and remove them from the df
        self.missing_obs, self.gaps = missing_timestamp_and_gap_check(
            df=self.df,
            gapsize_n=self.settings.gap["gaps_settings"]["gaps_finder"]["gapsize_n"],
        )

        # Create gaps and missing obs objects
        # self.gaps = gaps_list
        # self.missing_obs = Missingob_collection(missing_obs)

        # Perform QC checks on original observation frequencies
        self.df, dup_outl_df = duplicate_timestamp_check(
            df=self.df,
            checks_info=self.settings.qc["qc_checks_info"],
            checks_settings=self.settings.qc["qc_check_settings"],
        )
        if not dup_outl_df.empty:
            self.update_outliersdf(add_to_outliersdf=dup_outl_df)

        self.df, nan_outl_df = invalid_input_check(
            self.df, checks_info=self.settings.qc["qc_checks_info"]
        )
        if not nan_outl_df.empty:
            self.update_outliersdf(nan_outl_df)

        self.outliersdf = self.outliersdf.sort_index()

        # update the order and which qc is applied on which obstype
        checked_obstypes = [
            obs for obs in self.df.columns if obs in self.obstypes.keys()
        ]

        checknames = ["duplicated_timestamp", "invalid_input"]  # KEEP order

        self._applied_qc = concat_save(
            [
                self._applied_qc,
                conv_applied_qc_to_df(
                    obstypes=checked_obstypes, ordered_checknames=checknames
                ),
            ],
            ignore_index=True,
        )

    def _setup_of_obstypes_and_units(self):
        """Function to setup all attributes related to observation types and
        convert to standard units."""

        # Check if all present observation types are known.
        unknown_obs_cols = [
            obs_col
            for obs_col in self.df.columns
            if obs_col not in self.obstypes.keys()
        ]
        if len(unknown_obs_cols) > 0:
            sys.exit(f"The following observation types are unknown: {unknown_obs_cols}")

        for obs_col in self.df.columns:
            # Convert the units to the toolkit standards (if unit is known)

            self.df[obs_col] = self.obstypes[obs_col].convert_to_standard_units(
                input_data=self.df[obs_col],
                input_unit=self.template._get_input_unit_of_tlk_obstype(obs_col),
            )

            # Update the description of the obstype
            description = self.template._get_description_of_tlk_obstype(obs_col)
            if pd.isna(description):
                description = None
            self.obstypes[obs_col].set_description(desc=description)

            # Update the original column name and original units
            self.obstypes[obs_col].set_original_name(
                self.template._get_original_obstype_columnname(obs_col)
            )
            self.obstypes[obs_col].set_original_unit(
                self.template._get_input_unit_of_tlk_obstype(obs_col)
            )

        # subset the obstypes attribute
        self.obstypes = {
            name: obj for name, obj in self.obstypes.items() if name in self.df.columns
        }

    # =============================================================================
    # Physiography extractions
    # =============================================================================
    def get_lcz(self):
        """Extract local climate zones for all stations.

        Function to extract the Local CLimate zones (LCZ) from the
        wudapt global LCZ map on the Google engine for all stations.

        A 'LCZ' column will be added to the metadf, and series is returned.

        Returns
        -------
        lcz_series : pandas.Series()
            A series with the stationnames as index and the LCZ as values.

        Examples
        --------
        .. code-block:: python

             import metobs_toolkit

             # Import data into a Dataset
             dataset = metobs_toolkit.Dataset()
             dataset.update_settings(
                                     input_data_file=metobs_toolkit.demo_datafile,
                                     input_metadata_file=metobs_toolkit.demo_metadatafile,
                                     template_file=metobs_toolkit.demo_template,
                                     )
             dataset.import_data_from_file()

             # Get the local climate zones for all stations
             lcz_series = dataset.get_lcz()

             # in addition to the returned series, the metadf attribute is updated aswell
             print(dataset.metadf)


        """
        # connect to gee
        connect_to_gee()

        # Extract LCZ for all stations
        lcz_series = lcz_extractor(
            metadf=self.metadf,
            mapinfo=self.settings.gee["gee_dataset_info"]["global_lcz_map"],
        )

        # drop column if it was already present
        if "lcz" in self.metadf:
            self.metadf = self.metadf.drop(columns=["lcz"])

        # update metadata
        self.metadf = self.metadf.merge(
            lcz_series.to_frame(),
            how="left",
            left_index=True,
            right_index=True,
        )
        return lcz_series

    def get_altitude(self):
        """Extract Altitudes for all stations.

        Function to extract the Altitude from the SRTM Digital Elevation Data
        global map on the Google engine for all stations.

        A 'altitude' column will be added to the metadf, and series is returned.

        Returns
        -------
        altitude_series : pandas.Series()
            A series with the stationnames as index and the altitudes as values.

         Examples
         --------
         .. code-block:: python

              import metobs_toolkit

              # Import data into a Dataset
              dataset = metobs_toolkit.Dataset()
              dataset.update_settings(
                                      input_data_file=metobs_toolkit.demo_datafile,
                                      input_metadata_file=metobs_toolkit.demo_metadatafile,
                                      template_file=metobs_toolkit.demo_template,
                                      )
              dataset.import_data_from_file()

              # Get the altitude for all stations
              alt_series = dataset.get_altitude()

              # in addition to the returned series, the metadf attribute is updated aswell
              print(dataset.metadf)

        """
        # connect to gee
        connect_to_gee()

        # Extract LCZ for all stations
        altitude_series = height_extractor(
            metadf=self.metadf,
            mapinfo=self.settings.gee["gee_dataset_info"]["DEM"],
        )

        # drop column if it was already present
        if "altitude" in self.metadf:
            self.metadf = self.metadf.drop(columns=["altitude"])

        # update metadata
        self.metadf = self.metadf.merge(
            altitude_series.to_frame(),
            how="left",
            left_index=True,
            right_index=True,
        )
        return altitude_series

    def get_landcover(
        self,
        buffers=[100],
        aggregate=True,
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
        buffers : num, optional
            The list of buffer radia in dataset units (meters for ESA worldcover) . The default is 100.
        aggregate : bool, optional
            If True, the classes will be aggregated with the corresponding
            aggregation scheme. The default is True.
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

        Examples
        --------
        .. code-block:: python

             import metobs_toolkit

             # Import data into a Dataset
             dataset = metobs_toolkit.Dataset()
             dataset.update_settings(
                                     input_data_file=metobs_toolkit.demo_datafile,
                                     input_metadata_file=metobs_toolkit.demo_metadatafile,
                                     template_file=metobs_toolkit.demo_template,
                                     )
             dataset.import_data_from_file()

             # Get the landcover fractions for multiple buffers, for all stations
             lc_frac_series = dataset.get_landcover(buffers=[50, 100, 250, 500],
                                                    aggregate=False)

             # in addition to the returned dataframe, the metadf attribute is updated aswell
             print(dataset.metadf)

        """
        # connect to gee
        connect_to_gee()

        df_list = []
        for buffer in buffers:
            logger.info(
                f"Extracting landcover from {gee_map} with buffer radius = {buffer}"
            )
            # Extract landcover fractions for all stations
            lc_frac_df, buffer = lc_fractions_extractor(
                metadf=self.metadf,
                mapinfo=self.settings.gee["gee_dataset_info"][gee_map],
                buffer=buffer,
                agg=aggregate,
            )

            # add buffer to the index
            lc_frac_df["buffer_radius"] = buffer
            lc_frac_df = lc_frac_df.reset_index().set_index(["name", "buffer_radius"])
            lc_frac_df = lc_frac_df.sort_index()

            # add to the list
            df_list.append(lc_frac_df)

        # concat all df for different buffers to one
        frac_df = concat_save(df_list)
        frac_df = frac_df.sort_index()

        if overwrite:
            for buf in frac_df.index.get_level_values("buffer_radius").unique():
                buf_df = xs_save(frac_df, buf, level="buffer_radius")
                buf_df.columns = [col + f"_{int(buf)}m" for col in buf_df.columns]

                # overwrite the columns or add them if they did not exist
                self.metadf[buf_df.columns] = buf_df

        return frac_df
