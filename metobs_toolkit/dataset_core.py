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
    # invalid_input_check,
    #     toolkit_buddy_check,
    #     titan_buddy_check,
    #     titan_sct_resistant_check,
)


# from metobs_toolkit.qc_statistics import get_freq_statistics
from metobs_toolkit.writing_files import (
    # write_dataset_to_csv,
    write_df_to_csv,
    _does_trg_file_exist,
    _remove_file,
    MetobsWritingError,
)

# from metobs_toolkit.missingobs import Missingob_collection

from metobs_toolkit.gap import (
    # Gap,
    find_gaps,
    # remove_gaps_from_obs,
    # remove_gaps_from_outliers,
    # missing_timestamp_and_gap_check,
    # get_gaps_indx_in_obs_space,
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
    empty_outliers_df,
    # init_triple_multiindexdf,
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
        # self.outliersdf = init_triple_multiindexdf()

        # self.missing_obs = None  # becomes a Missingob_collection after import
        # self.gaps = None  # becomes a list of gaps

        self.gapfilldf = init_multiindexdf()
        self.missing_fill_df = init_multiindexdf()

        # # Template for mapping data and metadata
        # self.template = Template()

        self._istype = "Dataset"
        self._freqs = pd.Series(dtype=object)

        self._applied_qc = pd.DataFrame(columns=["obstype", "checkname"])
        self._qc_checked_obstypes = []  # list with qc-checked obstypes

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
            sta_outliers = empty_outliers_df()

        sta_gaps = get_station_gaps(self.gaps, stationname)

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
            # gapfilldf=sta_gapfill,
            # missing_fill_df=sta_missingfill,
            metadf=sta_metadf,
            obstypes=self.obstypes,
            template=self.template,
            settings=self.settings,
            # _qc_checked_obstypes=self._qc_checked_obstypes,
            # _applied_qc=self._applied_qc,
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
            mergedf = self.get_full_status_df()

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
        data_file,
        metadata_file=None,
        overwrite=False,
        all_outliers_as_nan=False,
        kwargs_to_csv={"index": False},
    ):
        """Write Dataset to a csv file.

        Write all present records to a file. This is done by combining the
        good records, outliers and gaps into one dataframe. An extra
        'label-column', for each observationtype is added.

        The metadata can be written to a seperate file.


        Parameters
        ----------
        data_file : str or None
            The path to a csv location to write the data to. If the path does
            not have an ".csv" extension, it will be added. If None, no data is
            written.
        metadata_file : str or None, optional
            The path to a csv location to write the metadata to. If the path does
            not have an ".csv" extension, it will be added. If None, no metadata is
            written.
        overwrite : bool, optional
            If True, the target files will be written even if they already exist.
            The default is False.
        all_outliers_as_nan : bool, optional
            If True, all records flagged as outlier are represented by Nan.
            The default is False.
        kwargs_to_csv : dict, optional
            Kwargs that are passed to the pandas.DataFrame.to_csv() method.
            The default is {'index': False}.

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

        # check path to target files
        if data_file is None:
            # do not write data
            write_data = False
        else:
            write_data = True
            # make sure it has a .csv filetype extension
            if not str(data_file).endswith(".csv"):
                # add .csv
                logger.warning(
                    f'{data_file} has no ".csv" extension, this will be added to the path.'
                )
                data_file = f"{data_file}.csv"
            # check if path exist
            if _does_trg_file_exist(str(data_file)):
                if not overwrite:
                    raise MetobsWritingError(
                        f"{data_file} already exists. Use overwrite=True, or change the path."
                    )
                else:
                    _remove_file(str(data_file))

        # check path to target files
        if metadata_file is None:
            # do not write data
            write_metadata = False
        else:
            write_metadata = True
            if not str(metadata_file).endswith(".csv"):
                # add .csv
                logger.warning(
                    f'{metadata_file} has no ".csv" extension, this will be added to the path.'
                )
                metadata_file = f"{metadata_file}.csv"
            # check if path exist
            if _does_trg_file_exist(str(metadata_file)):
                if not overwrite:
                    raise MetobsWritingError(
                        f"{metadata_file} already exists. Use overwrite=True, or change the path."
                    )
                else:
                    _remove_file(str(metadata_file))

        if (not write_data) & (not write_metadata):
            raise MetobsDatasetCoreError(
                f"Cannot write data to a file since no target datafile or target metadatafile is given."
            )

        # =============================================================================
        # Write data
        # =============================================================================

        if write_data:
            combdf = self.get_full_status_df()

            present_obs = combdf.columns.get_level_values(0).unique()
            # Too a single columnindex
            combdf.columns = combdf.columns.map("_".join)

            # drop the toolkit representation columns
            combdf = combdf.drop(
                columns=[
                    col
                    for col in combdf.columns
                    if col.endswith("_toolkit_representation")
                ]
            )

            for ob in present_obs:
                if all_outliers_as_nan:
                    combdf.loc[combdf[f"{ob}_label"] != "ok", f"{ob}_value"] = np.nan
            # write to file
            write_df_to_csv(
                df=combdf,
                trgfile=str(data_file),
                to_csv_kwargs=kwargs_to_csv,
            )

        if write_metadata:
            metadf = self.metadf

            # add present obstype info to the metadf
            for obs in self.df.index.get_level_values("obstype").unique():
                obstype = self.obstypes[obs]
                metadf[f"{obs}_unit"] = obstype.get_standard_unit()
                metadf[f"{obs}_description"] = obstype.get_description()

            write_df_to_csv(
                df=metadf,
                trgfile=str(metadata_file),
                to_csv_kwargs=kwargs_to_csv,
            )
        return

    def coarsen_time_resolution(
        self,
        origin=None,
        origin_tz=None,
        freq="60min",
        direction="nearest",
        timestamp_mapping_tolerance="4min",
        limit=1,
    ):
        """Resample the observations to coarser timeresolution.

        This method is used to convert the time resolution of the Dataset. This
        will affect the records (.df), the outliers (.outliersdf) and gaps (.gaps).



        The assumed dataset resolution (stored in the metadf attribute) will be
        updated.

        Parameters
        ----------
        origin : datetime.datetime, optional
            Define the origin (first timestamp) for the obervations. The origin
            is timezone naive, and is assumed to have the same timezone as the
            obervations. If None, the earliest occuring timestamp is used as
            origin. The default is None.
        origin_tz : str or None, optional
            Timezone string of the input observations. Element of
            pytz.all_timezones. If None, the timezone of the records is used. The default is None.
        freq : DateOffset, Timedelta or str, optional
            The offset string or object representing target conversion.
            Ex: '15min' is 15 minutes, '1h', is one hour. The default is '60min'.
        direction : 'backward', 'forward', or 'nearest'
            Whether to search for prior, subsequent, or closest matches for
            mapping to ideal timestamps. The default is 'nearest'.
        timestamp_mapping_tolerance : Timedelta or str
            The tolerance string or object representing the maximum translation
            (in time) to map a timestamp to a target timestamp.
            Ex: '5min' is 5 minutes. The default is '4min'.
        limit : int, optional
            Limit of how many values to fill with one original observations.  The default
            is 1.

        Returns
        -------
        None.

        Warning
        ---------
        Since the gaps depend on the records frequency and origin, all gaps are
        removed and re-located. All progress in gap(filling) will be lost.

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

        if origin_tz is None:
            origin_tz = str(self._get_tz())

        #  ------ Define ideal timestamprecords ----------

        # 1. Frequency
        self.metadf.loc[:, "dataset_resolution"] = pd.Timedelta(freq)

        # 2. Origin
        if origin is not None:
            origin_tz_aware = pytz.timezone(origin_tz).localize(origin)
            tstart = origin_tz_aware.astimezone(
                pytz.timezone(self.settings.time_settings["timezone"])
            )

            # update the metadf
            self.metadf.loc[:, "dt_start"] = tstart

        # TODO: implement limit thing

        # df = self.df.reset_index()

        if origin is None:
            fixed_origin = False
        else:
            origin_tz_aware = pytz.timezone(origin_tz).localize(origin)
            tstart = origin_tz_aware.astimezone(
                pytz.timezone(self.settings.time_settings["timezone"])
            )
            fixed_origin = True

        # TODO: IDEE: voeg een label column toe, zo kunnen de outliers ook gecoarsened worden

        self.construct_equi_spaced_records(
            timestamp_mapping_tolerance=timestamp_mapping_tolerance, direction=direction
        )
        # stadf_list = []
        # # Coarsen timeresolution (per station, because of different origins)
        # for (staname,obsname), groupdf in self.df.reset_index().groupby(['name', 'obstype']):
        #     #Get the origin for this station
        #     if fixed_origin:
        #         sta_origin = tstart
        #     else:
        #         #use sta origne from metadf
        #         sta_origin = self.metadf.loc[staname, 'dt_start']

        #     #Resample
        #     if method == "nearest":
        #         stadf = (
        #             groupdf.set_index("datetime")
        #             .resample(freq, origin=sta_origin)
        #             .nearest(limit=limit)
        #         )

        #     elif method == "bfill":
        #         stadf = (
        #             groupdf.set_index("datetime")
        #             .resample(freq, origin=tstart)
        #             .bfill(limit=limit)
        #         )
        #     else:
        #         raise MetobsDatasetCoreError(f'{method} is not a valid method for coarsening the time resolution.')

        #     stadf_list.append(stadf)

        #     # #update metadf
        #     # self.metadf.loc[staname, "dataset_resolution"] = pd.to_timedelta(freq)
        #     # self.metadf.loc[staname, "dt_start"] = stadf.index.min()
        #     # self.metadf.loc[staname, "dt_end"] = stadf.index.max()

        # #combine all to one dataframe
        # df = pd.concat(stadf_list)

        # # A little cleanup
        # df = (df.reset_index()
        #       .set_index(['name','obstype','datetime'])
        #       .sort_index()
        #       )

        # self._set_df(df)

        # Update the metadf (The records are ideally at this point, so no simplifications)
        self._get_timestamps_info(
            freq_estimation_method="highest",
            freq_simplify_tolerance="0min",  # no simplification
            origin_simplify_tolerance="0min",  # no simplification
        )

        # # Convert the records to clean equidistanced records
        # self.construct_equi_spaced_records(
        #     timestamp_mapping_tolerance=timestamp_tolerance
        # )

        # # clear the outliers them again on the new resolution.
        # if not self.outliersdf.empty:
        #     logger.warning('Because of the resampling, all outliers are removed. Recomput them again!')
        #     self._set_outliersdf(empty_outliers_df())
        # # Remove duplicates (needed in order to convert the units and find gaps)
        # df, outliersdf = duplicate_timestamp_check(df=self.df,
        #                                           checks_info = self.settings.qc['qc_checks_info'],
        #                                           checks_settings=self.settings.qc['qc_check_settings']
        #                                           )
        # self._set_df(df=df)
        # self._set_outliersdf(outliersdf)

        # clear the gaps and compute them again on the new resolution.
        if bool(self.gaps):
            logger.warning("Because of the resampling, all gaps are recomputed again!")

        gaps = find_gaps(
            df=self.df,
            metadf=self.metadf,
            outliersdf=self.outliersdf,
            obstypes=self.obstypes,
        )
        self._set_gaps(gaps)

        # # Remove gaps and missing from the observatios
        # # most gaps and missing are already removed but when increasing timeres,
        # # some records should be removed as well.
        # self.df = remove_gaps_from_obs(gaplist=self.gaps, obsdf=self.df)
        # self.df = self.missing_obs.remove_missing_from_obs(obsdf=self.df)

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
        self.outliersdf = empty_outliers_df()
        # self.gapfilldf = init_multiindexdf()
        # self.missing_obs = None
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
        freq_estimation_method="highest",
        freq_estimation_simplify_tolerance="2min",
        origin_simplify_tolerance="5min",
        timestamp_tolerance="4min",
        kwargs_data_read={},
        kwargs_metadata_read={},
    ):
        """Read observations from a csv file.

        The paths (data, metadata and template) are stored in the settings if
        Dataset.update_settings() is called on this object. These paths can be
        updated by adding them as argument to this method.

        The input data (and metadata) are interpreted by using a template
        (json file).

        In order to locate gaps, an ideal set of timestamps is exptected. This
        set of timestamps is computed for each station seperatly by:
            * Assuming a constant frequency. This frequency is estimated by using
            a freq_estimation_method. If multiple observationtypes are present,
            the assumed frequency is the highest of estimated frequency among
            the differnt observationtypes. To simplify the estimated frequency a
            freq_estimation_simplify_error can be specified.
            * A start timestamp (origin) is found for each station. If multiple observationtypes are present,
            the start timestamp is the first timestamp among
            the different observationtypes. The start
            timestamp can be simplified by specifying a origin_simplify_tolerance.
            * The last timestamp is found for each station by taking the timestamp
            which is closest and smaller then the latest timestamp found of a station,
            and is an element of the ideal set of timestamps.

        Each present observation record is linked to a timestamp of this ideal set,
        by using a 'nearest' merge. If the timediffernce is smaller than the
        timestamp_tolerance, the ideal timestamp is used. Else, the timestamp
        will be interpreted as a (part of a) gap.


        The Dataset attributes are set and the following checks are executed:
                * Duplicate check
                * Invalid input check
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
            median of the apearing frequencies is used. The default is 'highest'.
        freq_estimation_simplify_tolerance : Timedelta or str, optional
            The tolerance string or object representing the maximum translation
            in time to form a simplified frequency estimation.
            Ex: '5min' is 5 minutes, '1H', is one hour. A zero-tolerance (thus no
            simplification) can be set by '0min'. The default is '2min' (2 minutes).
        origin_simplify_tolerance : Timedelta or str, optional
            The tolerance string or object representing the maximum translation
            in time to apply on the start timestamp to create a simplified timestamp.
            Ex: '5min' is 5 minutes, '1H', is one hour. A zero-tolerance (thus no
            simplification) can be set by '0min'. The default is '5min' (5 minutes).
        timestamp_tolerance : Timedelta or str, optional
            The tolerance string or object representing the maximum translation
            in time to apply on a timestamp for conversion to an ideal set of timestamps.
            Ex: '5min' is 5 minutes, '1H', is one hour. A zero-tolerance (thus no
            simplification) can be set by '0min'. The default is '4min' (4 minutes).
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
        # Update paths to the input files, if given.
        if input_data_file is not None:
            self.update_settings(input_data_file=input_data_file)
        if input_metadata_file is not None:
            self.update_settings(input_metadata_file=input_metadata_file)
        if template_file is not None:
            self.update_settings(template_file=template_file)

        logger.info(f'Importing data from file: {self.settings.IO["input_data_file"]}')

        assert self.settings.templatefile is not None, "No templatefile is specified."

        # Read template
        self.template.read_template_from_file(jsonpath=self.settings.templatefile)

        # Read observations into pandas dataframe
        df = import_data_from_csv(
            input_file=self.settings.IO["input_data_file"],
            template=self.template,
            known_obstypes=list(self.obstypes.keys()),
            kwargs_data_read=kwargs_data_read,
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

        self._construct_dataset(
            df=df,
            freq_estimation_method=freq_estimation_method,
            freq_estimation_simplify_tolerance=freq_estimation_simplify_tolerance,
            origin_simplify_tolerance=origin_simplify_tolerance,
            timestamp_tolerance=timestamp_tolerance,
        )

    def _construct_dataset(
        self,
        df,
        freq_estimation_method,
        freq_estimation_simplify_tolerance,
        origin_simplify_tolerance,
        timestamp_tolerance,
        # fixed_freq_series=None,
        # update_full_metadf=True,
    ):
        """Construct the Dataset class from a IO dataframe.

        1. Set the dataframe and metadataframe attributes
        2. Drop stations that have Nan as name.
        2. Find the duplicates (remove them from observations +  add them to outliers)
        3. Convert the values to standard units + update the observationtypes (some template specific attribute)
        5. Find gaps in the records (duplicates are excluded from the gaps)
        6. Get a frequency estimate per station
        7. Initiate the gaps (find missing records)
        8. Add the missing records to the dataframe


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
            Ex: '5min' is 5 minutes, '1H', is one hour.
        fixed_freq_series : pandas.series or None, optional
            If you do not want the frequencies to be recalculated, one can pass the
            frequency series to update the metadf["dataset_resolution"]. If None, the frequencies will be estimated. The default is None.
        update_full_metadf : bool, optional
            If True, the full Dataset.metadf will be updated. If False, only the frequency columns in the Dataset.metadf will be updated. The default is True.


        Returns
        -------
        None.

        """
        # Set the df attribute
        self._construct_df(dataframe=df)
        # Set the metadf attribute
        self._construct_metadf(dataframe=df)

        # Apply QC on Nan and duplicates (needed before unit conversion and gapcreation)
        # Remove nan names
        self._remove_nan_names()

        # Convert to numeric --> "invalid check' will be triggered if not possible
        self._to_num_and_invalid_check()

        # Remove duplicates (needed in order to convert the units and find gaps)
        df, outliersdf = duplicate_timestamp_check(
            df=self.df,
            outlierlabel=self.settings.label_def["duplicated_timestamp"]["label"],
            checks_settings=self.settings.qc["qc_check_settings"],
        )
        self._set_df(df=df)
        self._update_outliersdf(outliersdf)

        # self._covert_timestamps_to_utc()

        # Check observation types and convert units if needed.
        self._setup_of_obstypes_and_units()

        # find the start, end timestamps and frequency for each station + write it to the metadf
        self._get_timestamps_info(
            freq_estimation_method=freq_estimation_method,
            freq_simplify_tolerance=freq_estimation_simplify_tolerance,
            origin_simplify_tolerance=origin_simplify_tolerance,
        )

        # Convert the records to clean equidistanced records for both the df and outliersdf
        self.construct_equi_spaced_records(
            timestamp_mapping_tolerance=timestamp_tolerance
        )

        # # Find gaps on Import resolution
        gaps = find_gaps(
            df=self.df,
            metadf=self.metadf,
            outliersdf=self.outliersdf,
            obstypes=self.obstypes,
        )
        self._set_gaps(gaps)

    def _to_num_and_invalid_check(self):
        # 8. map to numeric dtypes
        # When converting to numeric, this overrules the invalid check.

        df = self.df

        # 1 subset to the records with Nan values --> do not check these, just add them back in the end
        nandf = df[~df["value"].notnull()]
        # 2 Get the other subet with values not nan (can be numerics and strings) --> filter out the strings
        to_checkdf = df[df["value"].notnull()]

        # 3 Convert to numeric
        to_checkdf["value"] = pd.to_numeric(to_checkdf["value"], errors="coerce")

        # 4 All the Nan's in the to_checkdf are outliers triggerd as 'invalid'
        invalid_records = to_checkdf[~to_checkdf["value"].notnull()]
        # add the label of "invalid check' to it
        invalid_records["label"] = self.settings.label_def["invalid_input"]["label"]
        # special case: duplicates in the invalid records
        invalid_records = invalid_records[~invalid_records.index.duplicated()]

        # 5. Combine the df's back to one
        totaldf = pd.concat([nandf, to_checkdf]).sort_index()

        # Set attributes
        # Note that at this point, the duplicated check is not performed yet.
        self._set_df(totaldf, apply_dup_checks=False)
        self._update_outliersdf(invalid_records)
        return

    # def _initiate_df_attribute(self, dataframe, update_metadf=True):
    #     """Initialize dataframe attributes."""
    #     logger.info(f"Updating dataset by dataframe with shape: {dataframe.shape}.")

    #     # Create dataframe with fixed order of observational columns
    #     obs_col_order = [
    #         col for col in list(self.obstypes.keys()) if col in dataframe.columns
    #     ]

    #     self.df = dataframe[obs_col_order].sort_index()

    #     if update_metadf:
    #         metadf = dataframe

    #         # create metadataframe
    #         metacolumns = list(self.template._get_metadata_column_map().values())
    #         metacolumns.remove("name")  # This is the index
    #         metadf = metadf.reindex(columns=metacolumns)
    #         metadf.index = metadf.index.droplevel("datetime")  # drop datetimeindex

    #         # drop dubplicates due to datetime
    #         metadf = metadf[~metadf.index.duplicated(keep="first")]

    #         # "lat' and 'lon' column are required, so add them as empty if not present
    #         if "lat" not in metadf.columns:
    #             metadf["lat"] = np.nan
    #         if "lon" not in metadf.columns:
    #             metadf["lon"] = np.nan

    #         self.metadf = metadf_to_gdf(metadf)

    # def _apply_qc_on_import(self):
    #     # if the name is Nan, remove these records from df, and metadf (before)
    #     # they end up in the gaps and missing obs
    #     if np.nan in self.df.index.get_level_values("name"):
    #         logger.warning(
    #             f'Following observations are not linked to a station name and will be removed: {xs_save(self.df, np.nan, "name")}'
    #         )
    #         self.df = self.df[~self.df.index.get_level_values("name").isna()]
    #     if np.nan in self.metadf.index:
    #         logger.warning(
    #             f"Following station will be removed from the Dataset {self.metadf[self.metadf.index.isna()]}"
    #         )
    #         self.metadf = self.metadf[~self.metadf.index.isna()]

    #     # find missing obs and gaps, and remove them from the df
    #     self.missing_obs, self.gaps = missing_timestamp_and_gap_check(
    #         df=self.df,
    #         gapsize_n=self.settings.gap["gaps_settings"]["gaps_finder"]["gapsize_n"],
    #     )

    #     # Create gaps and missing obs objects
    #     # self.gaps = gaps_list
    #     # self.missing_obs = Missingob_collection(missing_obs)

    #     # Perform QC checks on original observation frequencies
    #     self.df, dup_outl_df = duplicate_timestamp_check(
    #         df=self.df,
    #         checks_info=self.settings.qc["qc_checks_info"],
    #         checks_settings=self.settings.qc["qc_check_settings"],
    #     )
    #     if not dup_outl_df.empty:
    #         self.update_outliersdf(add_to_outliersdf=dup_outl_df)

    #     self.df, nan_outl_df = invalid_input_check(
    #         self.df, checks_info=self.settings.qc["qc_checks_info"]
    #     )
    #     if not nan_outl_df.empty:
    #         self.update_outliersdf(nan_outl_df)

    #     self.outliersdf = self.outliersdf.sort_index()

    #     # update the order and which qc is applied on which obstype
    #     checked_obstypes = [
    #         obs for obs in self.df.columns if obs in self.obstypes.keys()
    #     ]

    #     checknames = ["duplicated_timestamp", "invalid_input"]  # KEEP order

    #     self._applied_qc = concat_save(
    #         [
    #             self._applied_qc,
    #             conv_applied_qc_to_df(
    #                 obstypes=checked_obstypes, ordered_checknames=checknames
    #             ),
    #         ],
    #         ignore_index=True,
    #     )

    def _setup_of_obstypes_and_units(self):
        """Function to setup all attributes related to observation types and
        convert to standard units."""
        # Check if all present observation types are known.
        present_obstypes = self._get_present_obstypes()

        # check if all present obstypes (in the df), are linked to a knonw Obstype
        self._are_all_present_obstypes_knonw()

        # Found that it is approx 70 times faster convert obstype per obstype,
        # add them to a list and concat them, than using the .loc method to
        # assign the converted values directly to the df attribute

        subdf_list = []
        for present_obs in present_obstypes:
            # Convert the units to the toolkit standards (if unit is known)
            input_unit = self.template._get_input_unit_of_tlk_obstype(present_obs)

            # locate the specific obstype records
            obstype_values = xs_save(self.df, present_obs, "obstype", drop_level=False)[
                "value"
            ]
            # Convert to standard unit and replace them in the df attribute
            subdf_list.append(
                self.obstypes[present_obs].convert_to_standard_units(
                    input_data=obstype_values, input_unit=input_unit
                )
            )
            # Update the description of the obstype
            self.obstypes[present_obs].set_description(
                desc=self.template._get_description_of_tlk_obstype(present_obs)
            )

            # Update the original name of the obstype (used for titles in plots)
            self.obstypes[present_obs].set_original_name(
                columnname=self.template._get_original_obstype_columnname(present_obs)
            )

            # Update the original unit of the obstype (not an application yet)
            self.obstypes[present_obs].set_original_unit(input_unit)

        df = pd.concat(subdf_list).to_frame().sort_index()
        self._set_df(df)

        # # subset the obstypes attribute
        # self.obstypes = {
        #     name: obj for name, obj in self.obstypes.items() if name in present_obstypes
        # }

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

    def _construct_df(self, dataframe):
        """fill the df attribute

        The dataframe is wide with data and metadata combined. This method will
        subset and format it to a long structured df to be set as the df attribute.
        """

        # subset to name, datetime and all obstypes
        df_cols = self.template._get_all_mapped_data_cols_in_tlk_space()
        # name and datetime are already index of dataframe, so drop them from df_cols
        df_cols = list(set(df_cols) - set(["name", "datetime"]))
        # subset the column
        df = dataframe.loc[:, df_cols]

        # convert the wide df to a long format
        triple_df = df.stack(future_stack=True)
        # rename the last level of the index to obstype
        triple_df.index.rename(names="obstype", level=-1, inplace=True)
        # rename the series to value
        triple_df.rename("value", inplace=True)
        # fix index order
        triple_df = triple_df.reorder_levels(
            ["name", "obstype", "datetime"]
        ).sort_index()
        # sort by index (by name --> datetime --> obstype)
        triple_df.sort_index(inplace=True)
        # convert to frame
        # TODO is this needed?
        triple_df = triple_df.to_frame()

        # set the attribute
        self._set_df(
            df=triple_df, apply_dup_checks=False
        )  # duplicate have yet to be removed

    def _construct_metadf(self, dataframe):
        """fill the metadf attribute

        The dataframe is wide with data and metadata combined. This method will
        subset and format the data and set the metadf attribute.


        """

        meta_cols = list(self.template._get_metadata_column_map().values())

        metadf = (
            dataframe.reset_index()
            .loc[:, meta_cols]
            .drop_duplicates()
            .set_index("name")
        )

        # Construct columns that are required
        if "lat" not in metadf.columns:
            metadf["lat"] = np.nan
        if "lon" not in metadf.columns:
            metadf["lon"] = np.nan

        # Convert to geopandas dataframe
        metadf = metadf_to_gdf(metadf)

        # set the attribute
        self._set_metadf(metadf=metadf)


def _construct_obsdf_and_metadf(dataframe, knownobstypes):

    # Create dataframe with fixed order of observational columns
    obs_col_order = [
        col for col in list(knownobstypes.keys()) if col in dataframe.columns
    ]

    df = dataframe[obs_col_order].sort_index()
    # convert the wide df to a long format
    triple_df = (
        # convert columns to an index level and reset the index
        df.stack().reset_index()
        # rename the default genereted columns
        .rename(columns={"level_2": "obstype", 0: "value"})
        # set and order the triple index
        .set_index(["name", "obstype", "datetime"])
        # sort by index
        .sort_index()
    )

    metacols = [col for col in dataframe.columns if col not in obs_col_order]
    metadf = dataframe[metacols]
    # create metadataframe with fixed number and order of columns
    # metadf = dataframe.reindex(columns=self.settings.app["location_info"])
    metadf.index = metadf.index.droplevel("datetime")  # drop datetimeindex
    # drop dubplicates due to datetime
    metadf = metadf[~metadf.index.duplicated(keep="first")]

    # Construct columns that are required
    if "lat" not in metadf.columns:
        metadf["lat"] = np.nan
    if "lon" not in metadf.columns:
        metadf["lon"] = np.nan

    # Convert to geopandas dataframe
    metadf = metadf_to_gdf(metadf)

    return triple_df, metadf


# =============================================================================
# Exceptions
# =============================================================================


class MetobsDatasetCoreError(Exception):
    """Exception raised for errors in the template."""

    pass
