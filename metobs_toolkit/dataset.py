import os
import copy
import pickle
import logging
import warnings
from typing import Literal, Union
from pathlib import Path

import pandas as pd
import numpy as np
from matplotlib.pyplot import Axes
import concurrent.futures

from metobs_toolkit.backend_collection.df_helpers import save_concat
from metobs_toolkit.template import (
    Template,
    update_known_obstype_with_original_data
    )
from metobs_toolkit.station import Station
from metobs_toolkit.io_collection.metadataparser import MetaDataParser
from metobs_toolkit.io_collection.dataparser import DataParser
from metobs_toolkit.io_collection.filereaders import (
    CsvFileReader, PickleFileReader
    )
from metobs_toolkit.site import Site
from metobs_toolkit.sensordata import SensorData
from metobs_toolkit.backend_collection.argumentcheckers import (
    fmt_timedelta_arg,
    fmt_datetime_arg,
)
from metobs_toolkit.timestampmatcher import simplify_time
from metobs_toolkit.obstypes import tlk_obstypes
from metobs_toolkit.obstypes import Obstype

import metobs_toolkit.plot_collection as plotting
import metobs_toolkit.backend_collection.printing_collection as printing

from metobs_toolkit.qc_collection import toolkit_buddy_check
from metobs_toolkit.backend_collection.dev_collection import copy_doc
from metobs_toolkit.backend_collection.errorclasses import (
    MetObsStationNotFound,
    MetObsDataAlreadyPresent,
    MetObsMissingFile,
    MetObsObstypeNotFound,
    MetObsMissingArgument,
    MetObsMetadataNotFound
)
from metobs_toolkit.modeltimeseries import ModelTimeSeries
from metobs_toolkit.settings_collection import label_def

from metobs_toolkit.geedatasetmanagers import (
    GEEStaticDatasetManager,
    GEEDynamicDatasetManager,
    default_datasets,
)
from metobs_toolkit.gee_api import connect_to_gee

logger = logging.getLogger("<metobs_toolkit>")


class Dataset:
    """
    Dataset class for managing and processing meteorological observation data.

    This class provides functionality to handle datasets containing observations from multiple stations.
    It includes methods for data synchronization, quality control, gap filling, and integration with
    Google Earth Engine (GEE) datasets. The Dataset class also supports importing and exporting data,
    as well as generating plots for visualization.

    Attributes
    ----------
    stations : list
        A list of Station objects representing the stations in the dataset.
    obstypes : dict
        A dictionary of observation types known to the dataset.
    template : Template
        The template instance used for data import and processing.
    """

    def __init__(self) -> None:
        """
        Initialize a Dataset instance.

        Returns
        -------
        None
        """
        self._stations = []  # stationname: Station
        self._obstypes = copy.copy(tlk_obstypes)  # init with all tlk obstypes
        self._template = Template()
        logger.debug("Dataset instance created.")

    # ------------------------------------------
    #    specials
    # ------------------------------------------

    def __eq__(self, other: object) -> bool:
        """Check equality with another Dataset instance."""
        logger.debug("Entering Dataset.__eq__")
        if not isinstance(other, Dataset):
            return False
        return self.stations == other.stations

    def __str__(self) -> str:
        """Return a string representation of the Dataset."""
        logger.debug("Entering Dataset.__str__")
        return "Dataset instance"

    def __repr__(self) -> str:
        """Return an info representation of the Dataset."""
        logger.debug("Entering Dataset.__repr__")
        class_name = type(self).__name__
        return f"Instance of {class_name} at {hex(id(self))}"

    def copy(self, deep: bool = True) -> "Dataset":
        """
        Return a copy of the Dataset.

        Parameters
        ----------
        deep : bool, optional
            If True, perform a deep copy. Default is True.

        Returns
        -------
        Dataset
            The copied Dataset.
        """
        logger.debug("Entering Dataset.copy")

        if deep:
            return copy.deepcopy(self)
        return copy.copy(self)

    @property
    def stations(self) -> list:
        """
        Get the list of Stations present in the Dataset.
        """
        return self._stations

    @property
    def obstypes(self) -> dict:
        """
        Get the dictionary of known Obstypes by the Dataset.

        """
        return self._obstypes

    @property
    def template(self) -> Template:
        """
        Get the Template instance used when the data was imported.

        """

        return self._template

    @stations.setter
    def stations(self, stationlist: list) -> None:
        """
        Set the list of stations.
        """

        self._stations = stationlist

    @obstypes.setter
    def obstypes(self, obstypesdict: dict) -> None:
        """
        Set the obstypes.

        """

        self._obstypes = obstypesdict

    @copy_doc(Station.df)
    @property
    def df(self) -> pd.DataFrame:
        concatlist = []
        for sta in self.stations:
            stadf = sta.df.reset_index()
            if stadf.empty:
                continue
            stadf["name"] = sta.name
            concatlist.append(stadf.set_index(["datetime", "obstype", "name"]))

        combdf = save_concat((concatlist))
        combdf.sort_index(inplace=True)
        if combdf.empty:
            combdf = pd.DataFrame(
                columns=["value", "label"],
                index=pd.MultiIndex(
                    levels=[[], [], []],
                    codes=[[], [], []],
                    names=["datetime", "obstype", "name"],
                ),
            )
        return combdf

    @copy_doc(Station.outliersdf)
    @property
    def outliersdf(self) -> pd.DataFrame:
        concatlist = []
        for sta in self.stations:
            stadf = sta.outliersdf.reset_index()
            stadf["name"] = sta.name
            concatlist.append(stadf.set_index(["datetime", "obstype", "name"]))

        combdf = save_concat((concatlist))
        combdf.sort_index(inplace=True)
        if combdf.empty:
            combdf = pd.DataFrame(
                columns=["value", "label"],
                index=pd.MultiIndex(
                    levels=[[], [], []],
                    codes=[[], [], []],
                    names=["datetime", "obstype", "name"],
                ),
            )
        return combdf

    @copy_doc(Station.gapsdf)
    @property
    def gapsdf(self) -> pd.DataFrame:
        concatlist = []
        for sta in self.stations:
            stadf = sta.gapsdf.reset_index()
            if stadf.empty:
                continue
            stadf["name"] = sta.name
            concatlist.append(stadf.set_index(["datetime", "obstype", "name"]))

        combdf = save_concat((concatlist))
        combdf.sort_index(inplace=True)
        if combdf.empty:
            combdf = pd.DataFrame(
                columns=["value", "label", "details"],
                index=pd.MultiIndex(
                    levels=[[], [], []],
                    codes=[[], [], []],
                    names=["datetime", "obstype", "name"],
                ),
            )
        return combdf

    @copy_doc(Station.modeldatadf)
    @property
    def modeldatadf(self) -> pd.DataFrame:
        concatlist = []
        for sta in self.stations:
            stadf = sta.modeldatadf.reset_index()
            if stadf.empty:
                continue
            stadf["name"] = sta.name
            concatlist.append(stadf.set_index(["datetime", "obstype", "name"]))

        combdf = save_concat((concatlist))
        combdf.sort_index(inplace=True)
        if combdf.empty:
            combdf = pd.DataFrame(
                columns=["value", "details"],
                index=pd.MultiIndex(
                    levels=[[], [], []],
                    codes=[[], [], []],
                    names=["datetime", "obstype", "name"],
                ),
            )
        return combdf

    @copy_doc(Station.metadf)
    @property
    def metadf(self) -> pd.DataFrame:
        concatlist = []
        for sta in self.stations:
            concatlist.append(sta.metadf)
        return save_concat((concatlist)).sort_index()

    @copy_doc(Station.start_datetime)
    @property
    def start_datetime(self) -> pd.Timestamp:
        return min([sta.start_datetime for sta in self.stations])

    @copy_doc(Station.end_datetime)
    @property
    def end_datetime(self) -> pd.Timestamp:
        return max([sta.end_datetime for sta in self.stations])

    @property
    def present_observations(self) -> list:
        """
        Get a list of all the present observation types.

        Returns
        -------
        list
            A list of all the present observations in the dataset.
        """
        allobs = set()
        for sta in self.stations:
            allobs.update(sta.present_observations)
        return sorted(list(allobs))

    # ------------------------------------------
    #   Extracting data
    # ------------------------------------------
    def subset_by_stations(
        self, stationnames: list, deepcopy: bool = False
    ) -> "Dataset":
        """
        Create a subset of the dataset by selecting specific stations.

        Parameters
        ----------
        stationnames : list
            A list of station names to filter the dataset by.
        deepcopy : bool, optional
            If True, creates a deep copy of the selected stations. Default is False.

        Returns
        -------
        Dataset
            A new Dataset instance containing only the selected stations.
        """
        logger.debug("Entering Dataset.subset_by_stations")

        new_dataset = Dataset()

        filtered_stations = [
            copy.deepcopy(sta) if deepcopy else sta
            for sta in self.stations
            if sta.name in stationnames
        ]
        if not filtered_stations:
            warnings.warn(
                "No stations matched the provided station names.", UserWarning
            )

        if len(filtered_stations) == 1:
            raise ValueError(
                "Only one station matches the provided names. Use the get_station() method instead."
            )

        new_dataset._stations = filtered_stations
        new_dataset._obstypes = self._obstypes
        new_dataset._template = self._template

        return new_dataset

    def get_station(self, stationname: str) -> Station:
        """
        Retrieve a Station by name.

        Parameters
        ----------
        stationname : str
            The name of the station to retrieve.

        Returns
        -------
        Station
            The station object corresponding to the given name.
        """
        logger.debug("Entering Dataset.get_station")

        stationlookup = {sta.name: sta for sta in self.stations}
        try:
            station = stationlookup[stationname]
        except KeyError:
            raise MetObsStationNotFound(f"{stationname} is not found in {self}.")

        return station

    @copy_doc(Station.get_info)
    def get_info(self, printout: bool = True) -> Union[str, None]:
        logger.debug("Entering Dataset.get_info")

        infostr = ""
        infostr += printing.print_fmt_title("General info of Dataset")

        # --- Observational info ---
        infostr += printing.print_fmt_section("Observational info")
        df = self.df
        if df.empty:
            infostr += printing.print_fmt_line(
                "Dataset instance without observation records."
            )
        else:
            present_obstypes = list(df.index.get_level_values("obstype").unique())

            infostr += printing.print_fmt_line("Dataset instance with:", 0)
            infostr += printing.print_fmt_line(
                f"{len(self.stations)} number of stations"
            )
            infostr += printing.print_fmt_line(f"{len(df.index)} number of records")
            infostr += printing.print_fmt_line(
                f"{len(present_obstypes)} types of sensor data are present."
            )
            infostr += printing.print_fmt_line(
                f"Observations from {self.start_datetime} -> {self.end_datetime}"
            )

            # -- outlier info --
            outldf = self.outliersdf
            infostr += printing.print_fmt_line("Outlier info:")
            if outldf.empty:
                infostr += printing.print_fmt_line(f"No QC outliers present.", 2)
            else:
                infostr += printing.print_fmt_line(
                    "A total of {outldf.shape[0]} outliers are present.", 2
                )
                infostr += printing.print_fmt_line(f"label counts:", 3)
                infostr += printing.print_fmt_dict(
                    outldf["label"].value_counts().to_dict(), identlvl=4
                )
                infostr += printing.print_fmt_line(
                    f"For these obstyes: {list(outldf.index.get_level_values('obstype').unique())}",
                    2,
                )
                unique_stations = list(outldf.index.get_level_values("name").unique())
                infostr += printing.print_fmt_line(
                    f"For {len(unique_stations)} stations: {unique_stations}", 2
                )

            # -- gap info --
            gapsdf = self.gapsdf
            infostr += printing.print_fmt_line("Gaps info:")
            if gapsdf.empty:
                infostr += printing.print_fmt_line(f"No gaps present.", 2)
            else:
                infostr += printing.print_fmt_line(
                    f"A total of {gapsdf.shape[0]} gaps are present.", 2
                )
                infostr += printing.print_fmt_line(f"label counts: ", 3)
                infostr += printing.print_fmt_dict(
                    gapsdf["label"].value_counts().to_dict(), identlvl=4
                )
                infostr += printing.print_fmt_line(
                    f"For these obstyes: {list(gapsdf.index.get_level_values('obstype').unique())}",
                    2,
                )
                unique_stations = list(gapsdf.index.get_level_values("name").unique())
                infostr += printing.print_fmt_line(
                    f"For {len(unique_stations)} stations: {unique_stations}", 2
                )

        # Meta data info
        infostr += printing.print_fmt_section("Metadata info")

        metadf = self.metadf
        if metadf.empty:
            infostr += printing.print_fmt_line("Dataset instance without metadata.")
        else:
            infostr += printing.print_fmt_line(
                f"{len(metadf.index)} number of stations"
            )
            infostr += printing.print_fmt_line(
                f"The following metadata is present: {list(metadf.columns)}"
            )

        # Modeldata info
        modeldf = self.modeldatadf
        infostr += printing.print_fmt_section("Modeldata info")
        if modeldf.empty:
            infostr += printing.print_fmt_line("Dataset instance without modeldata.")
        else:
            infostr += printing.print_fmt_line(
                f"Modeldata is present for {list(modeldf.index.get_level_values('obstype').unique())}"
            )
            infostr += printing.print_fmt_line(
                f"For period {modeldf.index.get_level_values('datetime').min()} -> {modeldf.index.get_level_values('datetime').max()}"
            )

        if printout:
            print(infostr)
        else:
            return infostr

    def sync_records(
        self,
        target_obstype: str = "temp",
        timestamp_shift_tolerance: Union[str, pd.Timedelta] = "2min",
        freq_shift_tolerance: Union[str, pd.Timedelta] = "1min",
        fixed_origin: Union[pd.Timedelta, None] = None,
    ) -> None:
        """Synchronize records of sensor data across stations.

                Synchronize records of sensor data across stations (for a specific observation type).
                This method aligns the sensor data of a specified observation type (`target_obstype`)
                across all stations by resampling the data to a common frequency and ensuring
                alignment errors of timestamps within specified tolerances.

                Parameters
                ----------
                target_obstype : str, optional
                    The observation type to synchronize (e.g., "temp" for temperature). Default is "temp".
                timestamp_shift_tolerance : str or pandas.Timedelta, optional
                    The maximum allowed time shift tolerance for aligning data during
                    resampling. Default is 2 minutes.
                freq_shift_tolerance : str or pandas.Timedelta, optional
                    The maximum allowed error in simplifying the target frequency. Default is "1min".
                fixed_origin : pandas.Timestamp, or None, optional
                    A fixed origin timestamp for resampling. If None, the origin is
                    determined automatically. Default is None.

                Returns
        -        -------
                None.

                Note
                ------
                In general, this method is a wrapper for `Dataset.resample()` but making sure that
                the target frequencies are naturel multiples of each other thus ensuring syncronisation
                accros stations.

                Warning
                ----------
                Since the gaps depend on the record’s frequency and origin, all gaps
                are removed and re-located. All progress in gap(filling) will be lost.

                Warning
                --------
                Cumulative tolerance errors can be introduced when this method is called multiple times.


        """

        # format arguments
        fixed_origin = fmt_timedelta_arg(fixed_origin)
        freq_shift_tolerance = fmt_timedelta_arg(freq_shift_tolerance)
        timestamp_shift_tolerance = fmt_timedelta_arg(timestamp_shift_tolerance)

        for sta in self.stations:
            if target_obstype in sta.sensordata.keys():
                sensor = sta.get_sensor(target_obstype)

                freq_target = simplify_time(
                    time=sensor.freq,
                    max_simplify_error=freq_shift_tolerance,
                    zero_protection=True,
                )

                # Note that simplify_time (tries to) simplifies to a target,
                # all targets are natural multipicates of each other --> ensuring syncing
                # over different sensors.

                sensor.resample(
                    target_freq=freq_target,
                    shift_tolerance=timestamp_shift_tolerance,
                    origin=fixed_origin,
                    origin_simplify_tolerance=timestamp_shift_tolerance,
                )
            else:
                logger.warning(
                    f"{sta} does not have {target_obstype} sensordata and is skipped in the synchronization."
                )

    @copy_doc(Station.resample)
    def resample(
        self,
        target_freq: Union[str, pd.Timedelta],
        target_obstype: Union[str, None] = None,
        shift_tolerance: Union[str, pd.Timedelta] = pd.Timedelta("4min"),
        origin: Union[str, pd.Timestamp, None] = None,
        origin_simplify_tolerance: Union[str, pd.Timedelta] = pd.Timedelta("4min"),
    ) -> None:
        logger.debug("Entering Dataset.resample")
        target_freq = fmt_timedelta_arg(target_freq)
        shift_tolerance = fmt_timedelta_arg(shift_tolerance)

        for sta in self.stations:
            sta.resample(
                target_freq=target_freq,
                target_obstype=target_obstype,
                shift_tolerance=shift_tolerance,
                origin=origin,
                origin_simplify_tolerance=origin_simplify_tolerance,
            )

    # ------------------------------------------
    #   Reading/writing data
    # ------------------------------------------

    def import_gee_data_from_file(
        self,
        filepath: str,
        geedynamicdatasetmanager: GEEDynamicDatasetManager,
        force_update: bool = True,
        _force_from_dataframe: Union[pd.DataFrame, None] = None,
    ) -> pd.DataFrame:
        """Import Google Earth Engine (GEE) data from a CSV file

        Import Google Earth Engine (GEE) data from a file, that was writen in your Google Drive,
        and integrate it into the dataset. Start by downloading the target CSV file from your
        Google Drive, then specify the path to the file in the `filepath` parameter.


        Parameters
        ----------
        filepath : str
            The path to the file containing the GEE data.
        geedynamicdatasetmanager : GEEDynamicDatasetManager
            An instance of `GEEDynamicDatasetManager` responsible for managing the GEE dataset
            and its metadata.
        force_update : bool, optional
            If True, forces the update of model data for the stations, by default True.
        _force_from_dataframe : pd.DataFrame, optional
            A pre-processed DataFrame to be used instead of reading from the file, by default None.

        Returns
        -------
        pd.DataFrame
            The processed DataFrame containing the GEE data.
        """
        logger.debug("Entering Dataset.import_gee_data_from_file")

        if _force_from_dataframe is None:
            reader = CsvFileReader(file_path=filepath)
            data = reader.read_as_local_file()
            force_update = True

            totaldf = geedynamicdatasetmanager._format_gee_df_structure(data)
            totaldf = geedynamicdatasetmanager._convert_units(totaldf)
        else:
            totaldf = _force_from_dataframe

        known_obstypes = list(geedynamicdatasetmanager.modelobstypes.keys())
        cols_to_skip = list(set(totaldf.columns) - set(known_obstypes))
        if bool(cols_to_skip):
            logger.warning(
                f"The following columns in the GEE datafile are not present in the known modelobstypes of {geedynamicdatasetmanager}: {cols_to_skip}"
            )

        known_and_present = set(totaldf.columns) & set(known_obstypes)
        df = totaldf[list(known_and_present)]

        for sta in self.stations:
            if sta.name not in df.index.get_level_values("name"):
                continue
            stadf = df.xs(key=sta.name, level="name", drop_level=True)

            for col in stadf.columns:
                modeltimeseries = ModelTimeSeries(
                    site=sta.site,
                    datarecords=stadf[col].to_numpy(),
                    timestamps=stadf.index.to_numpy(),
                    obstype=geedynamicdatasetmanager.modelobstypes[col],
                    timezone="UTC",
                    modelname=geedynamicdatasetmanager.name,
                    modelvariable=geedynamicdatasetmanager.modelobstypes[
                        col
                    ].model_band,
                )
                sta.add_to_modeldata(modeltimeseries, force_update=force_update)

        return totaldf

    def add_new_observationtype(self, obstype: Obstype) -> None:
        """
        Add a new observation type to the dataset known-obstypes.

        Parameters
        ----------
        obstype : Obstype
            An instance of the `Obstype` class representing the observation type
            to be added.

        Returns
        -------
        None
        """
        logger.debug("Entering Dataset.add_new_observationtype")
        if not isinstance(obstype, Obstype):
            raise TypeError("obstype must be an instance of Obstype.")

        if obstype.name in self.obstypes.keys():
            raise MetObsDataAlreadyPresent(
                f"An Obstype with {obstype.name} as name is already present in the obstypes: {self.obstypes}"
            )

        self.obstypes.update({obstype.name: obstype})

    def save_dataset_to_pkl(
        self,
        target_folder: Union[str, Path],
        filename: str = "saved_dataset.pkl",
        overwrite: bool = False,
    ) -> None:
        """
        Save the dataset to a pickle (.pkl) file.

        Parameters
        ----------
        target_folder : str or Path
            The folder where the pickle file will be saved. Must be an existing directory.
        filename : str, optional
            The name of the pickle file. Defaults to "saved_dataset.pkl". If the provided
            filename does not have a ".pkl" extension, it will be automatically appended.
        overwrite : bool, optional
            Whether to overwrite the file if it already exists. Defaults to False.

        Returns
        -------
        None
        """
        logger.debug("Entering Dataset.save_dataset_to_pkl")
        if not Path(target_folder).is_dir():
            raise FileNotFoundError(f"{target_folder} does not exist")

        if not filename.endswith(".pkl"):
            filename += ".pkl"

        target_path = Path(target_folder).joinpath(filename)
        if target_path.exists():
            if overwrite:
                target_path.unlink()
            else:
                raise FileExistsError(
                    f"The file {target_path} already exists. Remove it or set overwrite=True"
                )

        with open(target_path, "wb") as outp:
            pickle.dump(self, outp, pickle.HIGHEST_PROTOCOL)

    def import_data_from_file(
        self,
        template_file: Union[str, Path],
        input_data_file: Union[str, Path] = None,
        input_metadata_file: Union[str, Path] = None,
        freq_estimation_method: Literal["highest", "median"] = "median",
        freq_estimation_simplify_tolerance: Union[str, pd.Timedelta] = "2min",
        origin_simplify_tolerance: Union[str, pd.Timedelta] = "5min",
        timestamp_tolerance: Union[str, pd.Timedelta] = "4min",
        kwargs_data_read: dict = {},
        kwargs_metadata_read: dict = {},
        templatefile_is_url: bool = False,
    ) -> None:
        """Import observational data and metadata from files.

        Importing data requires a ´Template´ which is constructed from a template file (JSON).
        (Use ´´metobs_toolkit.build_template_prompt()´´ to create a template file).

        If `input_data_file` is provided, the method reads the raw observational data (CSV).
        A basic quality control (duplicate timestamps and invalid input) is performed, and
        a frequency estimation is made. Based on the estimated frequency, gaps are identified
        if present.

        The method performs the following steps:
        - Estimates the frequency of observations using the ´freq_estimation_method´.
        - Simplifies the estimated frequency and origin timestamps based on tolerances.
        - Alligns the raw timestamps with target timestamps (by origin, and freq) using
          a nearest merge, considering a specified timestamp tolerance.
        - Executes checks for duplicates and invalid input.
        - Identifies gaps in the data.


        if `input_metadata_file` is provided, the method reads the metadata (CSV).

        Parameters
        ------------
        template_file : str or Path
            Path to the template (JSON) file used to interpret the raw data/metadata files.
        input_data_file : str or Path, optional
            Path to the input data file containing observations. If None, no data is read.
        input_metadata_file : str or Path, optional
            Path to the input metadata file. If None, no metadata is read.
        freq_estimation_method : {'highest', 'median'}, optional
            Method to estimate the frequency of observations (per station per observation type).
        freq_estimation_simplify_tolerance : str or pd.Timedelta, optional
            The maximum allowed error in simplifying the target frequency.
        origin_simplify_tolerance : str or pd.Timedelta, optional
            For each time series, the origin (first occurring timestamp) is set and simplification
            is applied.
        timestamp_tolerance : str or pd.Timedelta, optional
            The maximum allowed time shift tolerance for aligning timestamps
            to target (perfect-frequency) timestamps.
        kwargs_data_read : dict, optional
            Additional keyword arguments to pass to `pandas.read_csv()` when
            reading the data file.
        kwargs_metadata_read : dict, optional
            Additional keyword arguments to pass to `pandas.read_csv()` when
            reading the metadata file.
        templatefile_is_url : bool, optional
            If True, the `template_file` is interpreted as a URL to an online
            template file. If False, it is interpreted as a local file path.

        Returns
        -------
        None
        """
        logger.debug("Entering Dataset.import_data_from_file")

        freq_estimation_simplify_tolerance = fmt_timedelta_arg(
            freq_estimation_simplify_tolerance
        )
        origin_simplify_tolerance = fmt_timedelta_arg(origin_simplify_tolerance)
        timestamp_tolerance = fmt_timedelta_arg(timestamp_tolerance)

        if (input_data_file is None) & (template_file is None):
            raise MetObsMissingFile(
                f"No input_data_file or input_metadata_file is provided"
            )

        assert template_file is not None, "No templatefile is specified."

        logger.info(f"Reading the templatefile")
        self.template.read_template_from_file(
            jsonpath=template_file, templatefile_is_url=templatefile_is_url
        )

        if input_data_file is not None:
            use_data = True
            dataparser = DataParser(
                datafilereader=CsvFileReader(file_path=input_data_file),
                template=self.template,
            )
            dataparser.parse(**kwargs_data_read)
        else:
            logger.info("No datafile is provided --> metadata-only mode")
            dataparser = None
            use_data = False

        if input_metadata_file is not None:
            use_metadata = True
            metadataparser = MetaDataParser(
                metadatafilereader=CsvFileReader(file_path=input_metadata_file),
                template=self.template,
            )
            metadataparser.parse(**kwargs_metadata_read)
        else:
            logger.info("No metadatafile is provided.")
            use_metadata = False
            metadataparser = None

        if use_data:
            self.obstypes = update_known_obstype_with_original_data(
                known_obstypes=self.obstypes, template=self.template
            )

        if (use_metadata) & (use_data) & (self.template.data_is_single_station):
            templatename = self.template.single_station_name
            metadataparser._overwrite_name(target_single_name=templatename)

        if use_data:
            stations = createstations(
                data_parser=dataparser,
                metadata_parser=metadataparser,
                use_metadata=use_metadata,
                known_obstypes=self.obstypes,
                timezone=self.template.tz,
                freq_estimation_method=freq_estimation_method,
                freq_estimation_simplify_tolerance=freq_estimation_simplify_tolerance,
                origin_simplify_tolerance=origin_simplify_tolerance,
                timestamp_tolerance=timestamp_tolerance,
            )
        else:
            stations = create_metadata_only_stations(metadata_parser=metadataparser)

        self.stations = stations

    # ------------------------------------------
    #    Plotting
    # ------------------------------------------

    def make_plot_of_modeldata(
        self,
        obstype: str = "temp",
        colormap: Union[dict, None] = None,
        title: Union[str, None] = None,
        linestyle: str = "--",
        ax: Union[Axes, None] = None,
        figkwargs: dict = {},
    ) -> Axes:
        """
        Generate a timeseries plot of model data for a specific observation type.

        Parameters
        ----------
        obstype : str, optional
            The type of observation to plot (e.g., "temp") modeldata for, by default "temp".
        colormap : dict | None,  optional
            The colormap for the lines per stationname. The keys must be the names of the stations,
            and the values the color. If None, a default color map is used.
        title : str or None, optional
            The title of the plot. If None, a default title is generated, by default None.
        linestyle : str, optional
            The style of the line in the plot, by default "--".
        ax : matplotlib.axes.Axes, optional
            The axes on which to plot. If None, a new axes object is created.
        figkwargs : dict, optional
            Additional keyword arguments passed to matplotlib.pyplot.subplots(), by default an empty dictionary.

        Returns
        -------
        matplotlib.axes.Axes
            The axes object containing the plot.
        """
        logger.debug("Entering Dataset.make_plot_of_modeldata")

        modeldatadf = self.modeldatadf
        if obstype not in modeldatadf.index.get_level_values("obstype"):
            raise MetObsObstypeNotFound(f"There is no modeldata present of {obstype}")

        for sta in self.stations:
            if obstype in sta.modeldata.keys():
                modelobstype = sta.modeldata[obstype].obstype
                modelname = sta.modeldata[obstype].modelname
                modelvar = sta.modeldata[obstype].modelname
                break

        if ax is None:
            ax = plotting.create_axes(**figkwargs)

        plotdf = (
            modeldatadf.xs(obstype, level="obstype", drop_level=False)
            .reset_index()
            .set_index(["name", "obstype", "datetime"])
            .sort_index()
        )

        plotdf = plotdf[["value"]]
        plotdf["label"] = label_def["goodrecord"]["label"]

        if colormap is None:
            colormap = plotting.create_categorical_color_map(
                catlist=plotdf.index.get_level_values("name").unique()
            )

        ax = plotting.plot_timeseries_color_by_station(
            plotdf=plotdf,
            colormap=colormap,
            show_outliers=False,
            show_gaps=False,
            ax=ax,
            linestyle=linestyle,
            legend_prefix=f"{modelname}:{modelvar}@",
        )

        if title is None:
            plotting.set_title(
                ax, f"{modelobstype.name} data of {modelname} at stations locations."
            )
        else:
            plotting.set_title(ax, title)

        plotting.set_ylabel(ax, modelobstype._get_plot_y_label())

        cur_tz = plotdf.index.get_level_values("datetime").tz
        plotting.set_xlabel(ax, f"Timestamps (in {cur_tz})")

        plotting.format_datetime_axes(ax)

        plotting.set_legend(ax)

        return ax

    def make_plot(
        self,
        obstype: str = "temp",
        colorby: Literal["station", "label"] = "label",
        show_modeldata: bool = False,
        show_outliers: bool = True,
        show_gaps: bool = True,
        title: Union[str, None] = None,
        ax: Union[Axes, None] = None,
        figkwargs: dict = {},
    ) -> Axes:
        """
        Generate a time series plot for observational data.

        Parameters
        ----------
        obstype : str, optional
            The type of observation to plot (e.g., "temp" for temperature). Default is "temp".
        colorby : {"station", "label"}, optional
            Determines how the data is colored in the plot.
            - "station": Colors by station.
            - "label": Colors by label (the labels refer to the status of a record).
            Default is "label".
        show_modeldata : bool, optional
            If True, includes model data (of the same obstype) if present, in the plot. Default is False.
        show_outliers : bool, optional
            If True, includes outliers (marked by the applied quality control) in the plot. Default is True.
        show_gaps : bool, optional
            If True, gaps are represented by vertical lines in the plot if the gap is unfilled.
            If the gap is filled, it is plotted as a line. Default is True.
        title : str or None, optional
            The title of the plot. If None, a default title is generated. Default is None.
        ax : matplotlib.axes.Axes or None, optional
            The axes on which to draw the plot. If None, a new axes is created. Default is None.
        figkwargs : dict, optional
            Additional keyword arguments passed to matplotlib.pyplot.subplots(), by default an empty dictionary.

        Returns
        -------
        matplotlib.axes.Axes
            The axes containing the plot.
        """
        logger.debug("Entering Dataset.make_plot")

        if ax is None:
            ax = plotting.create_axes(**figkwargs)

        plotdf = (
            self.df.xs(obstype, level="obstype", drop_level=False)
            .reset_index()
            .set_index(["name", "obstype", "datetime"])
            .sort_index()
        )

        if plotdf.empty:
            raise MetObsObstypeNotFound(
                f"There are no records of {obstype} for plotting."
            )

        colormap = None
        if show_modeldata:
            colormap = plotting.create_categorical_color_map(
                catlist=plotdf.index.get_level_values("name").unique()
            )
            ax = self.make_plot_of_modeldata(
                obstype=obstype,
                colormap=colormap,
                ax=ax,
                figkwargs=figkwargs,
                title=title,
            )
        if colorby == "station":
            if colormap is None:
                colormap = plotting.create_categorical_color_map(
                    catlist=plotdf.index.get_level_values("name").unique()
                )
            ax = plotting.plot_timeseries_color_by_station(
                plotdf=plotdf,
                colormap=colormap,
                show_outliers=show_outliers,
                show_gaps=show_gaps,
                ax=ax,
                linestyle="-",
            )

        elif colorby == "label":
            ax = plotting.plot_timeseries_color_by_label(
                plotdf=plotdf,
                show_outliers=show_outliers,
                show_gaps=show_gaps,
                ax=ax,
            )
        else:
            raise ValueError(
                f'colorby is either "station" or "label" but not {colorby}'
            )

        obstypeinstance = self.obstypes[obstype]

        if title is None:
            plotting.set_title(ax, f"{obstypeinstance.name} data.")
        else:
            plotting.set_title(ax, title)

        plotting.set_ylabel(ax, obstypeinstance._get_plot_y_label())

        cur_tz = plotdf.index.get_level_values("datetime").tz
        plotting.set_xlabel(ax, f"Timestamps (in {cur_tz})")

        plotting.format_datetime_axes(ax)

        plotting.set_legend(ax)

        return ax

    def make_gee_plot(
        self,
        geedatasetmanager: Union[
            GEEStaticDatasetManager, GEEDynamicDatasetManager
        ] = default_datasets["LCZ"],
        timeinstance: Union[pd.Timestamp, str, None] = None,
        modelobstype: str = None,
        save: bool = False,
        outputfolder: str = os.getcwd(),
        filename: str = "gee_plot.html",
        vmin: Union[int, float, None] = None,
        vmax: Union[int, float, None] = None,
        overwrite: bool = False,
    ):
        """Create an interactive spatial plot of the GEE dataset and stations.

        This method generates an interactive plot of the GEE dataset and the stations
        are added as markers on it. The interactive plot (a folium map) can be displayed
        in a Jupyter notebook (or Ipython console) or saved as an HTML file.

        If a GEEDynamicDatasetManager is used, the timeinstance and modelobstype
        parameters are required to specify the time and model variable for the plot.



        Parameters
        ----------
        geedatasetmanager : GEEStaticDatasetManager  |  GEEDynamicDatasetManager, optional
            The GEE dataset manager to plot. If a GEEDynamicDatasetManager is provided, a timinstance and modelobstype is required. The default is default_datasets["LCZ"]
        timeinstance :pandas.Timestamp or string or None, optional
            The timinstance to plot the GEE dataset for. This is only used and
            required when geedatasetmanager is a GEEDynamicDatasetManager. The default is None.
        modelobstype : str, optional
            The modelobstype to plot the GEE dataset for. This is only used and
            required when geedatasetmanager is a GEEDynamicDatasetManager. The default is None.
        save : bool, optional
            If True, the plot will be saved as a (HTML) file, that can be opened by a webbrowser. The default is False.
        outputfolder : str, optional
            The folder to save the file too. This is only used when save is True. The default is os.getcwd() (the current working directory).
        filename : str, optional
            The name of the file to save the plot to. This is only used when save is True. The default is "gee_plot.html".
        vmin : int | float | None, optional
            The minimum value for the color scale of the plot. If None, the scale is determined automatically. The default is None.
        vmax : int | float | None, optional
            The maximum value for the color scale of the plot. If None, the scale is determined automatically. The default is None.
        overwrite : bool, optional
            If True, the plot will be overwritten if it already exists. This is only relevant when save is True. The default is False.

        Returns
        -------
        geemap.foliummap.Map
            The interactive map.
        """
        logger.debug("Entering Dataset.make_gee_plot")
        if not isinstance(
            geedatasetmanager, (GEEStaticDatasetManager, GEEDynamicDatasetManager)
        ):
            raise TypeError(
                "geedatasetmanager must be a GEEStaticDatasetManager or GEEDynamicDatasetManager instance."
            )

        if isinstance(geedatasetmanager, GEEStaticDatasetManager):
            kwargs = dict(
                outputfolder=outputfolder,
                filename=filename,
                save=save,
                vmin=vmin,
                vmax=vmax,
                overwrite=overwrite,
            )

        elif isinstance(geedatasetmanager, GEEDynamicDatasetManager):
            if timeinstance is None:
                raise ValueError(
                    f"Timeinstance is None, but is required for a dynamic dataset like {geedatasetmanager}"
                )
            timeinstance = fmt_datetime_arg(timeinstance, tz_if_dt_is_naive="UTC")
            if modelobstype not in geedatasetmanager.modelobstypes:
                raise MetObsObstypeNotFound(
                    f"{modelobstype} is not a known modelobstype of {geedatasetmanager}. These are known: {geedatasetmanager.modelobstypes} "
                )

            kwargs = dict(
                timeinstance=timeinstance,
                modelobstype=modelobstype,
                outputfolder=outputfolder,
                filename=filename,
                save=save,
                vmin=vmin,
                vmax=vmax,
                overwrite=overwrite,
            )

        return geedatasetmanager.make_gee_plot(metadf=self.metadf, **kwargs)

    # ------------------------------------------
    #    Gee extracting
    # ------------------------------------------

    def get_static_gee_point_data(
        self,
        geestaticdatasetmanager: GEEStaticDatasetManager,
        overwrite: bool = True,
        initialize_gee: bool = True,
    ) -> pd.DataFrame:
        """Extract static data from GEE dataset at Stations locations.

        Retrieve Google Earth Engine (GEE) point data of a static dataset, for the stations.
        The retrieved data is stored if overwrite is True.

        Parameters
        ----------
        geestaticdataset : GEEStaticDatasetManager
            An instance of `GEEStaticDatasetManager` representing the static GEE dataset to query.
        overwrite : bool, optional
            If True, the retrieved data will overwrite existing data in the ´Station´'s metadata.
            Default is True.
        initialize_gee : bool, optional
            If True, initializes the GEE API before fetching the data. Default is True.

        Returns
        -------
        pandas.DataFrame
            A DataFrame with the name of the station as index and de extracted value as a column.

        Notes
        -----
        This method interacts with the GEE API to fetch metadata for the station's location.
        Ensure that the GEE API is properly authenticated and initialized before using this method.
        """

        # simple but slow option is to loop over all stations and get the LCZ,
        # but this requires N-API calls (N number of stations).

        # Faster: construct the metadf with all stations, and get the lcs from one api call

        if not isinstance(geestaticdatasetmanager, GEEStaticDatasetManager):
            raise TypeError(
                "geestaticdatasetmanager must be a GEEStaticDatasetManager instance."
            )

        if initialize_gee:
            connect_to_gee()

        geedf = geestaticdatasetmanager.extract_static_point_data(self.metadf)

        varname = geestaticdatasetmanager.name
        if overwrite:
            for staname, geedict in geedf.to_dict(orient="index").items():
                self.get_station(staname).site.set_geedata(varname, geedict[varname])
        return geedf

    def get_static_gee_buffer_fraction_data(
        self,
        geestaticdatasetmanager: GEEStaticDatasetManager,
        buffers: list = [100],
        aggregate: bool = False,
        overwrite: bool = True,
        initialize_gee: bool = True,
    ) -> pd.DataFrame:
        """Extract circular buffer fractions of a GEE dataset at Stations locations.

        Retrieve and optionally store static GEE buffer fraction data for the stations.
        This method interacts with a GEEStaticDatasetManager to fetch buffer fraction
        data for the station's location from Google Earth Engine (GEE). The results can be aggregated, stored, and
        optionally overwrite existing data.

        Parameters
        ----------
        geestaticdataset : GEEStaticDatasetManager
            An instance of GEEStaticDatasetManager used to retrieve static GEE data.
        buffers : list, optional
            A list of buffer radii (in meters) for which to compute the buffer fractions.
            Default is [100].
        aggregate : bool, optional
            If True, aggregate the buffer fraction data. Aggregation schemes are stored per ´GEEStaticDatasetManager´. Default is False.
        overwrite : bool, optional
            If True, overwrites existing fraction data stored in the ´Site´ attribute of the ´Station´s. Default is True.
        initialize_gee : bool, optional
            If True, initialize the GEE environment before retrieving data. Default is True.

        Returns
        -------
        pandas.DataFrame
            A DataFrame containing the buffer fraction data for the stations, with the station names as index.

        Warning
        --------
        This method makes use of GEE API. Make sure that you have access and user rights to use the GEE API.

        Warning
        --------
        It can happen that for stations located on small islands, or close to the coast, the sea-mask is not used as a landcover fraction.

        """
        # simple but slow option is to loop over all stations and get the LCZ,
        # but this requires N-API calls (N number of stations).

        # Faster: construct the metadf with all stations, and get the lcs from one api call

        if not isinstance(geestaticdatasetmanager, GEEStaticDatasetManager):
            raise TypeError(
                "geestaticdatasetmanager must be a GEEStaticDatasetManager instance."
            )

        if initialize_gee:
            connect_to_gee()

        dflist = []
        for bufferradius in buffers:
            geedf = geestaticdatasetmanager.extract_static_buffer_frac_data(
                metadf=self.metadf, bufferradius=bufferradius, agg_bool=aggregate
            )
            dflist.append(geedf)

        geedf = save_concat((dflist))

        if overwrite:
            for staname in geedf.index.get_level_values("name").unique():
                asdict = geedf.loc[staname].to_dict(orient="index")
                for radius, fractions in asdict.items():
                    self.get_station(staname).site.set_gee_buffered_frac_data(
                        buffer=radius, data=fractions
                    )

        return geedf

    def get_LCZ(
        self, overwrite: bool = True, initialize_gee: bool = True
    ) -> pd.DataFrame:
        """
        Retrieve Local Climate Zone (LCZ) for the stations using Google Earth Engine (GEE).

        Parameters
        ----------
        overwrite : bool, optional
            If True, overwrite existing LCZ data if stored in the ´Site´ instances. Default is True.
        initialize_gee : bool, optional
            If True, initialize the Google Earth Engine API before fetching data. Default is True.

        Returns
        -------
        pandas.DataFrame
            A DataFrame containing the LCZ data for the stations, with the station names as index.
        """
        logger.debug("Entering Dataset.get_LCZ")

        return self.get_static_gee_point_data(
            default_datasets["LCZ"],
            overwrite=overwrite,
            initialize_gee=initialize_gee,
        )

    def get_altitude(
        self, overwrite: bool = True, initialize_gee: bool = True
    ) -> pd.DataFrame:
        """
        Retrieve altitude for the stations using Google Earth Engine (GEE).

        Parameters
        ----------
        overwrite : bool, optional
            If True, overwrite existing altitude data if stored in the ´Site´ instances. Default is True.
        initialize_gee : bool, optional
            If True, initialize the Google Earth Engine API before fetching data. Default is True.

        Returns
        -------
        pandas.DataFrame
            A DataFrame containing the altitude data for the stations, with the station names as index.
        """
        logger.debug("Entering Dataset.get_altitude")

        return self.get_static_gee_point_data(
            default_datasets["altitude"],
            overwrite=overwrite,
            initialize_gee=initialize_gee,
        )

    def get_landcover_fractions(
        self,
        buffers: list = [100],
        aggregate: bool = False,
        update_stations: bool = True,
        initialize_gee: bool = True,
    ) -> pd.DataFrame:
        """
        Get landcover fractions for a circular buffer at the stations using GEE.

        Parameters
        ----------
        buffers : list of int, optional
            List of buffer sizes (in meters) to calculate land cover fractions for.
            Default is [100].
        aggregate : bool, optional
            If True, aggregates the data over the buffers. Default is False.
        overwrite : bool, optional
            If True, overwrites existing landcover fraction data stored in the ´Site´ attribute of the ´Station´s. Default is True.

        Returns
        -------
        dict
            A nested dictionary where the keys are buffer radii and the values are the
            corresponding (aggregated) landcoverclasses.

        Warning
        --------
        This method makes use of GEE API. Make sure that you have access and user rights to use the GEE API.

        """
        logger.debug("Entering Dataset.get_landcover_fractions")
        if not isinstance(buffers, list):
            raise TypeError("buffers must be a list.")

        return self.get_static_gee_buffer_fraction_data(
            geestaticdatasetmanager=default_datasets["worldcover"],
            buffers=buffers,
            aggregate=aggregate,
            overwrite=update_stations,
            initialize_gee=initialize_gee,
        )

    @copy_doc(Station.get_gee_timeseries_data)
    def get_gee_timeseries_data(
        self,
        geedynamicdatasetmanager: GEEDynamicDatasetManager,
        startdt_utc: Union[str, pd.Timestamp, None] = None,
        enddt_utc: Union[str, pd.Timestamp, None] = None,
        target_obstypes: list = ["temp"],
        get_all_bands: bool = False,
        drive_filename: Union[str, None] = None,
        drive_folder: str = "gee_timeseries_data",
        force_direct_transfer: bool = False,
        force_to_drive: bool = False,
    ) -> Union[pd.DataFrame, None]:
        logger.debug("Entering Dataset.get_gee_timeseries_data")
        if not isinstance(geedynamicdatasetmanager, GEEDynamicDatasetManager):
            raise TypeError(
                "geedynamicdatasetmanager must be a GEEDynamicDatasetManager instance."
            )
        if not isinstance(target_obstypes, list):
            raise TypeError("target_obstypes must be a list.")

        if startdt_utc is None:
            if self.df.empty:
                raise MetObsMissingArgument(
                    "No data is present in the dataset, thus a startdt_utc is required."
                )
            startdt_utc = self.start_datetime.tz_convert("UTC")
        else:
            startdt_utc = fmt_datetime_arg(startdt_utc, tz_if_dt_is_naive="UTC")

        if enddt_utc is None:
            if self.df.empty:
                raise MetObsMissingArgument(
                    "No data is present in the dataset, thus a enddt_utc is required."
                )
            enddt_utc = self.end_datetime.tz_convert("UTC")
        else:
            enddt_utc = fmt_datetime_arg(enddt_utc, tz_if_dt_is_naive="UTC")

        for obst in target_obstypes:
            if obst not in geedynamicdatasetmanager.modelobstypes.keys():
                raise MetObsMetadataNotFound(
                    f"{obst} is not a known modelobstype of {geedynamicdatasetmanager}."
                )

        if drive_filename is None:
            drive_filename = f"{geedynamicdatasetmanager.name}_timeseries_data_of_full_dataset_{len(self.stations)}_stations"

        df = geedynamicdatasetmanager.extract_timeseries_data(
            metadf=self.metadf,
            startdt_utc=startdt_utc,
            enddt_utc=enddt_utc,
            obstypes=target_obstypes,
            get_all_bands=get_all_bands,
            drive_filename=drive_filename,
            drive_folder=drive_folder,
            force_direct_transfer=force_direct_transfer,
            force_to_drive=force_to_drive,
        )
        if df is None:
            logger.warning("No data is returned by the GEE api request.")
            return

        # Import the data and add it to the Stations

        # Note: all functionallity of this part is already implemented in
        # import_gee_data_from_file, except for the actual reading/formatting and unit converting.
        # Thus we call that function, by providing it to _force_from_dataframe, these steps are skipped.
        _ = self.import_gee_data_from_file(
            filepath=None,
            geedynamicdatasetmanager=geedynamicdatasetmanager,
            force_update=True,
            _force_from_dataframe=df,
        )
        return df

    # ------------------------------------------
    #    QC
    # ------------------------------------------
    _use_mp_docargstr = """use_mp: bool, optional
        If True, the function will use multiprocessing to speed up the calculations.
        The default is False."""

    @copy_doc(copy_func=Station.gross_value_check, extra_param_desc=_use_mp_docargstr)
    def gross_value_check(
        self,
        target_obstype: str = "temp",
        lower_threshold: float = -15.0,
        upper_threshold: float = 39.0,
        use_mp: bool = True,
    ) -> None:
        logger.debug("Entering Dataset.gross_value_check")

        func_feed_list = _create_qc_arg_set(
            dataset=self,
            target_obstype=target_obstype,
            lower_threshold=lower_threshold,
            upper_threshold=upper_threshold,
        )
        if use_mp:
            with concurrent.futures.ProcessPoolExecutor() as executor:
                stationgenerator = executor.map(
                    _qc_grossvalue_generatorfunc, func_feed_list
                )
            self.stations = list(stationgenerator)
        else:
            self.stations = list(map(_qc_grossvalue_generatorfunc, func_feed_list))

    @copy_doc(copy_func=Station.persistence_check, extra_param_desc=_use_mp_docargstr)
    def persistence_check(
        self,
        target_obstype: str = "temp",
        timewindow: Union[str, pd.Timedelta] = pd.Timedelta("60min"),
        min_records_per_window: int = 5,
        use_mp: bool = True,
    ) -> None:
        logger.debug("Entering Dataset.persistence_check")
        timewindow = fmt_timedelta_arg(timewindow)
        func_feed_list = _create_qc_arg_set(
            dataset=self,
            target_obstype=target_obstype,
            timewindow=timewindow,
            min_records_per_window=min_records_per_window,
        )
        if use_mp:
            with concurrent.futures.ProcessPoolExecutor() as executor:
                stationgenerator = executor.map(
                    _qc_persistence_generatorfunc, func_feed_list
                )
            self.stations = list(stationgenerator)
        else:
            self.stations = list(map(_qc_persistence_generatorfunc, func_feed_list))

    @copy_doc(copy_func=Station.repetitions_check, extra_param_desc=_use_mp_docargstr)
    def repetitions_check(
        self,
        target_obstype: str = "temp",
        max_N_repetitions: int = 5,
        use_mp: bool = True,
    ) -> None:
        logger.debug("Entering Dataset.repetitions_check")

        func_feed_list = _create_qc_arg_set(
            dataset=self,
            target_obstype=target_obstype,
            max_N_repetitions=max_N_repetitions,
        )

        if use_mp:
            with concurrent.futures.ProcessPoolExecutor() as executor:
                stationgenerator = executor.map(
                    _qc_repetitions_generatorfunc, func_feed_list
                )
            self.stations = list(stationgenerator)
        else:
            self.stations = list(map(_qc_repetitions_generatorfunc, func_feed_list))

    @copy_doc(copy_func=Station.step_check, extra_param_desc=_use_mp_docargstr)
    def step_check(
        self,
        target_obstype: str = "temp",
        max_increase_per_second: Union[int, float] = 8.0 / 3600.0,
        max_decrease_per_second: Union[int, float] = -10.0 / 3600.0,
        use_mp: bool = True,
    ) -> None:
        logger.debug("Entering Dataset.step_check")

        func_feed_list = _create_qc_arg_set(
            dataset=self,
            target_obstype=target_obstype,
            max_increase_per_second=max_increase_per_second,
            max_decrease_per_second=max_decrease_per_second,
        )

        if use_mp:
            with concurrent.futures.ProcessPoolExecutor() as executor:
                stationgenerator = executor.map(_qc_step_generatorfunc, func_feed_list)
            self.stations = list(stationgenerator)
        else:
            self.stations = list(map(_qc_step_generatorfunc, func_feed_list))

    @copy_doc(
        copy_func=Station.window_variation_check, extra_param_desc=_use_mp_docargstr
    )
    def window_variation_check(
        self,
        target_obstype: str = "temp",
        timewindow: Union[str, pd.Timedelta] = pd.Timedelta("1h"),
        min_records_per_window: int = 3,
        max_increase_per_second: Union[int, float] = 0.0022,
        max_decrease_per_second: Union[int, float] = -0.0027,
        use_mp: bool = True,
    ) -> None:
        logger.debug("Entering Dataset.window_variation_check")

        timewindow = fmt_timedelta_arg(timewindow)

        func_feed_list = _create_qc_arg_set(
            dataset=self,
            target_obstype=target_obstype,
            timewindow=timewindow,
            min_records_per_window=min_records_per_window,
            max_increase_per_second=max_increase_per_second,
            max_decrease_per_second=max_decrease_per_second,
        )

        if use_mp:
            with concurrent.futures.ProcessPoolExecutor() as executor:
                stationgenerator = executor.map(
                    _qc_window_var_generatorfunc, func_feed_list
                )
            self.stations = list(stationgenerator)
        else:
            self.stations = list(map(_qc_window_var_generatorfunc, func_feed_list))

    def buddy_check(
        self,
        target_obstype: str = "temp",
        buddy_radius: Union[int, float] = 10000,
        min_sample_size: int = 4,
        max_alt_diff: Union[int, float, None] = None,
        min_std: Union[int, float] = 1.0,
        std_threshold: Union[int, float] = 3.1,
        N_iter: int = 2,
        instantanious_tolerance: Union[str, pd.Timedelta] = pd.Timedelta("4min"),
        lapserate: Union[float, None] = None,  # -0.0065 for temperature (in °C)
        use_mp: bool = True,
    ):
        """Spatial buddy check.

        The buddy check compares an observation against its neighbors (i.e. buddies). The check loops
        over all the groups, which are stations within a radius of each other. For each group, the absolute
        value of the difference with the groupmean, normalized by the standared deviation (with a defined minimum),
        is computed. If one (or more) exeeds the std_theshold, the most extreme (=baddest) observation of that group is labeled as an outlier.

        Multiple iterations of this checks can be done using the N_iter.

        A schematic step-by-step description of the buddy check:

        1. A distance matrix is constructed for all interdistances between the stations. This is done using the haversine approximation.
        2. Groups of buddies (neighbours) are created by using the buddy_radius. These groups are further filtered by:
            * removing stations from the groups that differ to much in altitude (based on the max_alt_diff)
            * removing groups of buddies that are too small (based on the min_sample_size)
        3. Observations per group are synchronized in time (using the max_shift as tolerance for allignment).
        4. If a lapsrate is specified, the observations are corrected for altitude differences.
        5. For each buddy group:
            * The mean, standard deviation (std), and sample size are computed.
            * If the std is lower than the minimum std, it is replaced by the minimum std.
            * Chi values are calculated for all records.
            * For each timestamp the record with the highest Chi is tested if it is larger then std_threshold.
            If so, that record is flagged as an outlier. It will be ignored in the next iteration.
            * This is repeated N_iter times.

        Parameters
        ----------
        target_obstype : str, optional
            The target observation to check. Default is "temp".
        buddy_radius : int | float, optional
            The radius to define spatial neighbors in meters. Default is 10000.
        min_sample_size : int, optional
            The minimum sample size to calculate statistics on. Default is 4.
        max_alt_diff : int | float | None, optional
            The maximum altitude difference allowed for buddies. Default is None.
        min_std : int | float, optional
            The minimum standard deviation for sample statistics. Default is 1.0.
        std_threshold : int | float, optional
            The threshold (std units) for flagging observations as outliers. Default is 3.1.
        N_iter : int, optional
            The number of iterations to perform the buddy check. Default is 2.
        instantanious_tolerance : str | pd.Timedelta, optional
            The maximum time difference allowed for synchronizing observations. Default is pd.Timedelta("4min").
        lapserate : int | float | None, optional
            Describe how the obstype changes with altitude (in meters). Default is None.
        use_mp : bool, optional
            Use multiprocessing to speed up the buddy check. Default is True.

        Returns
        -------
        None

        Notes
        ------
        - This method modifies the outliers in place and does not return anything.
          You can use the `outliersdf` property to view all flagged outliers.
        - The altitude of the stations can be extracted from GEE by using the `Dataset.get_altitude()` method.

        """
        logger.debug("Entering Dataset.buddy_check")

        instantanious_tolerance = fmt_timedelta_arg(instantanious_tolerance)
        if (lapserate is not None) | (max_alt_diff is not None):
            if not all([sta.site.flag_has_altitude() for sta in self.stations]):
                raise MetObsMetadataNotFound(
                    "Not all stations have altitude data, lapserate correction and max_alt_diff filtering could not be applied."
                )

        qc_kwargs = dict(
            obstype=target_obstype,
            buddy_radius=buddy_radius,
            min_sample_size=min_sample_size,
            max_alt_diff=max_alt_diff,
            min_std=min_std,
            std_threshold=std_threshold,
            N_iter=N_iter,
            instantanious_tolerance=instantanious_tolerance,
            lapserate=lapserate,
            use_mp=use_mp,
        )

        outlierslist, timestamp_map = toolkit_buddy_check(dataset=self, **qc_kwargs)
        # outlierslist is a list of tuples (stationname, datetime, msg) that are outliers
        # timestamp_map is a dict with keys the stationname and values a series to map the syncronized
        # timestamps to the original timestamps

        # convert to a dataframe
        df = pd.DataFrame(data=outlierslist, columns=["name", "datetime", "detail_msg"])

        # update all the sensordata
        for station in self.stations:
            if target_obstype in station.sensordata.keys():
                # Get the sensordata object
                sensorddata = station.get_sensor(target_obstype)

                # get outlier datetimeindex
                outldt = pd.DatetimeIndex(df[df["name"] == station.name]["datetime"])

                if not outldt.empty:
                    # convert to original timestamps
                    dtmap = timestamp_map[station.name]
                    outldt = outldt.map(dtmap)

                # update the sensordata
                sensorddata._update_outliers(
                    qccheckname="buddy_check",
                    outliertimestamps=outldt,
                    check_kwargs=qc_kwargs,
                    extra_columns={
                        "detail_msg": df[df["name"] == station.name][
                            "detail_msg"
                        ].to_numpy()
                    },
                )

    @copy_doc(Station.get_qc_stats)
    def get_qc_stats(
        self, target_obstype: str = "temp", make_plot: bool = True
    ) -> Union[pd.DataFrame, None]:
        logger.debug("Entering Dataset.get_qc_stats")
        freqdf_list = [
            sta.get_qc_stats(target_obstype=target_obstype, make_plot=False)
            for sta in self.stations
        ]

        dfagg = (
            pd.concat(freqdf_list)
            .reset_index()
            .groupby(["qc_check"])
            .sum()
            .drop(columns=["name"])
        )

        if make_plot:
            fig = plotting.qc_overview_pies(df=dfagg)
            fig.suptitle(
                f"QC frequency statistics of {target_obstype} on Dataset level."
            )
            return fig
        else:
            return dfagg

    # ------------------------------------------
    #    Other methods
    # ------------------------------------------

    def rename_stations(self, renamedict: dict) -> None:
        """
        Rename stations in the dataset.

        Parameters
        ----------
        renamedict : dict
            A dictionary where keys are the original station names and values are the new station names.

        Returns
        -------
        None

        Warnings
        --------
        - If a station name in `renamedict` does not exist in the dataset, it will be skipped.
        - If a target station name in `renamedict` already exists in the dataset, the renaming
          operation for that station will be skipped.

        """
        logger.debug("Entering Dataset.rename_stations")

        for origname, trgname in renamedict.items():
            if origname not in [sta.name for sta in self.stations]:
                logger.warning(f"{origname} is not present, skipped.")
                continue

            if trgname in [sta.name for sta in self.stations]:
                logger.warning(
                    f"{trgname} is already present, renaming {origname} --> {trgname} is skipped."
                )
                continue
            self.get_station(origname)._rename(targetname=trgname)

    @copy_doc(Station.convert_outliers_to_gaps)
    def convert_outliers_to_gaps(
        self, all_observations: bool = True, obstype: str = "temp"
    ) -> None:
        logger.debug("Entering Dataset.convert_outliers_to_gaps")

        for sta in self.stations:
            sta.convert_outliers_to_gaps(
                all_observations=all_observations, obstype=obstype
            )

    # ------------------------------------------
    #    Gapfilling
    # ------------------------------------------
    @copy_doc(Station.interpolate_gaps)
    def interpolate_gaps(
        self,
        target_obstype: str,
        method: str = "time",
        max_consec_fill: int = 10,
        n_leading_anchors: int = 1,
        n_trailing_anchors: int = 1,
        max_lead_to_gap_distance: Union[pd.Timedelta, None] = None,
        max_trail_to_gap_distance: Union[pd.Timedelta, None] = None,
        overwrite_fill: bool = False,
        method_kwargs: dict = {},
    ) -> None:
        logger.debug("Entering Dataset.interpolate_gaps")

        max_lead_to_gap_distance = fmt_timedelta_arg(max_lead_to_gap_distance)
        max_trail_to_gap_distance = fmt_timedelta_arg(max_trail_to_gap_distance)

        for sta in self.stations:
            sta.interpolate_gaps(
                target_obstype=target_obstype,
                method=method,
                max_consec_fill=max_consec_fill,
                n_leading_anchors=n_leading_anchors,
                n_trailing_anchors=n_trailing_anchors,
                max_lead_to_gap_distance=max_lead_to_gap_distance,
                max_trail_to_gap_distance=max_trail_to_gap_distance,
                overwrite_fill=overwrite_fill,
                method_kwargs=method_kwargs,
            )

    @copy_doc(Station.fill_gaps_with_raw_modeldata)
    def fill_gaps_with_raw_modeldata(
        self, target_obstype: str, overwrite_fill: bool = False
    ) -> None:
        logger.debug("Entering Dataset.fill_gaps_with_raw_modeldata")

        for sta in self.stations:
            sta.fill_gaps_with_raw_modeldata(
                target_obstype=target_obstype, overwrite_fill=overwrite_fill
            )

    @copy_doc(Station.fill_gaps_with_debiased_modeldata)
    def fill_gaps_with_debiased_modeldata(
        self,
        target_obstype: str,
        leading_period_duration: Union[str, pd.Timedelta] = pd.Timedelta("24h"),
        min_leading_records_total: int = 60,
        trailing_period_duration: Union[str, pd.Timedelta] = pd.Timedelta("24h"),
        min_trailing_records_total: int = 60,
        overwrite_fill: bool = False,
    ) -> None:
        logger.debug("Entering Dataset.fill_gaps_with_debiased_modeldata")

        leading_period_duration = fmt_timedelta_arg(leading_period_duration)
        trailing_period_duration = fmt_timedelta_arg(trailing_period_duration)

        for sta in self.stations:
            sta.fill_gaps_with_debiased_modeldata(
                target_obstype=target_obstype,
                leading_period_duration=leading_period_duration,
                min_leading_records_total=min_leading_records_total,
                trailing_period_duration=trailing_period_duration,
                min_trailing_records_total=min_trailing_records_total,
                overwrite_fill=overwrite_fill,
            )

    @copy_doc(Station.fill_gaps_with_diurnal_debiased_modeldata)
    def fill_gaps_with_diurnal_debiased_modeldata(
        self,
        target_obstype: str,
        leading_period_duration: Union[str, pd.Timedelta] = pd.Timedelta("24h"),
        trailing_period_duration: Union[str, pd.Timedelta] = pd.Timedelta("24h"),
        min_debias_sample_size: int = 6,
        overwrite_fill: bool = False,
    ) -> None:
        logger.debug("Entering Dataset.fill_gaps_with_diurnal_debiased_modeldata")

        leading_period_duration = fmt_timedelta_arg(leading_period_duration)
        trailing_period_duration = fmt_timedelta_arg(trailing_period_duration)

        for sta in self.stations:
            sta.fill_gaps_with_diurnal_debiased_modeldata(
                target_obstype=target_obstype,
                leading_period_duration=leading_period_duration,
                trailing_period_duration=trailing_period_duration,
                min_debias_sample_size=min_debias_sample_size,
                overwrite_fill=overwrite_fill,
            )

    @copy_doc(Station.fill_gaps_with_weighted_diurnal_debiased_modeldata)
    def fill_gaps_with_weighted_diurnal_debiased_modeldata(
        self,
        target_obstype: str,
        leading_period_duration: Union[str, pd.Timedelta] = pd.Timedelta("24h"),
        trailing_period_duration: Union[str, pd.Timedelta] = pd.Timedelta("24h"),
        min_lead_debias_sample_size: int = 2,
        min_trail_debias_sample_size: int = 2,
        overwrite_fill: bool = False,
    ) -> None:
        logger.debug(
            "Entering Dataset.fill_gaps_with_weighted_diurnal_debiased_modeldata"
        )

        leading_period_duration = fmt_timedelta_arg(leading_period_duration)
        trailing_period_duration = fmt_timedelta_arg(trailing_period_duration)

        for sta in self.stations:
            sta.fill_gaps_with_weighted_diurnal_debiased_modeldata(
                target_obstype=target_obstype,
                leading_period_duration=leading_period_duration,
                trailing_period_duration=trailing_period_duration,
                min_lead_debias_sample_size=min_lead_debias_sample_size,
                min_trail_debias_sample_size=min_trail_debias_sample_size,
                overwrite_fill=overwrite_fill,
            )


def _qc_grossvalue_generatorfunc(input: list) -> Station:
    station, kwargs = input
    station.gross_value_check(**kwargs)
    return station


def _qc_persistence_generatorfunc(input: list) -> Station:
    station, kwargs = input
    station.persistence_check(**kwargs)
    return station


def _qc_repetitions_generatorfunc(input: list) -> Station:
    station, kwargs = input
    station.repetitions_check(**kwargs)
    return station


def _qc_step_generatorfunc(input: list) -> Station:
    station, kwargs = input
    station.step_check(**kwargs)
    return station


def _qc_window_var_generatorfunc(input: list) -> Station:
    station, kwargs = input
    station.window_variation_check(**kwargs)
    return station


# ------------------------------------------
#    Helping functions
# ------------------------------------------


def _create_qc_arg_set(dataset: Dataset, **qckwargs: dict) -> list:
    return [([sta, dict(**qckwargs)]) for sta in dataset.stations]


def create_metadata_only_stations(metadata_parser: MetaDataParser) -> list:
    """
    Create a list of Station objects from metadata.

    Parameters
    ----------
    metadata_parser : MetaDataParser
        An object that provides access to the metadata for the stations.

    Returns
    -------
    list
        A list of Station objects, each representing a station with its associated metadata.
    """
    logger.debug("Entering create_metadata_only_stations")

    stations = []

    for stationname in metadata_parser.get_df().index:
        stationsite = Site(
            stationname=stationname,
            latitude=metadata_parser.get_station_lat(stationname),
            longitude=metadata_parser.get_station_lon(stationname),
            extradata=metadata_parser.get_station_extra_metadata(stationname),
        )

        station = Station(
            stationname=stationname,
            site=stationsite,
            all_sensor_data=[],
        )

        stations.append(station)
    return stations


def createstations(
    data_parser: DataParser,
    metadata_parser: Union[MetaDataParser, None],
    use_metadata: bool,
    known_obstypes: dict,
    timezone: str,
    freq_estimation_method: str,
    freq_estimation_simplify_tolerance: Union[pd.Timedelta, str],
    origin_simplify_tolerance: Union[pd.Timedelta, str],
    timestamp_tolerance: Union[pd.Timedelta, str],
) -> list:
    """
    Create a list of Station objects from parsed data and metadata.

    Parameters
    ----------
    data_parser : DataParser
        An object that provides access to the observational data.
    metadata_parser : MetadataParser or None
        An object that provides access to the metadata for the stations.
    use_metadata : bool
        Whether to use metadata for creating stations.
    known_obstypes : dict
        A dictionary mapping observation type names to their corresponding ObsType objects.
    timezone : str
        A pytz equivalent string indicating the timezone of the timestamps.
    freq_estimation_method : str
        Method to estimate the frequency of observations (per station per observation type).
    freq_estimation_simplify_tolerance : pd.Timedelta or str
        The maximum allowed error in simplifying the target frequency.
    origin_simplify_tolerance : pd.Timedelta or str
        For each time series, the origin (first occurring timestamp) is set and simplification
        is applied.
    timestamp_tolerance : pd.Timedelta or str
        The maximum allowed time shift tolerance for aligning timestamps
        to target (perfect-frequency) timestamps.

    Returns
    -------
    list
        A list of Station objects, each representing a station with its associated sensor data and metadata.
    """
    logger.debug("Entering createstations")
    datadf = data_parser.get_df()

    not_an_obstype = ["name", "datetime"]
    stations = []
    for stationname, stationdata in datadf.groupby("name"):
        all_station_sensor_data = []
        for obstypename in stationdata.columns:
            if obstypename in not_an_obstype:
                continue

            obstype = known_obstypes[obstypename]

            # Formatting on raw data
            # 1. Skip stations with nan as name (caused by string casting errors)
            if (stationname == str(np.nan)) | (pd.isnull(stationname)):
                logger.warning(
                    "Skipping the records belonging to station with Nan as name. This could be the result from stringcasting the stationnames."
                )
                continue
            # 2. Drop NAT datetimes if present
            stationdata = stationdata.loc[pd.notnull(stationdata["datetime"])]

            # 3. Skip stations if there are less than 2 records (no freq can be estimated)
            if stationdata.shape[0] < 2:
                logger.warning(
                    f"Station {stationname} is skipped because it has only one record."
                )
                continue

            # Get dataseries:
            sensordata = SensorData(
                stationname=stationname,
                datarecords=stationdata[obstype.name].to_numpy(),
                timestamps=stationdata["datetime"].to_numpy(),
                obstype=obstype,
                timezone=timezone,
                freq_estimation_method=freq_estimation_method,
                freq_estimation_simplify_tolerance=freq_estimation_simplify_tolerance,
                origin_simplify_tolerance=origin_simplify_tolerance,
                timestamp_tolerance=timestamp_tolerance,
            )

            # Add Sensordata:
            all_station_sensor_data.append(sensordata)

        # Create a Site object for the station
        if use_metadata:
            stationsite = Site(
                stationname=stationname,
                latitude=metadata_parser.get_station_lat(stationname),
                longitude=metadata_parser.get_station_lon(stationname),
                extradata=metadata_parser.get_station_extra_metadata(stationname),
            )
        else:
            # no metafile is provided
            stationsite = Site(
                stationname=stationname, latitude=np.nan, longitude=np.nan, extradata={}
            )

        # Combine into a Station
        station = Station(
            stationname=stationname,
            site=stationsite,
            all_sensor_data=all_station_sensor_data,
        )

        stations.append(station)

    if use_metadata:
        # Check the stations present in the metadata but not in the data
        missing_in_data = list(
            set(metadata_parser.get_df().index) - set(datadf["name"].unique())
        )
        if bool(missing_in_data):
            logger.warning(
                f"The following stations are defined in the metadatafile but no records are found in the data:\n {missing_in_data}"
            )

    return stations


def import_dataset_from_pkl(target_path: Union[str, Path]) -> Dataset:
    """
    Import a Dataset instance from a pickle file.

    Parameters
    ----------
    target_path : str or Path
        The path to the pickle file.

    Returns
    -------
    Dataset
        The Dataset instance.
    """
    logger.debug("Entering import_dataset_from_pkl")

    picklereader = PickleFileReader(file_path=target_path)
    return picklereader.read_as_local_file()
