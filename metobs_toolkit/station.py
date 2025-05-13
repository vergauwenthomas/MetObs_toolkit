import logging
import numpy as np
import pandas as pd
from typing import Literal, Union
from datetime import datetime
from matplotlib.pyplot import Axes

from metobs_toolkit.site import Site
from metobs_toolkit.backend_collection.argumentcheckers import (
    fmt_datetime_arg,
    fmt_timedelta_arg,
)
import metobs_toolkit.plot_collection as plotting

from metobs_toolkit.backend_collection.errorclasses import *
import metobs_toolkit.backend_collection.printing_collection as printing
from metobs_toolkit.backend_collection.df_helpers import save_concat
from metobs_toolkit.settings_collection import label_def
from metobs_toolkit.backend_collection.dev_collection import copy_doc
import metobs_toolkit.qc_collection as qc_collection
from metobs_toolkit.geedatasetmanagers import (
    GEEStaticDatasetManager,
    GEEDynamicDatasetManager,
)
from metobs_toolkit.geedatasetmanagers import default_datasets as default_gee_datasets
from metobs_toolkit.modeltimeseries import ModelTimeSeries

logger = logging.getLogger("<metobs_toolkit>")


class Station:
    """
    Represents a weather station, holding metadata, sensor data, and model data.

    Parameters
    ----------
    stationname : str
        Name of the station.
    site : Site
        Site instance containing metadata and location.
    all_sensor_data : list
        List of SensorData instances for the station.
    """

    def __init__(self, stationname: str, site: Site, all_sensor_data: list):
        # dimension attributes
        self._name = str(stationname)
        self._site = site
        self.obsdata = self._set_stationdata(
            all_sensor_data
        )  # obstypename : SensorData

        # Extra extracted data
        self._modeldata = {}  # dict of ModelTimeSeries

    def __eq__(self, other):
        """Check equality with another Station object."""
        if not isinstance(other, Station):
            return False
        return (
            self.name == other.name
            and self.site == other.site
            and self.obsdata == other.obsdata
            and self._modeldata == other._modeldata
        )

    def __repr__(self):
        """Return a string representation for debugging."""
        return f"Station instance of {self.name}"

    @property
    def name(self) -> str:
        """The name of the station."""
        return str(self._name)

    @property
    def site(self) -> "Site":
        """The Site instance of the station."""
        return self._site

    @property
    def sensordata(self) -> dict:
        """The SensorData related to the station, as a dictionary."""
        return dict(self.obsdata)

    def get_sensor(self, obstype: str) -> "SensorData":
        """Get the SensorData instance for a specific observation type.

        Parameters
        ----------
        obstype : str
            The observation type to retrieve.

        Returns
        -------
        SensorData
            The SensorData instance for the specified observation type.
        """
        self._obstype_is_known_check(obstype)
        return self.obsdata[obstype]

    @property
    def df(self) -> pd.DataFrame:
        """Construct a DataFrame representation of the observations.

        Returns
        --------
            pd.DataFrame: A pandas DataFrame with a single column 'value'.
        """

        # return dataframe with ['datetime', 'obstype'] as index and 'value' as single column.
        concatdf = save_concat(([sensor.df for sensor in self.sensordata.values()]))

        # sort by datetime
        concatdf.sort_index(inplace=True)
        return concatdf

    @property
    def outliersdf(self) -> pd.DataFrame:
        """Construct a DataFrame representation of all the outliers.

        Outliers are the observations that are flagged by the performed quality control.


        Returns
        --------
        pandas.DataFrame:
            A DataFrame with two columns ['value', 'label'], representing
            the value and details of the flagged observation.

        """

        concatlist = []
        for sensordata in self.sensordata.values():
            stadf = sensordata.outliersdf[["value", "label"]].reset_index()
            stadf["obstype"] = sensordata.obstype.name
            concatlist.append(stadf.set_index(["datetime", "obstype"]))

        combdf = save_concat((concatlist))
        combdf.sort_index(inplace=True)
        if combdf.empty:

            combdf = pd.DataFrame(
                columns=["value", "label"],
                index=pd.MultiIndex(
                    levels=[[], []], codes=[[], []], names=["datetime", "obstype"]
                ),
            )
        return combdf

    @property
    def gapsdf(self) -> pd.DataFrame:
        """Construct a DataFrame representation of all the gaps.

        Returns
        --------
        pandas.DataFrame:
            A DataFrame with columns ['value', 'label', 'details'], representing
            the value, the gap-label and details of the gap record.

        """
        concatlist = []
        for sensordata in self.sensordata.values():
            stadf = sensordata.gapsdf.reset_index()
            stadf["obstype"] = sensordata.obstype.name
            concatlist.append(stadf.set_index(["datetime", "obstype"]))

        combdf = save_concat(concatlist)
        combdf.sort_index(inplace=True)
        if combdf.empty:
            combdf = pd.DataFrame(
                columns=["value", "label", "details"],
                index=pd.MultiIndex(
                    levels=[[], []], codes=[[], []], names=["datetime", "obstype"]
                ),
            )

        return combdf

    @property
    def metadf(self) -> pd.DataFrame:
        """Construct a DataFrame representation of metadata.

        Metadata is the information related to the sensors, that does not change over time.
        The metadata is extracted from the site instance.

        Returns
        --------
        pandas.DataFrame:
            A DataFrame with the station names as index, and the metadata as columns.

        """

        return self.site.metadf

    @property
    def modeldatadf(self) -> pd.DataFrame:
        """Construct a DataFrame representation of all the present modeldata.

        Modeldata is stored as ´ModelTimeSeries´ instances, and is set as an attribute of
        a ´Station´.

        Returns
        --------
        pandas.DataFrame:
            A DataFrame with columns ['value', 'details'], representing
            the value, and details of the corresponding modeldata.

        """
        concatlist = []
        for modeldata in self.modeldata.values():
            df = (
                modeldata.df.assign(obstype=modeldata.obstype.name)
                .assign(
                    details=f"{modeldata.modelname}:{modeldata.modelvariable} converted from {modeldata.obstype.model_unit} -> {modeldata.obstype.std_unit}"
                )
                .reset_index()
                .set_index(["datetime", "obstype"])
            )
            concatlist.append(df)
        combdf = save_concat(concatlist)
        if combdf.empty:
            combdf = pd.DataFrame(
                columns=["value", "details"],
                index=pd.MultiIndex(
                    levels=[[], []], codes=[[], []], names=["datetime", "obstype"]
                ),
            )
        # formatting
        combdf = combdf[["value", "details"]]
        combdf.sort_index(inplace=True)
        return combdf

    @property
    def start_datetime(self) -> pd.Timestamp:
        """Get the earliest start datetime from the observation data.

        Returns
        -------
        pd.Timestamp
            The earliest start datetime among all sensor data.
        """

        return min([sensdata.start_datetime for sensdata in self.sensordata.values()])

    @property
    def end_datetime(self) -> pd.Timestamp:
        """Get the latest end datetime from the observation data.

        Returns
        -------
        pd.Timestamp
            The maximum end datetime across all sensor data in the observation data.
        """
        return max([sensdata.end_datetime for sensdata in self.sensordata.values()])

    @property
    def modeldata(self) -> dict:
        """
        Retrieve the model data associated with the station.

        Returns
        -------
        dict
            A dictionary with the observation type as key and the corresponding
            ModelTimeSeries as values.
        """
        return self._modeldata

    @property
    def present_observations(self) -> list:
        """Get a list of all the present observation types.

        Returns
        -------
        list
            A list of all the present observations in the station.
        """
        return sorted(list(self.sensordata.keys()))

    def add_to_modeldata(
        self, new_modeltimeseries: ModelTimeSeries, force_update: bool = False
    ) -> None:
        """
        Add a new ModelTimeSeries to the Station.

        Parameters
        ----------
        new_modeltimeseries : ModelTimeSeries
            The new model time series to be added. Must be an instance of `ModelTimeSeries`.
        force_update : bool, optional
            If True, overwrite existing model data for the same observation type. Default is False.

        Returns
        -------
        None
        """
        logger.debug("Entering add_to_modeldata for %s", self)
        # Validate argument types
        if not isinstance(new_modeltimeseries, ModelTimeSeries):
            raise TypeError(
                "new_modeltimeseries must be an instance of ModelTimeSeries."
            )
        if not isinstance(force_update, bool):
            raise TypeError("force_update must be a boolean.")

        # Test if there is already modeldata for the same obstype available
        if (new_modeltimeseries.obstype.name in self.sensordata) & (not force_update):
            raise MetObsDataAlreadyPresent(
                f"There is already an modeltimeseries instance represinting {new_modeltimeseries.obstype.name}, and force_update is False."
            )

        self._modeldata.update({new_modeltimeseries.obstype.name: new_modeltimeseries})

    def get_info(self, printout: bool = True) -> Union[str, None]:
        """
        Retrieve and optionally print detailed information about the station.

        Parameters
        ----------
        printout : bool, optional
            If True, prints the information to the console. If False, returns
            the information as a string. Default is True.

        Returns
        -------
        str or None
            A string containing the station information if `printout` is False.
            Otherwise, returns None.
        """
        logger.debug("Entering get_info for %s", self)

        infostr = ""
        infostr += printing.print_fmt_title("General info of Station")

        # --- Observational info ---
        infostr += printing.print_fmt_section("Observational info")
        df = self.df
        if df.empty:
            infostr += printing.print_fmt_line(
                "Station instance without observation records."
            )
        else:
            present_obstypes = list(df.index.get_level_values("obstype").unique())

            infostr += printing.print_fmt_line("Station instance with:", 0)
            for obstype in present_obstypes:
                infostr += printing.print_fmt_line(f"{obstype}:", 1)
                infostr += self.get_sensor(obstype)._get_info_core(nident_root=2)

        # Meta data info
        infostr += printing.print_fmt_section("Metadata info")
        infostr += self.site._get_info_core(nident_root=1)

        # Modeldata info
        infostr += printing.print_fmt_section("Modeldata info")
        if not bool(self.modeldata):
            infostr += printing.print_fmt_line("Station instance without modeldata.")
        else:

            for obstype, modeldata in self.modeldata.items():
                infostr += printing.print_fmt_line(f"{obstype}:", 1)
                infostr += modeldata._get_info_core(nident_root=2)

        if printout:
            print(infostr)
        else:
            return infostr

    def resample(
        self,
        target_freq: Union[str, pd.Timedelta],
        target_obstype: Union[str, None] = None,
        shift_tolerance: Union[str, pd.Timedelta] = pd.Timedelta("4min"),
        origin: Union[None, pd.Timestamp] = None,
        origin_simplify_tolerance: Union[str, pd.Timedelta] = pd.Timedelta("4min"),
    ) -> None:
        """
        Resample observation data to a specified frequency.

        Resampling is done by creating target timestamps, aligning them with the present timestamps,
        and transfer the nearest (in time) value to the target timestamp. There is no interpolation in time.

        Aligning restrictions can be specified by setting the shif_tolerance, which indicates the maximum shift
        (in time) that alligned timestamps can have.

        Since the origin (the first timestamp in the current observations) is in general not a suitable candidate
        to serve as origin for the target frequency, a suitable origin is deduced. This can be done by specifying
        the origin argument, or if None it is done automatically. This is done by simplifying the
        current origin by rounding it down, and testing if the shift in the new tolerance is smaller than origin_simplify_tolerance.


        Parameters
        ----------
        target_freq : str or pandas.Timedelta
            The target frequency to which the data should be resampled. Can be
            specified as a pandas frequency string (e.g., '5T' for 5 minutes)
            or a pandas.Timedelta object.
        target_obstype : str or None, optional
            The observation type (sensor) to resample. If None, all sensors
            will be resampled to the same frequency. Default is None.
        shift_tolerance : str or pandas.Timedelta, optional
            The maximum allowed time shift tolerance for aligning data during
            resampling. Default is 4 minutes.
        origin_utc : pandas.Timestamp, or None, optional
            The origin timestamp (=first timestamp) of the target records. Can be a pandas.Timestamp or None.
            If a timezone-naive timestamp is given, it is assumed to be in UTC.
            If None, then a suitable origin is used by simplifying the current origing and making
            sure that origin_simplify_tolerance is met. Default is None.
        origin_simplify_tolerance : pandas.Timedelta, optional
            This is only used when origin is None. The tolerance for simplifying the origin
            alignment during resampling. Default is 4 minutes.


        Returns
        -------
        None

        Notes
        -----
        - If `target_obstype` is None, all sensors in `self.sensordata` will be
          resampled to the same frequency.
        - If `target_obstype` is specified, it must be a known observation type.
        Warning
        -------
        Since the gaps depend on the record’s frequency and origin, all gaps
        are removed and re-located. All progress in gap(filling) will be lost.

        Warning
        -------
        Cumulative tolerance errors can be introduced when this method is called multiple times.
        """
        logger.debug("Entering resample for %s", self)
        # format arguments
        target_freq = fmt_timedelta_arg(target_freq)
        shift_tolerance = fmt_timedelta_arg(shift_tolerance)
        origin = fmt_datetime_arg(origin, none_is_none=True)

        if target_obstype is None:
            for sensor in self.sensordata.values():
                sensor.resample(
                    target_freq=target_freq,
                    shift_tolerance=shift_tolerance,
                    origin=origin,
                    origin_simplify_tolerance=origin_simplify_tolerance,
                )
        else:
            # check if target obstype is knonw
            self._obstype_is_known_check(target_obstype)

            # resample
            self.sensordata[target_obstype].resample(
                target_freq=target_freq,
                shift_tolerance=shift_tolerance,
                origin=origin,
                origin_simplify_tolerance=origin_simplify_tolerance,
            )

    def convert_outliers_to_gaps(
        self, all_observations: bool = False, obstype: str = "temp"
    ) -> None:
        """
        Convert outlier values in the observation data to gaps.

        This method replaces outlier values with gaps, effectively marking them as missing.
        In practice, this method is used so that generated gaps can be filled, to obtain a continious timeseries.

        The operation can be applied to all observation types or a specific observation type.

        Parameters
        ----------
        all_observations : bool, optional
            If True, convert outliers to gaps for all observation types.
            If False, only convert outliers specified by the obstype argument. The default is False.
        obstype : str, optional
            The type of observation to convert outliers for. This parameter is
            only used if `all_observations` is False. The default is "temp".

        Returns
        -------
        None

        Warning
        -------
        QC labels are lost when outliers are converted to gaps.
        """
        logger.debug("Entering convert_outliers_to_gaps for %s", self)
        if all_observations:
            for sensor in self.sensordata.values():
                sensor.convert_outliers_to_gaps()

        else:
            self._obstype_is_known_check(obstype)
            self.get_sensor(obstype).convert_outliers_to_gaps()

    def _rename(self, targetname):
        # Note: Not for users, one could accidentaly rename to another station in the dataset.
        # So --> only accecible as method in the dataset, that checks this possible error.

        # rename all
        self._name = str(targetname)
        self._site._stationname = str(targetname)
        for sensordat in self.sensordata.values():
            sensordat._rename(targetname)

    def get_static_gee_point_data(
        self,
        geestaticdatasetmanager: GEEStaticDatasetManager,
        overwrite: bool = True,
        initialize_gee: bool = True,
    ):
        """Extract static data from GEE dataset at Station locations.

        Retrieve Google Earth Engine (GEE) point data of a static dataset, for the station.
        The retrieved data is stored if overwrite is True.

        Parameters
        ----------
        geestaticdatasetmanager : GEEStaticDatasetManager
            An instance of `GEEStaticDatasetManager` representing the static GEE dataset to query.
        overwrite : bool, optional
            If True, the retrieved data will overwrite existing data in the Station's metadata.
            Default is True.
        initialize_gee : bool, optional
            If True, initializes the GEE API before fetching the data. Default is True.

        Returns
        -------
        str or float
            The retrieved metadata from the specified GEE dataset.

        Notes
        -----
        This method interacts with the GEE API to fetch metadata for the station's location.
        Ensure that the GEE API is properly authenticated and initialized before using this method.
        """
        logger.debug("Entering get_static_gee_point_data for %s", self)
        if not isinstance(geestaticdatasetmanager, GEEStaticDatasetManager):
            raise ValueError(
                f"geestaticdataset should be an isntance of GeeStaticDataset, not {type(geestaticdatasetmanager)}"
            )

        value = self.site.get_gee_point_metadata(
            geestaticdataset=geestaticdatasetmanager, initialize_gee=initialize_gee
        )
        if overwrite:
            self.site.set_geedata(geestaticdatasetmanager.name, value)

        return value

    def get_LCZ(self, overwrite: bool = True, initialize_gee: bool = True) -> str:
        """
        Retrieve Local Climate Zone (LCZ) for the stations using Google Earth Engine (GEE).

        Parameters
        ----------
        overwrite : bool, optional
            If True, overwrite existing LCZ data if stored in the Site attribute. Default is True.
        initialize_gee : bool, optional
            If True, initialize the Google Earth Engine API before fetching data. Default is True.

        Returns
        -------
        str
            The LCZ of the station.

        Notes
        -----
        This method relies on the `get_static_gee_point_data` function and the
        `default_gee_datasets` dictionary to fetch the LCZ data.
        """
        logger.debug("Entering get_LCZ for %s", self)
        return self.get_static_gee_point_data(
            geestaticdatasetmanager=default_gee_datasets["LCZ"],
            overwrite=overwrite,
            initialize_gee=initialize_gee,
        )

    def get_altitude(
        self, overwrite: bool = True, initialize_gee: bool = True
    ) -> float:
        """
        Retrieve altitude for the stations using Google Earth Engine (GEE).

        Parameters
        ----------
        overwrite : bool, optional
            If True, overwrite existing altitude data if stored in the Site attribute. Default is True.
        initialize_gee : bool, optional
            If True, initialize the Google Earth Engine API before fetching data. Default is True.

        Returns
        -------
        float
            The altitude of the station.

        Notes
        -----
        This method relies on the `get_static_gee_point_data` function and the
        `default_gee_datasets` dictionary to fetch the altitude data.
        """
        logger.debug("Entering get_altitude for %s", self)
        return self.get_static_gee_point_data(
            geestaticdatasetmanager=default_gee_datasets["altitude"],
            overwrite=overwrite,
            initialize_gee=initialize_gee,
        )

    def get_static_gee_buffer_fraction_data(
        self,
        geestaticdataset: GEEStaticDatasetManager,
        buffers: list = [100],
        aggregate: bool = False,
        overwrite: bool = True,
        initialize_gee: bool = True,
    ) -> dict:
        """
        Extract circular buffer fractions of a GEE dataset at Station locations.

        Parameters
        ----------
        geestaticdataset : GEEStaticDatasetManager
            An instance of GEEStaticDatasetManager used to retrieve static GEE data.
        buffers : list, optional
            A list of buffer radii (in meters) for which to compute the buffer fractions.
            Default is [100].
        aggregate : bool, optional
            If True, aggregate the buffer fraction data. Default is False.
        overwrite : bool, optional
            If True, overwrite the existing buffer fraction data in the station's Site attribute.
            Default is True.
        initialize_gee : bool, optional
            If True, initialize the GEE environment before retrieving data. Default is True.

        Returns
        -------
        dict
            A nested dictionary where the keys are buffer radii and the values are the
            corresponding buffer fraction data.

        Warning
        -------
        This method makes use of GEE API. Make sure that you have access and user rights to use the GEE API.

        Warning
        -------
        It can happen that for stations located on small islands, or close to the coast, the sea-mask is not used as a landcover fraction.
        """
        logger.debug("Entering get_static_gee_buffer_fraction_data for %s", self)
        if not isinstance(geestaticdataset, GEEStaticDatasetManager):
            raise ValueError(
                f"geestaticdataset should be an isntance of GeeStaticDataset, not {type(geestaticdataset)}"
            )

        nesteddict = self.site.get_gee_point_buffer_fractions(
            geestaticdataset=geestaticdataset,
            buffers=buffers,
            aggregate=aggregate,
            initialize_gee=initialize_gee,
        )
        if overwrite:
            for bufferrad, fractions in nesteddict.items():
                self.site.set_gee_buffered_frac_data(buffer=bufferrad, data=fractions)
        return nesteddict

    def get_landcover_fractions(
        self, buffers: list = [100], aggregate: bool = False, overwrite: bool = True
    ) -> dict:
        """
        Get landcover fractions for a circular buffer at the station using GEE.

        Wrapper method for `get_static_gee_buffer_fraction_data` to retrieve land cover fractions
        based on the ESA worldcoverV200 dataset.
        Parameters
        ----------
        buffers : list of int, optional
            List of buffer sizes (in meters) to calculate land cover fractions for.
            Default is [100].
        aggregate : bool, optional
            If True, aggregates the data over the buffers. Default is False.
        overwrite : bool, optional
            If True, overwrites existing data. Default is True.

        Returns
        -------
        dict
            A nested dictionary where the keys are buffer radii and the values are the
            corresponding (aggregated) landcover classes.

        Warning
        -------
        This method makes use of GEE API. Make sure that you have access and user rights to use the GEE API.
        """
        logger.debug("Entering get_landcover_fractions for %s", self)
        return self.get_static_gee_buffer_fraction_data(
            geestaticdataset=default_gee_datasets["worldcover"],
            buffers=buffers,
            aggregate=aggregate,
            overwrite=overwrite,
        )

    def get_gee_timeseries_data(
        self,
        geedynamicdatasetmanager: GEEDynamicDatasetManager,
        startdt_utc: Union[datetime, pd.Timestamp, str, None] = None,
        enddt_utc: Union[datetime, pd.Timestamp, str, None] = None,
        target_obstypes: list = ["temp"],
        get_all_bands: bool = False,
        drive_filename: str = None,
        drive_folder: str = "gee_timeseries_data",
        force_direct_transfer: bool = False,
        force_to_drive: bool = False,
    ) -> Union[pd.DataFrame, None]:
        """Extract time series data from GEE.

        Extracts time series (extraction at station locations) data from a Google Earth Engine (GEE) dynamic dataset
        for a specified time range and observation types.

        If the data request is small, GEE sends the data directly. If not, the data will be
        writen to a CSV file and saved on your Google Drive. In this case, you can import the modeldata
        by using the ´Datset.import_gee_data_from_file()´ method.


        Parameters
        ----------
        geedynamicdatasetmanager : GEEDynamicDatasetManager
            The dynamic dataset manager instance describing the target GEE dataset.
        startdt_utc : datetime or str, optional
            The start datetime in UTC for the time series data. If None, the
            station's start datetime (converted to UTC) is used.
        enddt_utc : datetime or str, optional
            The end datetime in UTC for the time series data. If None, the
            station's end datetime (converted to UTC) is used.
        target_obstypes : list of str, optional
            List of observation types to extract. Defaults to ["temp"].
        get_all_bands : bool, optional
            If True, extracts all bands from the dataset. Defaults to False.
        drive_filename : str, optional
            The filename of the CSV file to use when saving the data to Google Drive. If None, a
            default name is generated.
        drive_folder : str, optional
            The folder name in Google Drive where the file will be saved. Defaults
            to "gee_timeseries_data".
        force_direct_transfer : bool, optional
            If True, forces direct data transfer instead of using Google Drive.
            Defaults to False.
        force_to_drive : bool, optional
            If True, forces saving the data to Google Drive. Defaults to False.

        Returns
        -------
        pandas.DataFrame or None
            A DataFrame containing the extracted time series data. Returns None if
            no data is retrieved.

        Note
        -----
        If a timezone unaware datetime is given as an argument, it is interpreted as if it has the same timezone as the observations.

        Warning
        -------
        This method makes use of GEE API. Make sure that you have access and user rights to use the GEE API.

        Warning
        -------
        When extracting large amounts of data,
        the timeseries data will be written to a file and saved
        on your Google Drive. In this case, you can import the modeldata
        by using the Dataset.import_gee_data_from_file() method.

        Notes
        -----
        - The method creates `ModelTimeSeries` instances for each valid observation
          type in the extracted data and appends them to the station's model data.
        - If no data is returned by the GEE API request, a warning is logged and
          the method returns None.
        """
        logger.debug("Entering get_gee_timeseries_data for %s", self)
        # Check geedynamic dataset
        if not isinstance(geedynamicdatasetmanager, GEEDynamicDatasetManager):
            raise ValueError(
                f"geedynamicdataset should be an isntance of GeeDynamicDataset, not {type(geedynamicdatasetmanager)}"
            )

        # Format datetime arguments
        if startdt_utc is None:
            startdt_utc = self.start_datetime.tz_convert("UTC")
        else:
            startdt_utc = fmt_datetime_arg(startdt_utc, tz_if_dt_is_naive="UTC")

        if enddt_utc is None:
            enddt_utc = self.end_datetime.tz_convert("UTC")
        else:
            enddt_utc = fmt_datetime_arg(enddt_utc, tz_if_dt_is_naive="UTC")

        # check if target_obstypes are mapped to bands
        for obst in target_obstypes:
            if obst not in geedynamicdatasetmanager.modelobstypes.keys():
                raise MetObsMetadataNotFound(
                    f"{obst} is not a known modelobstype of {geedynamicdatasetmanager}."
                )

        # create specific name for the file that might be written to Drive
        if drive_filename is None:
            drive_filename = (
                f"{geedynamicdatasetmanager.name}_timeseries_data_of_{self.name}.csv"
            )

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

        # Create ModelTimeSeries instances
        for modelobscol in df.columns:
            if modelobscol in geedynamicdatasetmanager.modelobstypes.keys():
                modeltimeseries = ModelTimeSeries(
                    site=self.site,
                    datarecords=df[modelobscol].to_numpy(),
                    timestamps=df.index.get_level_values("datetime").to_numpy(),
                    obstype=geedynamicdatasetmanager.modelobstypes[modelobscol],
                    datadtype=np.float32,
                    timezone="UTC",
                    modelname=geedynamicdatasetmanager.name,
                    modelvariable=geedynamicdatasetmanager.modelobstypes[
                        modelobscol
                    ].model_band,
                )
                # todo: duplicacy check
                self._modeldata[modeltimeseries.obstype.name] = modeltimeseries
            else:
                logger.info(
                    f"Skip {modelobscol} for creating a ModelTimeeries because of unknown obstype."
                )

        return df

    def gross_value_check(
        self,
        target_obstype: str = "temp",
        lower_threshold: float = -15.0,
        upper_threshold: float = 39.0,
    ) -> None:
        """
        Identify outliers based on thresholds.

        Parameters
        ----------
        target_obstype : str, optional
            The target observation to check. by default "temp"
        lower_threshold : float, optional
           Thresholds to flag records below as outliers. The default is -15.0.
        upper_threshold : float, optional
            Thresholds to flag records above as outliers. The default is 39.0.

        Returns
        -------
        None

        Note
        -----
        This method modifies the outliers in place and does not return anything.
        You can use the `outliersdf` property to view all flagged outliers.
        """
        logger.debug("Entering gross_value_check for %s", self)
        # argument validity checks
        self._obstype_is_known_check(target_obstype)

        self.get_sensor(target_obstype).gross_value_check(
            lower_threshold=lower_threshold, upper_threshold=upper_threshold
        )

    def persistence_check(
        self,
        target_obstype: str = "temp",
        timewindow: Union[str, pd.Timedelta] = pd.Timedelta("60min"),
        min_records_per_window: int = 5,
    ) -> None:
        """Check if values are not constant in a moving time window.

        Perform a persistence check on a time series to identify periods where observations remain constant
        within a specified time window. If the values are constant, all records in the moving window are
        flagged as outliers.
        Parameters
        ----------
        target_obstype : str, optional
            The target observation to check. by default "temp"
        timewindow : str or pd.Timedelta
            The size of the rolling time window to check for persistence.
            The default is pd.Timedelta("60min")
        min_records_per_window : int
            The minimum number of non-NaN records required within the time window for the check to be valid.
            The default is 5

        Returns
        -------
        None

        Notes
        -----
        - This method modifies the outliers in place and does not return anything.
          You can use the `outliersdf` property to view all flagged outliers.
        - If the minimum number of records per window is locally not met, the function logs a warning and skips
          the persistence check.
        - This function can be computationally expensive for large datasets or small time windows.
        - The repetitions check is similar to the persistence check, but not identical.
          The persistence check uses thresholds that are meteorologically based (i.e. the moving window is defined by a duration),
          in contrast to the repetitions check whose thresholds are instrumentally based (i.e. the "window" is defined by a number of records.)

        Warnings
        -------
        If the minimum number of records per window is not met over the full time series, a warning is logged, and the function
        returns an empty DatetimeIndex.
        """
        logger.debug("Entering persistence_check for %s", self)
        # argument checks
        self._obstype_is_known_check(target_obstype)
        timewindow = fmt_timedelta_arg(timewindow)

        # apply check on the sensordata
        self.get_sensor(target_obstype).persistence_check(
            timewindow=timewindow, min_records_per_window=min_records_per_window
        )

    def repetitions_check(
        self, target_obstype: str = "temp", max_N_repetitions: int = 5
    ) -> None:
        """Test if an observation changes after a number of repetitions.

        Perform a check that tests if the observation changes after a number of repetitions.
        If a value is repeated more than the specified number of times, all the repeated
        records are flagged as outliers.

        Be aware that the performance of this check depends on the `max_N_repetitions`
        and the time resolution of the observations.
        Parameters
        ----------
        target_obstype : str, optional
            The target observation to check. by default "temp"
        max_N_repetitions : int
            The maximum number of repetitions allowed before the records are flagged as outliers.
            If the number of repetitions exceeds this value, all repeated records are flagged as outliers. The default is 5.

        Returns
        -------
        None

        Notes
        -----
        - This method modifies the outliers in place and does not return anything.
          You can use the `outliersdf` property to view all flagged outliers.
        - The repetitions check is similar to the persistence check, but not identical.
          The persistence check uses thresholds that are meteorologically based (i.e. the moving window is defined by a duration),
          in contrast to the repetitions check whose thresholds are instrumentally based (i.e. the "window" is defined by a number of records.)

        Warnings
        -------
        If the minimum number of records per window is not met over the full time series, a warning is logged, and the function
        returns an empty DatetimeIndex.
        """
        logger.debug("Entering repetitions_check for %s", self)
        # argument checks
        self._obstype_is_known_check(target_obstype)

        # apply check on the sensordata
        self.get_sensor(target_obstype).repetitions_check(
            max_N_repetitions=max_N_repetitions
        )

    def step_check(
        self,
        target_obstype: str = "temp",
        max_increase_per_second: Union[int, float] = 8.0 / 3600.0,
        max_decrease_per_second: Union[int, float] = -10.0 / 3600.0,
    ) -> None:
        """Check for 'spikes' and 'dips' in a timeseries.

        Test if observations do not produce spikes in timeseries. The maximum
        allowed increase and decrease per second is set in the argument,
        and is tested to each record (with respect to the previous record).

        If the difference between two consecutive records (i.e., the spike/dip) is larger than the
        threshold, the record is flagged as an outlier.
        Parameters
        ----------
        target_obstype : str, optional
            The target observation to check. by default "temp"
        max_increase_per_second : int or float, >0, optional
            The maximum allowed increase (per second). This value is extrapolated to the time resolution of records.
            This value must be positive!  The default is 8.0/3600.0
        max_decrease_per_second : int or float, <0, optional
            The maximum allowed decrease (per second). This value is extrapolated to the time resolution of records.
            This value must be negative!, The default is -10.0/3600.0

        Returns
        -------
        None

        Notes
        -----
        - This method modifies the outliers in place and does not return anything.
          You can use the `outliersdf` property to view all flagged outliers.
        - In general, for temperatures, the decrease threshold is set less stringent than the increase
          threshold. This is because a temperature drop is meteorologically more
          common than a sudden increase which is often the result of a radiation error.
        """
        logger.debug("Entering step_check for %s", self)
        # argument checks
        self._obstype_is_known_check(target_obstype)

        # apply check on the sensordata
        self.get_sensor(target_obstype).step_check(
            max_increase_per_second=max_increase_per_second,
            max_decrease_per_second=max_decrease_per_second,
        )

    def window_variation_check(
        self,
        target_obstype: str = "temp",
        timewindow: pd.Timedelta = pd.Timedelta("1h"),
        min_records_per_window: int = 3,
        max_increase_per_second: Union[int, float] = 8.0 / 3600,
        max_decrease_per_second: Union[int, float] = -10.0 / 3600,
    ) -> None:
        """
        Test if the increase/decrease in a timewindow exceeds a threshold.

        Checks if the variation of observations in time,
        does not exceed a threshold. This is done by applying a moving window
        over the time series. The moving window is defined by a duration (timewindow),
        and tested if the window contains at least a minimum number of records.

        If the observations in the window increase/decrease more than a threshold, all
        observations in the window are flagged as outliers. The threshold is defined by the
        maximum increase/decrease per second multiplied by the window size in seconds.
        Parameters
        ----------
        target_obstype : str, optional
            The target observation to check. by default "temp"
        timewindow : pd.Timedelta
            The duration of the moving window. This should be a pandas Timedelta object. The default is pd.Timedelta("1h")
        min_records_per_window : int
            The minimum number of non-NaN records required within the time window for the check to be valid.
            This is dependent on the time resolution of the records. The default is 3
        max_increase_per_second : int or float, >0
            The maximum allowed increase (per second). This value is extrapolated to the window duration.
            This value must be positive! The default is 8.0/3600
        max_decrease_per_second : int or float, <0
            The maximum allowed decrease (per second). This value is extrapolated to the window duration.
            This value must be negative! The default is -10.0/3600

        Returns
        -------
        None

        Notes
        -----
        - This method modifies the outliers in place and does not return anything.
          You can use the `outliersdf` property to view all flagged outliers.
        - In general, for temperatures, the decrease threshold is set less stringent than the increase
          threshold. This is because a temperature drop is meteorologically more
          common than a sudden increase which is often the result of a radiation error.
        - A suitable value for the min_records_per_window depends on the time resolution of the records and the window size.
        - This check is similar to the step check, but not identical. The step check a maximum allowed increase/decrease
          with respect to the previous value. The window variation check uses a moving window to test the maximum allowed variation.
        """
        logger.debug("Entering window_variation_check for %s", self)
        # argument checks
        self._obstype_is_known_check(target_obstype)
        timewindow = fmt_timedelta_arg(timewindow)

        # apply check on the sensordata
        self.get_sensor(target_obstype).window_variation_check(
            timewindow=timewindow,
            min_records_per_window=min_records_per_window,
            max_increase_per_second=max_increase_per_second,
            max_decrease_per_second=max_decrease_per_second,
        )

    def get_qc_stats(
        self, target_obstype: str = "temp", make_plot: bool = True
    ) -> pd.DataFrame:
        """Generate quality control (QC) frequency statistics.

        This method calculates the frequency statistics for various QC checks
        applied, and other labels. The order of checks is taken into
        account.

        Frequency of labels is combuted based on the set of all labels (for all
        records including gaps). The effectiveness of a check is shown by
        the frequency of outliers wrt the numer of records that were given
        to the check (thus taking into account the order of checks).

        The frequencies are retured in a dataframe, and can be plotted
        as pie charts.

        Parameters
        ----------
        target_obstype : str, optional
            The target observationtype for which to compute frequency statistics, by default "temp".
        make_plot : bool, optional
            If True, a figure with pie charts representing the frequencies is generated. The default is True.

        Returns
        -------
        pandas.DataFrame
            A DataFrame containing the QC frequency statistics. The DataFrame
            has a multi-index with the station name and QC check label, and
            includes the following columns:
            - `N_all`: Total number of records in the dataset (including gaps).
            - `N_labeled`: Number of records with the specific label.
            - `N_checked`: Number of records checked for the specific QC check.
        """
        logger.debug("Entering get_qc_stats for %s", self)
        # argument checks
        self._obstype_is_known_check(target_obstype)
        # get freq statistics
        qc_df = self.get_sensor(target_obstype).get_qc_freq_statistics()

        if make_plot:
            plotdf = qc_df.reset_index().drop(columns=["name"]).set_index("qc_check")

            fig = plotting.qc_overview_pies(df=plotdf)
            fig.suptitle(
                f"QC frequency statistics of {target_obstype} on Station level: {self.stationname}."
            )
            return fig
        else:
            return qc_df

    def _set_stationdata(self, all_sensor_data: list) -> dict:
        """Set the station's sensor data dictionary."""
        return {
            sensor_data.obstype.name: sensor_data for sensor_data in all_sensor_data
        }

    def _set_tz(self, current_tz):
        """Set the timezone (not implemented)."""
        pass

    def make_plot_of_modeldata(
        self,
        obstype: str = "temp",
        linecolor: Union[str, None] = None,
        title: Union[str, None] = None,
        linestyle: str = "--",
        ax: Union[None, Axes] = None,
        figkwargs: dict = {},
    ) -> Axes:
        """
        Generate a timeseries plot of model data for a specific observation type.

        Parameters
        ----------
        obstype : str, optional
            The type of observation to plot modeldata for, by default "temp".
        linecolor : str or None, optional
            The color of the line in the plot. If None, a default color map is used.
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
        logger.debug("Entering make_plot_of_modeldata for %s", self)
        # test if the obstype has modeldata
        self._obstype_has_modeldata_check(obstype)

        # Create new axes if needed
        if ax is None:
            ax = plotting.create_axes(**figkwargs)

        plotdf = (
            self.modeldatadf.xs(obstype, level="obstype", drop_level=False)
            .assign(name=self.name)
            .reset_index()
            .set_index(["name", "obstype", "datetime"])
            .sort_index()
        )

        plotdf = plotdf[["value"]]
        plotdf["label"] = label_def["goodrecord"][
            "label"
        ]  # Just so that they are plotted as lines

        # Define linecolor (needed here if modeldata is added )
        if linecolor is None:
            colormap = plotting.create_categorical_color_map([self.name])
        else:
            colormap = {self.name: linecolor}
        ax = plotting.plot_timeseries_color_by_station(
            plotdf=plotdf,
            colormap=colormap,
            show_outliers=False,  # will not be used,
            show_gaps=False,  # will not be used
            ax=ax,
            linestyle=linestyle,
            legend_prefix=f"{self.modeldata[obstype].modelname}:{self.modeldata[obstype].modelvariable}@",
        )
        # Styling
        obstypeinstance = self.modeldata[obstype].obstype

        # Set title:
        if title is None:
            plotting.set_title(
                ax, f"{obstypeinstance.name} data for station {self.name}"
            )
        else:
            plotting.set_title(ax, title)

        # Set ylabel
        plotting.set_ylabel(ax, obstypeinstance._get_plot_y_label())

        # Set xlabel
        cur_tz = plotdf.index.get_level_values("datetime").tz
        plotting.set_xlabel(ax, f"Timestamps (in {cur_tz})")

        # Format timestamp ticks
        plotting.format_datetime_axes(ax)

        # Add legend
        plotting.set_legend(ax)

        return ax

    def make_plot(
        self,
        obstype: str = "temp",
        colorby: Literal["station", "label"] = "label",
        show_modeldata: bool = False,
        linecolor: Union[str, None] = None,
        show_outliers=True,
        show_gaps=True,
        title: Union[str, None] = None,
        ax: Union[Axes, None] = None,
        figkwargs: dict = {},
    ) -> Axes:
        """
        Generate a time series plot for observational data.

        Parameters
        ----------
        obstype : str, optional
            The type of observation to plot. Default is "temp".
        colorby : {"station", "label"}, optional
            Determines how the data is colored in the plot.
            - "station": Colors by station.
            - "label": Colors by label (the labels refer to the status of a record).
            Default is "label".
        show_modeldata : bool, optional
            If True, includes model data (of the same obstype) if present, in the plot. Default is False.
        linecolor : str or None, optional
            The color of the line for the model data. If None, a default categorical color map is used. Default is None.
        show_outliers : bool, optional
            If True, includes outliers in the plot. Default is True.
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

        Notes
        -----
        - The method checks if the specified `obstype` is known before proceeding.
        - The plot can include observational data, model data, or both.
        - The x-axis timestamps are formatted according to the timezone of the data.
        """
        logger.debug("Entering make_plot for %s", self)
        # test if obstype have sensordata
        self._obstype_is_known_check(obstype)

        # Create new axes if needed
        if ax is None:
            ax = plotting.create_axes(**figkwargs)

        if show_modeldata:
            if linecolor is None:
                colormap = plotting.create_categorical_color_map([self.name])
            else:
                colormap = {self.name: linecolor}

            ax = self.make_plot_of_modeldata(
                obstype=obstype,
                linecolor=linecolor,
                ax=ax,
                figkwargs=figkwargs,
                title=title,
            )

        # Create plotdf
        plotdf = (
            self.df.xs(obstype, level="obstype", drop_level=False)
            .assign(name=self.name)
            .reset_index()
            .set_index(["name", "obstype", "datetime"])
            .sort_index()
        )

        if colorby == "station":
            # Define linecolor (needed here if modeldata is added )
            if linecolor is None:
                colormap = plotting.create_categorical_color_map([self.name])
            else:
                colormap = {self.name: linecolor}
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

        # Styling

        obstypeinstance = self.get_sensor(obstype).obstype

        # Set title:
        if title is None:
            plotting.set_title(
                ax, f"{obstypeinstance.name} data for station {self.name}"
            )
        else:
            plotting.set_title(ax, title)

        # Set ylabel
        plotting.set_ylabel(ax, obstypeinstance._get_plot_y_label())

        # Set xlabel
        cur_tz = plotdf.index.get_level_values("datetime").tz
        plotting.set_xlabel(ax, f"Timestamps (in {cur_tz})")

        # Format timestamp ticks
        plotting.format_datetime_axes(ax)

        # Add legend
        plotting.set_legend(ax)

        return ax

    def fill_gaps_with_raw_modeldata(
        self, target_obstype: str, overwrite_fill: bool = False
    ) -> None:
        """Fill the gap(s) using model data without correction.

        This method fills all the gaps of a specific *target_obstype*, by directly interpolating
        the model data to the missing records.

        Parameters
        ----------
        target_obstype :  str
            The target obstype to fill the gaps for.
        overwrite_fill : bool, optional
            If True, the status of a `gap` and present gapfill info will be ignored and overwritten.
            If False, only gaps without gapfill data are filled. Defaults to False.

        Returns
        -------
        None

        Notes
        -----
        A schematic description of the raw model data gap fill:

        1. Check if the target_obstype is knonw, and if the corresponding modeldata is present.
        2. Iterate over the gaps of the target_obstype.
        3. Check the compatibility of the `ModelTimeSeries` with the `gap`.
        4. Ensure both the `ModelTimeSeries` and `gap` have the same timezone.
        5. Interpolate the model data to match the missing records in the gap.
        6. Update the `gap` attributes with the interpolated values, labels, and details.
        """
        logger.debug("Entering fill_gaps_with_raw_modeldata for %s", self)
        # obstype check
        self._obstype_is_known_check(obstype=target_obstype)

        # Get modeltimeseries
        if target_obstype not in self.modeldata:
            raise MetObsModelDataError(
                f"No Modeldata found for {target_obstype} in {self}"
            )

        modeltimeseries = self.modeldata[target_obstype]

        # fill the gaps
        self.get_sensor(target_obstype).fill_gap_with_modeldata(
            modeltimeseries=modeltimeseries, method="raw", overwrite_fill=overwrite_fill
        )

    def fill_gaps_with_debiased_modeldata(
        self,
        target_obstype: str,
        leading_period_duration: Union[pd.Timedelta, str] = pd.Timedelta("24h"),
        min_leading_records_total: int = 60,
        trailing_period_duration: Union[pd.Timedelta, str] = pd.Timedelta("24h"),
        min_trailing_records_total: int = 60,
        overwrite_fill: bool = False,
    ) -> None:
        """Fill the gaps using modeldata corrected for the bias.

        This method fills the gaps using model data corrected for bias.
        The bias is estimated using a leading (before the gap)
        and trailing (after the gap) period. The bias is computed by combining the
        leading and trailing period, and comparing the model with the observations
        (not labeled as outliers). The model data is then interpolated to the missing
        records, and corrected with the estimated bias.

        Parameters
        ----------
        target_obstype :  str
            The target obstype to fill the gaps for.
        leading_period_duration : str or pd.Timedelta, optional
            The duration of the leading period. The default is "24h".
        min_leading_records_total : int, optional
            The minimum number of records required in the leading period. The default is 60.
        trailing_period_duration : str or pd.Timedelta, optional
            The duration of the trailing period. The default is "24h".
        min_trailing_records_total : int, optional
            The minimum number of records required in the trailing period. The default is 60.
        overwrite_fill : bool, optional
            If True, the status of a `gap` and present gapfill info will be ignored and overwritten.
            If False, only gaps without gapfill data are filled. The default is False.

        Returns
        ----------
        None

        Notes
        -----
        A schematic description of the debiased modeldata gap fill:

        1. Check if the target_obstype is knonw, and if the corresponding modeldata is present.
        2. Iterate over the gaps of the target_obstype.
        3. Check the compatibility of the `ModelTimeSeries` with the `gap`.
        4. Construct a leading and trailing sample, and test if they meet the required conditions.
        5. Compute the bias of the modeldata (combine leading and trailing samples).
        6. Fill the gap records by using raw (interpolated) modeldata that is corrected by subtracting the bias.
        7. Update the `gap` attributes with the interpolated values, labels, and details.

        """
        logger.debug("Entering fill_gaps_with_debiased_modeldata for %s", self)

        # special formatters
        leading_period_duration = fmt_timedelta_arg(leading_period_duration)
        trailing_period_duration = fmt_timedelta_arg(trailing_period_duration)

        # obstype check
        self._obstype_is_known_check(obstype=target_obstype)

        # Get modeltimeseries
        if target_obstype not in self.modeldata:
            raise MetObsModelDataError(
                f"No Modeldata found for {target_obstype} in {self}"
            )

        modeltimeseries = self.modeldata[target_obstype]

        # fill the gaps
        self.get_sensor(target_obstype).fill_gap_with_modeldata(
            modeltimeseries=modeltimeseries,
            method="debiased",
            method_kwargs={
                "leading_period_duration": leading_period_duration,
                "min_leading_records_total": min_leading_records_total,
                "trailing_period_duration": trailing_period_duration,
                "min_trailing_records_total": min_trailing_records_total,
            },
            overwrite_fill=overwrite_fill,
        )

    def fill_gaps_with_diurnal_debiased_modeldata(
        self,
        target_obstype: str,
        leading_period_duration: Union[pd.Timedelta, str] = pd.Timedelta("24h"),
        trailing_period_duration: Union[pd.Timedelta, str] = pd.Timedelta("24h"),
        min_debias_sample_size: int = 6,
        overwrite_fill: bool = False,
    ) -> None:
        """Fill the gaps using modeldata corrected for the diurnal-bias.

        This method fills the gap using model data corrected for its diuranal-bias.
        The diurnal bias is a bias that is estimated for each timestamp in the leading
        and trailing period. All biases are averaged over hour, minute and second, to
        obtain a diurnal bias (for each timestamp).


        Parameters
        ----------
        target_obstype :  str
            The target obstype to fill the gaps for.
        leading_period_duration :  str or pd.Timedelta, optional
            The duration of the leading period. That is the period before the gap, used
            for bias estimation. The default is "24h".
        trailing_period_duration :  str or pd.Timedelta, optional
            The duration of the trailing period. That is the period after the gap, used
            for bias estimation. The default is "24h".
        min_debias_sample_size : int, optional
            The minimum number of samples required for bias estimation. The default is 6.
        overwrite_fill : bool, optional
            If True, the status of a `gap` and present gapfill info will be ignored and overwritten.
            If False, only gaps without gapfill data are filled. The default is False.

        Returns
        -------
        None

        Notes
        -----
        A schematic description of the diurnal debiased modeldata gap fill:

        1. Check if the target_obstype is knonw, and if the corresponding modeldata is present.
        2. Iterate over the gaps of the target_obstype.
        3. Check the compatibility of the `ModelTimeSeries` with the `gap`.
        4. Construct a leading and trailing sample, and test if they meet the required conditions.
          The required conditions are tested by testing the samplesizes per hour, minute and second for the leading + trailing periods.
        5. A diurnal bias is computed by grouping to hour, minute and second, and averaging the biases.
        6. Fill the gap records by using raw (interpolated) modeldata that is corrected by subtracting the coresponding diurnal bias.
        7. Update the `gap` attributes with the interpolated values, labels, and details.

        Notes
        --------
        Note that a suitable `min_debias_sample_size` depends on the sizes of the
        leading- and trailing periods, and also on the time resolution gap (=time resolution of the corresponding SensorData).

        References
        -----------
        Jacobs .A, et. al. (2024) `Filling gaps in urban temperature observations by debiasing ERA5 reanalysis data <https://doi.org/10.1016/j.uclim.2024.102226>`_
        """
        logger.debug("Entering fill_gaps_with_diurnal_debiased_modeldata for %s", self)
        # special formatters
        leading_period_duration = fmt_timedelta_arg(leading_period_duration)
        trailing_period_duration = fmt_timedelta_arg(trailing_period_duration)

        # obstype check
        self._obstype_is_known_check(obstype=target_obstype)

        # Get modeltimeseries
        if target_obstype not in self.modeldata:
            raise MetObsModelDataError(
                f"No Modeldata found for {target_obstype} in {self}"
            )

        modeltimeseries = self.modeldata[target_obstype]

        # fill the gaps
        self.get_sensor(target_obstype).fill_gap_with_modeldata(
            modeltimeseries=modeltimeseries,
            method="diurnal_debiased",
            method_kwargs={
                "leading_period_duration": leading_period_duration,
                "trailing_period_duration": trailing_period_duration,
                "min_debias_sample_size": min_debias_sample_size,
            },
            overwrite_fill=overwrite_fill,
        )

    def fill_gaps_with_weighted_diurnal_debiased_modeldata(
        self,
        target_obstype: str,
        leading_period_duration: Union[pd.Timedelta, str] = pd.Timedelta("24h"),
        trailing_period_duration: Union[pd.Timedelta, str] = pd.Timedelta("24h"),
        min_lead_debias_sample_size: int = 2,
        min_trail_debias_sample_size: int = 2,
        overwrite_fill=False,
    ):
        """Fill the gaps using a weighted sum of modeldata corrected for the diurnal-bias and weights wrt the start of the gap.

        This method fills the gaps using model data corrected for its diuranal-bias.
        The diurnal bias is a bias that is estimated for each timestamp in the leading
        and trailing period (seperatly). For both periods seperatly, all biases are averaged over hour, minute and second, to
        obtain a diurnal bias (for each timestamp).

        In addition, a normalized weight is computed for each gap-record indicating the distance (in time) to
        the start and end of the gap. The correction applied on the interpolated (in time) modeldata, is
        thus a weighted sum of corrections comming from both the leading and trailing period.


        Parameters
        ----------
        target_obstype :  str
            The target obstype to fill the gaps for.
        leading_period_duration : str or pd.Timedelta, optional
            The duration of the leading period. That is the period before the gap, used
            for bias estimation. The default is "24h".
        trailing_period_duration : str or pd.Timedelta, optional
            The duration of the trailing period. That is the period after the gap, used
            for bias estimation. The default is "24h".
        min_lead_debias_sample_size : int, optional
            The minimum number of leading-samples required for bias estimation. If this condition is not met, the gap
            is not filled. The default is 2.
        min_trail_debias_sample_size : int, optional
            The minimum number of trailing-samples required for bias estimation. If this condition is not met, the gap
            is not filled. The default is 2.
        overwrite_fill : bool, optional
            If True, the status of a `gap` and present gapfill info will be ignored and overwritten.
            If False, only gaps without gapfill data are filled. The default is False.

        Returns
        --------
        None.


        Notes
        -----
        A schematic description of the weighted diurnal debiased modeldata gap fill:

        1. Check if the target_obstype is knonw, and if the corresponding modeldata is present.
        2. Iterate over the gaps of the target_obstype.
        3. Check the compatibility of the `ModelTimeSeries` with the `gap`.
        4. Construct a leading and trailing sample, and test if they meet the required conditions. The required conditions are tested by testing the samplesizes per hour, minute and second for the leading and trailing periods (seperatly).
        5. A leading and trailing set of diurnal biases are computed by grouping to hour, minute and second, and averaging the biases.
        6. A weight is computed for each gap record, that is the normalized distance to the start and end of the gap.
        7. Fill the gap records by using raw (interpolated) modeldata is corrected by a weighted sum the coresponding diurnal bias for the lead and trail periods.
        8. Update the `gap` attributes with the interpolated values, labels, and details.

        Notes
        --------
        Note that a suitable `min_debias_sample_size` depends on the sizes of the
        leading- and trailing periods, and also on the time resolution gap (=time resolution of the corresponding SensorData).


        References
        -----------
        Jacobs .A, et. al. (2024) `Filling gaps in urban temperature observations by debiasing ERA5 reanalysis data <https://doi.org/10.1016/j.uclim.2024.102226>`_

        """
        logger.debug(
            "Entering fill_gaps_with_weighted_diurnal_debiased_modeldata for %s", self
        )
        # special formatters
        leading_period_duration = fmt_timedelta_arg(leading_period_duration)
        trailing_period_duration = fmt_timedelta_arg(trailing_period_duration)

        # obstype check
        self._obstype_is_known_check(obstype=target_obstype)

        # Get modeltimeseries
        if target_obstype not in self.modeldata:
            raise MetObsModelDataError(
                f"No Modeldata found for {target_obstype} in {self}"
            )

        modeltimeseries = self.modeldata[target_obstype]

        # fill the gaps
        self.get_sensor(target_obstype).fill_gap_with_modeldata(
            modeltimeseries=modeltimeseries,
            method="weighted_diurnal_debiased",
            method_kwargs={
                "leading_period_duration": leading_period_duration,
                "trailing_period_duration": trailing_period_duration,
                "min_lead_debias_sample_size": min_lead_debias_sample_size,
                "min_trail_debias_sample_size": min_trail_debias_sample_size,
            },
            overwrite_fill=overwrite_fill,
        )

    def interpolate_gaps(
        self,
        target_obstype: str,
        method: str = "time",
        max_consec_fill: int = 10,
        n_leading_anchors: int = 1,
        n_trailing_anchors: int = 1,
        max_lead_to_gap_distance: Union[pd.Timedelta, None] = None,
        max_trail_to_gap_distance: Union[pd.Timedelta, None] = None,
        overwrite_fill=False,
        method_kwargs={},
    ):
        """Fill the gap(s) using interpolation of SensorData.

        This method fills all the gaps of a specific *target_obstype*, by directly interpolating
        corresponding SensorData. Each gap is interpolated using the leading and trailing periods of the gap. One can select different
        interpolation methods. By using restrictions on the leading and trailing periods, one can
        ensure that the interpolation is only done when there are enough leading and trailing data available.

        Parameters
        ----------
        target_obstype :  str
            The target obstype to fill the gaps for.
         method:  str, optional
            Interpolation technique to use. See pandas.DataFrame.interpolate
            method argument for possible values. Make sure that
            `n_leading_anchors`, `n_trailing_anchors` and `method_kwargs` are
            set accordingly to the method (higher order interpolation techniques require more leading and trailing anchors). The default is "time".. Defaults to "time".
        max_consec_fill: int, optional
            The maximum number of consecutive timestamps to fill. The result is
            dependent on the time-resolution of the gap (=equal to that of the related SensorData). Defaults to 10.
        n_leading_anchors : int, optional
            The number of leading anchors to use for the interpolation. A leading anchor is
            a near record (not rejected by QC) just before the start of the gap, that is used for interpolation.
            Higher-order interpolation techniques require multiple leading anchors. Defaults to 1.
        n_trailing_anchors : int, optional
            The number of trailing anchors to use for the interpolation. A trailing anchor is
            a near record (not rejected by QC) just after the end of the gap, that is used for interpolation.
            Higher-order interpolation techniques require multiple leading anchors. Defaults to 1.
        max_lead_to_gap_distance:  pd.Timedelta or None, optional
            The maximum time difference between the start of the gap and a
            leading anchor(s). If None, no time restriction is applied on the leading anchors. The default is None.
        max_trail_to_gap_distance : pd.Timedelta or None, optional
            The maximum time difference between the end of the gap and a
            trailing anchor(s). If None, no time restriction is applied on the trailing anchors. Defaults to None.
        overwrite_fill : bool, optional
            If True, the status of a `gap` and present gapfill info will be ignored and overwritten.
            If False, only gaps without gapfill data are filled. Defaults to False.
        method_kwargs : dict, optional
            Extra arguments that are passed to pandas.Dataframe.interpolate() structured in a dict. Defaults to {}.

        Returns
        -------
        None

        Notes
        -----
        A schematic description:

        1. Iterate over all gaps related to the target obstype.
        2. Get the leading and trailing periods of the gap.
        3. Check if the leading and trailing periods are valid.
        4. Create a combined DataFrame with the leading, trailing, and gap data.
        5. Interpolate the missing records using the specified method.
        6. Update the gap attributes with the interpolated values, labels, and details.

        Note
        -------
        The impact of `max_consec_fill` is highly dependent on the resolution
        of your records.

        Note
        ------
        If you want to use a higher-order method of interpolation, make sure to
        increase the `n_leading_anchors` and `n_trailing_anchors` accordingly.
        For example, for a cubic interpolation, you need at least 2 leading and 2 trailing anchors.

        """
        logger.debug("Entering interpolate_gaps for %s", self)
        # special formatters
        max_lead_to_gap_distance = fmt_timedelta_arg(max_lead_to_gap_distance)
        max_trail_to_gap_distance = fmt_timedelta_arg(max_trail_to_gap_distance)
        # obstype check
        self._obstype_is_known_check(obstype=target_obstype)

        # interpolate all the gaps
        self.get_sensor(target_obstype).interpolate_gaps(
            overwrite_fill=overwrite_fill,
            method=method,
            max_consec_fill=max_consec_fill,
            n_leading_anchors=n_leading_anchors,
            n_trailing_anchors=n_trailing_anchors,
            max_lead_to_gap_distance=max_lead_to_gap_distance,
            max_trail_to_gap_distance=max_trail_to_gap_distance,
            method_kwargs=method_kwargs,
        )

    def _obstype_is_known_check(self, obstype: str) -> None:
        """Raise error if obstype is not present in obsdata."""
        if obstype not in self.obsdata.keys():
            raise MetObsSensorDataNotFound(
                f"{self} does not hold {obstype} sensordata. The present sensordata is: {list(self.obsdata.keys())}"
            )

    def _obstype_has_modeldata_check(self, obstype: str) -> None:
        """Raise error if obstype is not present in modeldata."""
        if obstype not in self.modeldata.keys():
            raise MetObsObstypeNotFound(
                f"There is no {obstype} - modeldata present for {self}"
            )
