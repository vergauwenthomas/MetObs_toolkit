import logging
import copy
import numpy as np
import pandas as pd
from typing import Literal, Union
from datetime import datetime
from matplotlib.pyplot import Axes
from pathlib import Path

from metobs_toolkit.site import Site
from metobs_toolkit.backend_collection.argumentcheckers import (
    fmt_datetime_arg,
    fmt_timedelta_arg,
)
from metobs_toolkit.backend_collection.uniqueness import join_collections
from metobs_toolkit.backend_collection.dev_collection import copy_doc
import metobs_toolkit.plot_collection as plotting
from metobs_toolkit.xrconversions import station_to_xr
from metobs_toolkit.backend_collection.errorclasses import (
    MetObsDataAlreadyPresent,
    MetObsMetadataNotFound,
    MetObsModelDataError,
    MetObsSensorDataNotFound,
    MetObsObstypeNotFound,
    MetObsAdditionError,
)
import metobs_toolkit.backend_collection.printing_collection as printing
from metobs_toolkit.backend_collection.df_helpers import save_concat
from metobs_toolkit.settings_collection import label_def
from metobs_toolkit.geedatasetmanagers import (
    GEEStaticDatasetManager,
    GEEDynamicDatasetManager,
)
from metobs_toolkit.geedatasetmanagers import default_datasets as default_gee_datasets
from metobs_toolkit.sensordata import SensorData
from metobs_toolkit.modeltimeseries import ModelTimeSeries

from metobs_toolkit.backend_collection.loggingmodule import log_entry
from metobs_toolkit.backend_collection.dev_collection import copy_doc
from metobs_toolkit.backend_collection.dataframe_constructors import station_df

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
        self.obsdata = {
            sensor_data.obstype.name: sensor_data for sensor_data in all_sensor_data
        }

        # Extra extracted data
        self._modeldata = []  # list of ModelTimeSeries

    def _id(self) -> str:
        """A physical unique id.

        In the __add__ methods, if the id of two instances differs, adding is
        a regular concatenation.
        """
        return f"{self.name}"

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

    def __add__(self, other: "Station") -> "Station":
        """
        Combine two Station instances with the same _id.

        Joining takes the _id() of underlying metobs objects into account. When
        a combination of two objects with the same _id is encountered, the
        addition is handled by the __add__ method of that class.

        Parameters
        ----------
        other : Station
            The Station instance to add to the current Station.

        Returns
        -------
        Station
            A new Station instance containing merged Sensordata, Site and ModelDataTimeseries.

        Warning
        -------
        All progress on outliers and gaps will be lost! Outliers and gaps are reset.
        This is necessary to be able to join SensorData with other time resolutions.

        Warning
        -------
        When two Stations are joined with an overlap in sensortype, and
        timestamps, the values (if not-NaN in other) are taken from *other*.

        Examples
        --------
        >>> # Assume sta1_A and sta1_B Stations, representing the same station.
        >>> sta1 = sta1_A + sta1_B
        """
        if not isinstance(other, Station):
            raise MetObsAdditionError("Can only add Station to Station.")
        if self.name != other.name:
            raise MetObsAdditionError("Cannot add Station with different names.")

        #  ----  Merge site ----
        merged_site = self.site + other.site

        # --- Merge sensordata ----

        # use collection merge
        merged_sensorlist = join_collections(
            col_A=self.sensordata.values(), col_B=other.sensordata.values()
        )

        # --- Merge Modeldata ----
        merged_modeldatalist = join_collections(
            col_A=self.modeldata,
            col_B=other.modeldata,
        )

        # Construct a new station
        new_sta = Station(
            stationname=self.name, site=merged_site, all_sensor_data=merged_sensorlist
        )

        for moddata in merged_modeldatalist:
            new_sta.add_to_modeldata(new_modeltimeseries=moddata, force_update=True)

        return new_sta

    def __repr__(self):
        """Return a string representation for debugging."""
        return f"{type(self).__name__}(id={self._id()})"

    @log_entry
    def copy(self, deep: bool = True) -> "Station":
        """
        Return a copy of the Station.

        Parameters
        ----------
        deep : bool, optional
            If True, perform a deep copy. Default is True.

        Returns
        -------
        Station
            The copied station.
        """

        if deep:
            return copy.deepcopy(self)
        return copy.copy(self)

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

    @log_entry
    def get_sensor(self, obstype: str) -> "SensorData":  # type: ignore #noqa: F821
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

    @copy_doc(station_to_xr)
    @log_entry
    def to_xr(self) -> "xarray.Dataset":
        return station_to_xr(self, fmt_datetime_coordinate=True)

    @log_entry
    def to_netcdf(self, filepath: str, **kwargs) -> None:
        """
        Save the Station as a netCDF file.

        This method converts the Station to an xarray Dataset and saves it as a
        netCDF file.

        Parameters
        ----------
        filepath : str
            Path where the netCDF file will be saved.
        **kwargs
            Additional keyword arguments passed to xarray.Dataset.to_netcdf().
            Common options include:
            - format : str, netCDF format ('NETCDF4', 'NETCDF4_CLASSIC', 'NETCDF3_64BIT', 'NETCDF3_CLASSIC')
            - engine : str, netCDF engine to use ('netcdf4', 'scipy', 'h5netcdf')
            - encoding : dict, variable-specific encoding parameters

        Examples
        --------
        >>> station.to_netcdf('station_data.nc')
        >>> station.to_netcdf('data.nc', format='NETCDF4_CLASSIC')

        Notes
        -----
        This method is an export method. It is not possible to convert a netCDF
        to a metobs_toolkit.Station object.
        """

        # Convert to xarray Dataset
        ds = self.to_xr()
        # Save to netCDF
        ds.to_netcdf(filepath, **kwargs)

    @log_entry
    def to_parquet(self, target_file: Union[str, Path], **kwargs) -> None:
        """
        Save the station observations to a parquet file.

        The DataFrame returned by the `.df` property is written to a parquet file.
        This includes all observations with their QC labels (or gapfill labels) for this station.

        Parameters
        ----------
        target_file : str or Path
            The file path where the parquet file will be saved.
        **kwargs
            Additional keyword arguments to pass to pandas.DataFrame.to_parquet().

        Returns
        -------
        None

        See Also
        --------
        Station.df : The DataFrame property that is written to file.
        Station.to_csv : Save station data to CSV format.
        """
        df = self.df
        df.to_parquet(target_file, **kwargs)

    @log_entry
    def to_csv(self, target_file: Union[str, Path], **kwargs) -> None:
        """
        Save the station observations to a CSV file.

        The DataFrame returned by the `.df` property is written to a CSV file.
        This includes all observations with their QC labels (or gapfill labels) for this station.

        Parameters
        ----------
        target_file : str or Path
            The file path where the CSV file will be saved.
        **kwargs
            Additional keyword arguments to pass to pandas.DataFrame.to_csv().

        Returns
        -------
        None

        See Also
        --------
        Station.df : The DataFrame property that is written to file.
        Station.to_parquet : Save station data to parquet format.
        """
        df = self.df
        df.to_csv(target_file, **kwargs)

    @copy_doc(station_df)
    @property
    def df(self) -> pd.DataFrame:
        return station_df(self)

    @property
    def outliersdf(self) -> pd.DataFrame:
        """
        Construct a DataFrame representation of all the outliers.

        Outliers are the observations that are flagged by the performed quality control.

        Returns
        -------
        pandas.DataFrame
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
        """
        Construct a DataFrame representation of all the gaps.

        Returns
        -------
        pandas.DataFrame
            A DataFrame with columns ['value', 'label', 'details'], representing
            the value, the gap label, and details of the gap record.
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
        """
        Construct a DataFrame representation of metadata.

        Metadata is the information related to the sensors, that does not change over time.
        The metadata is extracted from the site instance.

        Returns
        -------
        pandas.DataFrame
            A DataFrame with the station names as index, and the metadata as columns.
        """

        return self.site.metadf

    @property
    def modeldatadf(self) -> pd.DataFrame:
        """
        Construct a DataFrame representation of all the present model data.

        Model data is stored as `ModelTimeSeries` instances, and is set as an attribute of
        a `Station`.

        Returns
        -------
        pandas.DataFrame
            A DataFrame with columns ['value', 'details'], representing
            the value, and details of the corresponding model data.
        """
        concatlist = []
        for modeldata in self.modeldata:
            df = (
                modeldata.df.assign(obstype=modeldata.modelobstype.name)
                .assign(
                    details=f"{modeldata.modelname}:{modeldata.modelvariable} converted from {modeldata.modelobstype.model_unit} -> {modeldata.modelobstype.std_unit}"
                )
                .assign(modelname=modeldata.modelname)
                .assign(modelvariable=modeldata.modelvariable)
                .reset_index()
                .set_index(["datetime", "obstype"])
            )
            concatlist.append(df)
        combdf = save_concat(concatlist)
        if combdf.empty:
            combdf = pd.DataFrame(
                columns=["value", "details", "modelname", "modelvariable"],
                index=pd.MultiIndex(
                    levels=[[], []], codes=[[], []], names=["datetime", "obstype"]
                ),
            )
        # formatting
        combdf = combdf[["value", "details", "modelname", "modelvariable"]]
        combdf.sort_index(inplace=True)
        return combdf

    @log_entry
    def get_modeltimeseries(
        self,
        obstype: str,
        modelname: str | None = None,
        modelvariable: str | None = None,
    ) -> "ModelTimeSeries":  # type: ignore #noqa: F821
        """Get the ModelTimeSeries instance for a specific observation type.

        Parameters
        ----------
        obstype : str
            The observation type to retrieve.
        modelname : str, optional
            The model name to filter by. Use this parameter when multiple model
            data sources exist for the same observation type. If None, no
            filtering by model name is applied. The default is None.
        modelvariable : str, optional
            The model variable to filter by. Use this parameter when multiple
            model variables exist for the same observation type and model. If None,
            no filtering by model variable is applied. The default is None.

        Returns
        -------
        ModelTimeSeries
            The ModelTimeSeries instance for the specified observation type.

        Raises
        ------
        MetObsModelDataError
            If no model data is found for the specified parameters, or if
            multiple model data instances are found and additional filtering
            parameters (modelname, modelvariable) are needed to uniquely
            identify the target model data.
        """

        if not bool(self.modeldata):
            raise MetObsModelDataError(f"No model data found in {self}. ")

        # Get all candidates of modeldata for target obstype
        all_candidates = self.modeldata

        # 1. obstype filter
        candidates = [
            modeldata
            for modeldata in all_candidates
            if modeldata.modelobstype.name == str(obstype)
        ]

        if len(candidates) == 0:
            raise MetObsModelDataError(
                f"No model data found for {obstype} in {self}. The following modeldata is available: \n{all_candidates}"
            )

        # 2. modelname filter
        if modelname is not None:
            candidates = [
                modeldata
                for modeldata in candidates
                if modeldata.modelname == str(modelname)
            ]
            if len(candidates) == 0:
                raise MetObsModelDataError(
                    f"No model data found for {obstype} and {modelname} in {self}. The following modeldata is available: \n{all_candidates}"
                )

        # 3. modelvariable filter
        if modelvariable is not None:
            candidates = [
                modeldata
                for modeldata in candidates
                if str(modeldata.modelvariable) == str(modelvariable)
            ]
            if len(candidates) == 0:
                if modelname is None:
                    raise MetObsModelDataError(
                        f"No model data found for {obstype} with model variable: {modelvariable} in {self}. The following modeldata is available: \n{all_candidates}"
                    )
                else:
                    raise MetObsModelDataError(
                        f"No model data found for {obstype}, {modelname} and {modelvariable} in {self}. The following modeldata is available: \n{all_candidates}"
                    )

        if len(candidates) > 1:
            raise MetObsModelDataError(
                f"Multiple model data found for {obstype} in {self}. Please specify modelname and/or modelvariable. The following modeldata is available for {obstype}: \n{candidates}"
            )
        target = candidates[0]
        return target

    @property
    def start_datetime(self) -> pd.Timestamp:
        """
        Get the earliest start datetime from the observation data.

        Returns
        -------
        pd.Timestamp
            The earliest start datetime among all sensor data.
        """
        if bool(self.sensordata):
            mindt = min(
                [sensdata.start_datetime for sensdata in self.sensordata.values()]
            )
        else:
            # no sensordata, metadata only station
            mindt = pd.NaT
        return mindt

    @property
    def end_datetime(self) -> pd.Timestamp:
        """
        Get the latest end datetime from the observation data.

        Returns
        -------
        pd.Timestamp
            The maximum end datetime across all sensor data in the observation data.
        """

        if bool(self.sensordata):
            mindt = max(
                [sensdata.end_datetime for sensdata in self.sensordata.values()]
            )
        else:
            # no sensordata, metadata only station
            mindt = pd.NaT
        return mindt

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
        """
        Get a list of all the present observation types.

        Returns
        -------
        list
            A list of all the present observations in the station.
        """
        return sorted(list(self.sensordata.keys()))

    @log_entry
    def add_to_sensordata(
        self, new_sensordata: SensorData, force_update: bool = False
    ) -> None:
        """
        Add a new SensorData to the Station. This can only be done when
        the new_sensordata is not already present in self based on its id,
        or if force_update is true.

        Parameters
        ----------
        new_sensordata : SensorData
            The new sensor data to be added. Must be an instance of `SensorData`.
        force_update : bool, optional
            If True, overwrite existing sensor data for the same id. Default is False.

        Returns
        -------
        None
        """
        # Validate argument types
        if not isinstance(new_sensordata, SensorData):
            raise TypeError("new_sensordata must be an instance of SensorData.")
        present_sensordata_ids = [
            sensordat._id() for sensordat in self.sensordata.values()
        ]
        # Test if there is already sensordata for the same id available
        if (new_sensordata._id() in present_sensordata_ids) & (not force_update):
            raise MetObsDataAlreadyPresent(
                f"There is already a SensorData instance with id {new_sensordata._id()}, and force_update is False."
            )

        self.obsdata.update({new_sensordata.obstype.name: new_sensordata})

    @log_entry
    def add_to_modeldata(
        self, new_modeltimeseries: ModelTimeSeries, force_update: bool = False
    ) -> None:
        """
        Add a new ModelTimeSeries to the Station. This can only be done when
        the new_modeltieseries is not already present in self based on its id,
        or if force_update is true.

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
        # Validate argument types
        if not isinstance(new_modeltimeseries, ModelTimeSeries):
            raise TypeError(
                "new_modeltimeseries must be an instance of ModelTimeSeries."
            )
        present_modeldata_ids = [modeldat._id() for modeldat in self.modeldata]

        # Test if there is already model data for the same ID available
        if new_modeltimeseries._id() in present_modeldata_ids:
            if force_update:
                # Subset to all other id's
                logger.warning(
                    f"ModelTimeSeries with id {new_modeltimeseries._id()} already exists and will be overwritten due to force_update=True."
                )
                self._modeldata = [
                    md
                    for md in self._modeldata
                    if md._id() != new_modeltimeseries._id()
                ]
                self._modeldata.append(new_modeltimeseries)
            else:
                raise MetObsDataAlreadyPresent(
                    f"There is already a modeltimeseries instance with id {new_modeltimeseries._id()}, and force_update is False."
                )
        else:
            self._modeldata.append(new_modeltimeseries)

    @log_entry
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

        # Model data info
        infostr += printing.print_fmt_section("Modeldata info")
        if not bool(self.modeldata):
            infostr += printing.print_fmt_line("Station instance without model data.")
        else:
            for modeldata in self.modeldata:
                infostr += printing.print_fmt_line(f"{modeldata._id()}:", 1)
                infostr += modeldata._get_info_core(nident_root=2)

        if printout:
            print(infostr)
        else:
            return infostr

    @log_entry
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
        and transferring the nearest (in time) value to the target timestamp. There is no interpolation in time.

        Alignment restrictions can be specified by setting the shift_tolerance, which indicates the maximum shift
        (in time) that aligned timestamps can have.

        Since the origin (the first timestamp in the current observations) is in general not a suitable candidate
        to serve as origin for the target frequency, a suitable origin is deduced. This can be done by specifying the
        origin argument, or if None it is done automatically. This is done by simplifying the
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
        origin : pandas.Timestamp, or None, optional
            The origin timestamp (=first timestamp) of the target records. Can be a pandas.Timestamp or None.
            If a timezone-naive timestamp is given, it is assumed to be in UTC.
            If None, then a suitable origin is used by simplifying the current origin and making
            sure that origin_simplify_tolerance is met. Default is None.
        origin_simplify_tolerance : pandas.Timedelta, optional
            This is only used when origin is None. The tolerance for simplifying the origin
            alignment during resampling. Default is 4 minutes.

        Returns
        -------
        None

        Notes
        -----

        * If `target_obstype` is None, all sensors in `self.sensordata` will be
          resampled to the same frequency.
        * If `target_obstype` is specified, it must be a known observation type.
        Warning
        -------
        Since the gaps depend on the record’s frequency and origin, all gaps
        are removed and re-located. All progress in gap filling will be lost.

        Warning
        -------
        Cumulative tolerance errors can be introduced when this method is called multiple times.
        """
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
            # check if target obstype is known
            self._obstype_is_known_check(target_obstype)

            # resample
            self.sensordata[target_obstype].resample(
                target_freq=target_freq,
                shift_tolerance=shift_tolerance,
                origin=origin,
                origin_simplify_tolerance=origin_simplify_tolerance,
            )

    @log_entry
    def convert_outliers_to_gaps(
        self, all_observations: bool = False, obstype: str = "temp"
    ) -> None:
        """
        Convert outlier values in the observation data to gaps.

        This method replaces outlier values with gaps, effectively marking them as missing.
        In practice, this method is used so that generated gaps can be filled, to obtain a continuous time series.

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
        if all_observations:
            for sensor in self.sensordata.values():
                sensor.convert_outliers_to_gaps()

        else:
            self._obstype_is_known_check(obstype)
            self.get_sensor(obstype).convert_outliers_to_gaps()

    def _rename(self, targetname):
        # Note: Not for users, one could accidentally rename to another station in the dataset.
        # So --> only accessible as method in the dataset, that checks this possible error.

        # rename all
        self._name = str(targetname)
        self._site._stationname = str(targetname)
        for sensordat in self.sensordata.values():
            sensordat._rename(targetname)

    @log_entry
    def get_static_gee_point_data(
        self,
        geestaticdatasetmanager: GEEStaticDatasetManager,
        overwrite: bool = True,
        initialize_gee: bool = True,
    ):
        """
        Extract static data from GEE dataset at Station locations.

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
        if not isinstance(geestaticdatasetmanager, GEEStaticDatasetManager):
            raise ValueError(
                f"geestaticdataset should be an instance of GeeStaticDataset, not {type(geestaticdatasetmanager)}"
            )

        value = self.site.get_gee_point_metadata(
            geestaticdataset=geestaticdatasetmanager, initialize_gee=initialize_gee
        )
        if overwrite:
            self.site.set_geedata(geestaticdatasetmanager.name, value)

        return value

    @log_entry
    def get_LCZ(
        self,
        overwrite: bool = True,
        initialize_gee: bool = True,
        apply_seamask_fix: bool = True,
    ) -> str:
        """
        Retrieve Local Climate Zone (LCZ) for the station using Google Earth Engine (GEE).

        Parameters
        ----------
        overwrite : bool, optional
            If True, overwrite existing LCZ data if stored in the Site attribute. Default is True.
        initialize_gee : bool, optional
            If True, initialize the Google Earth Engine API before fetching data. Default is True.
        apply_seamask_fix: bool, optional
            The LCZ map is only defined over land, and thus locations in sea
            will have a LCZ of Nan. If this argument is set to True, Nan values
            return by the GEE call are converted to the LCZ-G (water) category.

        Returns
        -------
        str
            The LCZ of the station.

        Notes
        -----
        This method relies on the `get_static_gee_point_data` function and the
        `default_gee_datasets` dictionary to fetch the LCZ data.
        """
        lcz = self.get_static_gee_point_data(
            geestaticdatasetmanager=default_gee_datasets["LCZ"],
            overwrite=False,  # overwrite is done in this method
            initialize_gee=initialize_gee,
        )

        if apply_seamask_fix:
            if isinstance(lcz, str):
                pass  # already a valid LCZ class
            elif np.isnan(lcz):
                logger.warning(f"Seamask fix for LCZ applied for {self}")
                lcz = default_gee_datasets["LCZ"].class_map[17]  # LCZ-G water
            else:
                raise ValueError("Unexpected LCZ value")

        # update the Site instance
        if overwrite:
            self.site.set_LCZ(lcz)
        else:
            if self.site.flag_has_LCZ():
                # LCZ is present and overwrite == False (-> do not update)
                pass
            else:
                self.site.set_LCZ(lcz)
        return lcz

    @log_entry
    def get_altitude(
        self, overwrite: bool = True, initialize_gee: bool = True
    ) -> float:
        """
        Retrieve altitude for the station using Google Earth Engine (GEE).

        Parameters
        ----------
        overwrite : bool, optional
            If True, overwrite existing altitude data if stored in the Site attribute. Default is True.
        initialize_gee : bool, optional
            If True, initialize the Google Earth Engine API before fetching data. Default is True.

        apply_seamask_fix: bool, optional
            The LCZ map is only defined over land, and thus locations in sea
            will have a LCZ of Nan. If this argument is set to True, Nan values
            returned by the GEE call are converted to the LCZ-G (water) category.

        Returns
        -------
        float
            The altitude of the station.

        Notes
        -----
        This method relies on the `get_static_gee_point_data` function and the
        `default_gee_datasets` dictionary to fetch the altitude data.
        """
        altitude = self.get_static_gee_point_data(
            geestaticdatasetmanager=default_gee_datasets["altitude"],
            overwrite=False,  # overwrite is done in this method
            initialize_gee=initialize_gee,
        )

        # update the Site instance
        if overwrite:
            self.site.set_altitude(altitude)
        else:
            if self.site.flag_has_altitude():
                # altitude is present and overwrite == False (-> do not update)
                pass
            else:
                self.site.set_altitude(altitude)
        return altitude

    @log_entry
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
        This method makes use of the GEE API. Make sure that you have access and user rights to use the GEE API.

        Warning
        -------
        It can happen that for stations located on small islands, or close to the coast, the sea-mask is not used as a landcover fraction.
        """
        if not isinstance(geestaticdataset, GEEStaticDatasetManager):
            raise ValueError(
                f"geestaticdataset should be an instance of GeeStaticDataset, not {type(geestaticdataset)}"
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

    @log_entry
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
        This method makes use of the GEE API. Make sure that you have access and user rights to use the GEE API.
        """
        return self.get_static_gee_buffer_fraction_data(
            geestaticdataset=default_gee_datasets["worldcover"],
            buffers=buffers,
            aggregate=aggregate,
            overwrite=overwrite,
        )

    @log_entry
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
        """
        Extract time series data from GEE.

        Extracts time series (extraction at station locations) data from a Google Earth Engine (GEE) dynamic dataset
        for a specified time range and observation types.

        If the data request is small, GEE sends the data directly. If not, the data will be
        written to a CSV file and saved on your Google Drive. In this case, you can import the model data
        by using the `Dataset.import_gee_data_from_file()` method.

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
        If a timezone-unaware datetime is given as an argument, it is interpreted as if it has the same timezone as the observations.

        Warning
        -------
        This method makes use of the GEE API. Make sure that you have access and user rights to use the GEE API.

        Warning
        -------
        When extracting large amounts of data,
        the time series data will be written to a file and saved
        on your Google Drive. In this case, you can import the model data
        by using the Dataset.import_gee_data_from_file() method.

        Notes
        -----

        * The method creates `ModelTimeSeries` instances for each valid observation
          type in the extracted data and appends them to the station's model data.
        * If no data is returned by the GEE API request, a warning is logged and
          the method returns None.

        """
        # Check geedynamic dataset
        if not isinstance(geedynamicdatasetmanager, GEEDynamicDatasetManager):
            raise ValueError(
                f"geedynamicdataset should be an instance of GeeDynamicDataset, not {type(geedynamicdatasetmanager)}"
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
                    modelobstype=geedynamicdatasetmanager.modelobstypes[modelobscol],
                    datadtype=np.float32,
                    timezone="UTC",
                    modelname=geedynamicdatasetmanager.name,
                    modelvariable=geedynamicdatasetmanager.modelobstypes[
                        modelobscol
                    ].model_band,
                )
                # todo: duplicacy check
                self.add_to_modeldata(modeltimeseries, force_update=True)
            else:
                logger.info(
                    f"Skip {modelobscol} for creating a ModelTimeeries because of unknown obstype."
                )

        return df

    @log_entry
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
            The target observation to check. By default "temp"
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
        # argument validity checks
        self._obstype_is_known_check(target_obstype)

        self.get_sensor(target_obstype).gross_value_check(
            lower_threshold=lower_threshold, upper_threshold=upper_threshold
        )

    @log_entry
    def persistence_check(
        self,
        target_obstype: str = "temp",
        timewindow: Union[str, pd.Timedelta] = pd.Timedelta("60min"),
        min_records_per_window: int = 5,
    ) -> None:
        """
        Check if values are not constant in a moving time window.

        Perform a persistence check on a time series to identify periods where observations remain constant
        within a specified time window. If the values are constant, all records in the moving window are
        flagged as outliers.

        Parameters
        ----------
        target_obstype : str, optional
            The target observation to check. By default "temp"
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

        * This method modifies the outliers in place and does not return anything.
          You can use the `outliersdf` property to view all flagged outliers.
        * If the minimum number of records per window is locally not met, the function logs a warning and skips
          the persistence check.
        * This function can be computationally expensive for large datasets or small time windows.
        * The repetitions check is similar to the persistence check, but not identical.
          The persistence check uses thresholds that are meteorologically based (i.e. the moving window is defined by a duration),
          in contrast to the repetitions check whose thresholds are instrumentally based (i.e. the "window" is defined by a number of records.)

        Warnings
        -------
        If the minimum number of records per window is not met over the full time series, a warning is logged, and the function
        returns an empty DatetimeIndex.
        """
        # argument checks
        self._obstype_is_known_check(target_obstype)
        timewindow = fmt_timedelta_arg(timewindow)

        # apply check on the sensordata
        self.get_sensor(target_obstype).persistence_check(
            timewindow=timewindow, min_records_per_window=min_records_per_window
        )

    @log_entry
    def repetitions_check(
        self, target_obstype: str = "temp", max_N_repetitions: int = 5
    ) -> None:
        """
        Test if an observation changes after a number of repetitions.

        Perform a check that tests if the observation changes after a number of repetitions.
        If a value is repeated more than the specified number of times, all the repeated
        records are flagged as outliers.

        Be aware that the performance of this check depends on the `max_N_repetitions`
        and the time resolution of the observations.

        Parameters
        ----------
        target_obstype : str, optional
            The target observation to check. By default "temp"
        max_N_repetitions : int
            The maximum number of repetitions allowed before the records are flagged as outliers.
            If the number of repetitions exceeds this value, all repeated records are flagged as outliers. The default is 5.

        Returns
        -------
        None

        Notes
        -----

        * This method modifies the outliers in place and does not return anything.
          You can use the `outliersdf` property to view all flagged outliers.
        * The repetitions check is similar to the persistence check, but not identical.
          The persistence check uses thresholds that are meteorologically based (i.e. the moving window is defined by a duration),
          in contrast to the repetitions check whose thresholds are instrumentally based (i.e. the "window" is defined by a number of records.)

        Warnings
        -------
        If the minimum number of records per window is not met over the full time series, a warning is logged, and the function
        returns an empty DatetimeIndex.
        """
        # argument checks
        self._obstype_is_known_check(target_obstype)

        # apply check on the sensordata
        self.get_sensor(target_obstype).repetitions_check(
            max_N_repetitions=max_N_repetitions
        )

    @log_entry
    def step_check(
        self,
        target_obstype: str = "temp",
        max_increase_per_second: Union[int, float] = 8.0 / 3600.0,
        max_decrease_per_second: Union[int, float] = -10.0 / 3600.0,
    ) -> None:
        """
        Check for 'spikes' and 'dips' in a time series.

        Test if observations do not produce spikes in the time series. The maximum
        allowed increase and decrease per second is set in the argument,
        and is tested for each record (with respect to the previous record).

        If the difference between two consecutive records (i.e., the spike/dip) is larger than the
        threshold, the record is flagged as an outlier.

        Parameters
        ----------
        target_obstype : str, optional
            The target observation to check. By default "temp"
        max_increase_per_second : int or float, >0, optional
            The maximum allowed increase (per second). This value is extrapolated to the time resolution of records.
            This value must be positive!  The default is 8.0/3600.0
        max_decrease_per_second : int or float, <0, optional
            The maximum allowed decrease (per second). This value is extrapolated to the time resolution of records.
            This value must be negative! The default is -10.0/3600.0

        Returns
        -------
        None

        Notes
        -----

        * This method modifies the outliers in place and does not return anything.
          You can use the `outliersdf` property to view all flagged outliers.
        * In general, for temperatures, the decrease threshold is set less stringent than the increase
          threshold. This is because a temperature drop is meteorologically more
          common than a sudden increase which is often the result of a radiation error.

        """
        # argument checks
        self._obstype_is_known_check(target_obstype)

        # apply check on the sensordata
        self.get_sensor(target_obstype).step_check(
            max_increase_per_second=max_increase_per_second,
            max_decrease_per_second=max_decrease_per_second,
        )

    @log_entry
    def window_variation_check(
        self,
        target_obstype: str = "temp",
        timewindow: pd.Timedelta = pd.Timedelta("1h"),
        min_records_per_window: int = 3,
        max_increase_per_second: Union[int, float] = 8.0 / 3600,
        max_decrease_per_second: Union[int, float] = -10.0 / 3600,
    ) -> None:
        """
        Test if the increase/decrease in a time window exceeds a threshold.

        Checks if the variation of observations in time
        does not exceed a threshold. This is done by applying a moving window
        over the time series. The moving window is defined by a duration (timewindow),
        and tested if the window contains at least a minimum number of records.

        If the observations in the window increase/decrease more than a threshold, all
        observations in the window are flagged as outliers. The threshold is defined by the
        maximum increase/decrease per second multiplied by the window size in seconds.

        Parameters
        ----------
        target_obstype : str, optional
            The target observation to check. By default "temp"
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

        * This method modifies the outliers in place and does not return anything.
          You can use the `outliersdf` property to view all flagged outliers.
        * In general, for temperatures, the decrease threshold is set less stringent than the increase
        * In general, for temperatures, the decrease threshold is set less stringent than the increase
          threshold. This is because a temperature drop is meteorologically more
          common than a sudden increase which is often the result of a radiation error.
        * A suitable value for the min_records_per_window depends on the time resolution of the records and the window size.
        * This check is similar to the step check, but not identical. The step check a maximum allowed increase/decrease
          with respect to the previous value. The window variation check uses a moving window to test the maximum allowed variation.

        """
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

    @log_entry
    def get_qc_stats(
        self, target_obstype: str = "temp", make_plot: bool = True
    ) -> pd.DataFrame:
        """
        Generate quality control (QC) frequency statistics.

        This method calculates the frequency statistics for various QC checks
        applied, and other labels. The order of checks is taken into
        account.

        Frequency of labels is computed based on the set of all labels (for all
        records including gaps). The effectiveness of a check is shown by
        the frequency of outliers with respect to the number of records that were given
        to the check (thus taking into account the order of checks).

        The frequencies are returned in a dataframe, and can be plotted
        as pie charts.

        Parameters
        ----------
        target_obstype : str, optional
            The target observation type for which to compute frequency statistics, by default "temp".
        make_plot : bool, optional
            If True, a figure with pie charts representing the frequencies is generated. The default is True.

        Returns
        -------
        pandas.DataFrame
            A DataFrame containing the QC frequency statistics. The DataFrame
            has a multi-index with the station name and QC check label, and
            includes the following columns:

            * `N_all`: Total number of records in the dataset (including gaps).
            * `N_labeled`: Number of records with the specific label.
            * `N_checked`: Number of records checked for the specific QC check.

        """
        # argument checks
        self._obstype_is_known_check(target_obstype)
        # get freq statistics
        qc_df = self.get_sensor(target_obstype).get_qc_freq_statistics()

        if make_plot:
            plotdf = qc_df.reset_index().drop(columns=["name"]).set_index("qc_check")

            fig = plotting.qc_overview_pies(df=plotdf)
            fig.suptitle(
                f"QC frequency statistics of {target_obstype} on Station level: {self.name}."
            )
            return fig
        else:
            return qc_df

    @log_entry
    def make_plot_of_modeldata(
        self,
        obstype: str = "temp",
        modelname: str | None = None,
        modelvariable: str | None = None,
        linecolor: Union[str, None] = None,
        title: Union[str, None] = None,
        linestyle: str = "--",
        ax: Union[None, Axes] = None,
        figkwargs: dict = {},
    ) -> Axes:
        """
        Generate a time series plot of model data for a specific observation type.

        Parameters
        ----------
        obstype : str, optional
            The type of observation to plot model data for, by default "temp".
        modelname : str, optional
            The model name to filter by when multiple model data sources exist
            for the same observation type. If None, no filtering by model name
            is applied. The default is None.
        modelvariable : str, optional
            The model variable to filter by when multiple model variables exist
            for the same observation type and model. If None, no filtering by
            model variable is applied. The default is None.
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
        # test if the obstype has model data
        trg_modeltimeseries = self.get_modeltimeseries(
            obstype=obstype, modelname=modelname, modelvariable=modelvariable
        )

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
            legend_prefix=f"{trg_modeltimeseries.modelname}:{trg_modeltimeseries.modelvariable}@",
        )
        # Styling
        obstypeinstance = trg_modeltimeseries.modelobstype

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

    @log_entry
    def make_plot(
        self,
        obstype: str = "temp",
        colorby: Literal["station", "label"] = "label",
        show_modeldata: bool = False,
        modelobstype: str = None,
        modeldata_kwargs: dict = {},
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
        modelobstype: str, optional
            The name of the ModelObstype to plot. It is only used if show_modeldata is True. If None, it is set equal to obstype.
            The default is None.
        colorby : {"station", "label"}, optional
            Determines how the data is colored in the plot.

            * "station": Colors by station.
            * "label": Colors by label (the labels refer to the status of a record).

            Default is "label".
        show_modeldata : bool, optional
            If True, includes model data (of the same obstype) if present, in the plot. Default is False.
        modeldata_kwargs: dict, optional
            Additional keyword arguments passed to make_plot_of_modeldata(), by default an empty dictionary. Use it for example to specify modelname if multiple model data is available.
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

        * The method checks if the specified `obstype` is known before proceeding.
        * The plot can include observational data, model data, or both.
        * The x-axis timestamps are formatted according to the timezone of the data.

        """
        # test if obstype have sensordata
        self._obstype_is_known_check(obstype)

        # Create new axes if needed
        if ax is None:
            ax = plotting.create_axes(**figkwargs)

        if show_modeldata:
            if modelobstype is None:
                modelobstype = obstype
            if linecolor is None:
                colormap = plotting.create_categorical_color_map([self.name])
            else:
                colormap = {self.name: linecolor}

            ax = self.make_plot_of_modeldata(
                obstype=modelobstype,
                linecolor=linecolor,
                ax=ax,
                figkwargs=figkwargs,
                title=title,
                **modeldata_kwargs,
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

    @log_entry
    def fill_gaps_with_raw_modeldata(
        self,
        target_obstype: str,
        overwrite_fill: bool = False,
        modelname: str | None = None,
        modelvariable: str | None = None,
    ) -> None:
        """
        Fill the gap(s) using model data without correction.

        This method fills all the gaps of a specific *target_obstype*, by directly interpolating
        the model data to the missing records.

        Parameters
        ----------
        target_obstype :  str
            The target obstype to fill the gaps for.
        overwrite_fill : bool, optional
            If True, the status of a `gap` and present gapfill info will be ignored and overwritten.
            If False, only gaps without gapfill data are filled. Defaults to False.
        modelname : str, optional
            The model name to filter by when multiple model data sources exist
            for the same observation type. If None, no filtering by model name
            is applied. The default is None.
        modelvariable : str, optional
            The model variable to filter by when multiple model variables exist
            for the same observation type and model. If None, no filtering by
            model variable is applied. The default is None.

        Returns
        -------
        None

        Notes
        -----
        A schematic description of the raw model data gap fill:

        #. Check if the target_obstype is knonw, and if the corresponding modeldata is present.
        #. Iterate over the gaps of the target_obstype.
        #. Check the compatibility of the `ModelTimeSeries` with the `gap`.
        #. Ensure both the `ModelTimeSeries` and `gap` have the same timezone.
        #. Interpolate the model data to match the missing records in the gap.
        #. Update the `gap` attributes with the interpolated values, labels, and details.

        """
        # obstype check
        self._obstype_is_known_check(obstype=target_obstype)

        # Get modeltimeseries
        modeltimeseries = self.get_modeltimeseries(
            target_obstype, modelname=modelname, modelvariable=modelvariable
        )

        # fill the gaps
        self.get_sensor(target_obstype).fill_gap_with_modeldata(
            modeltimeseries=modeltimeseries, method="raw", overwrite_fill=overwrite_fill
        )

    @log_entry
    def fill_gaps_with_debiased_modeldata(
        self,
        target_obstype: str,
        leading_period_duration: Union[pd.Timedelta, str] = pd.Timedelta("24h"),
        min_leading_records_total: int = 60,
        trailing_period_duration: Union[pd.Timedelta, str] = pd.Timedelta("24h"),
        min_trailing_records_total: int = 60,
        overwrite_fill: bool = False,
        modelname: str | None = None,
        modelvariable: str | None = None,
    ) -> None:
        """
        Fill the gaps using model data corrected for the bias.

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
        modelname : str, optional
            The model name to filter by when multiple model data sources exist
            for the same observation type. If None, no filtering by model name
            is applied. The default is None.
        modelvariable : str, optional
            The model variable to filter by when multiple model variables exist
            for the same observation type and model. If None, no filtering by
            model variable is applied. The default is None.

        Returns
        -------
        None

        Notes
        -----
        A schematic description of the debiased model data gap fill:

        #. Check if the target_obstype is knonw, and if the corresponding modeldata is present.
        #. Iterate over the gaps of the target_obstype.
        #. Check the compatibility of the `ModelTimeSeries` with the `gap`.
        #. Construct a leading and trailing sample, and test if they meet the required conditions.
        #. Compute the bias of the modeldata (combine leading and trailing samples).
        #. Fill the gap records by using raw (interpolated) modeldata that is corrected by subtracting the bias.
        #. Update the `gap` attributes with the interpolated values, labels, and details.

        """

        # special formatters
        leading_period_duration = fmt_timedelta_arg(leading_period_duration)
        trailing_period_duration = fmt_timedelta_arg(trailing_period_duration)

        # obstype check
        self._obstype_is_known_check(obstype=target_obstype)

        # Get modeltimeseries
        modeltimeseries = self.get_modeltimeseries(
            target_obstype, modelname=modelname, modelvariable=modelvariable
        )

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

    @log_entry
    def fill_gaps_with_diurnal_debiased_modeldata(
        self,
        target_obstype: str,
        leading_period_duration: Union[pd.Timedelta, str] = pd.Timedelta("24h"),
        trailing_period_duration: Union[pd.Timedelta, str] = pd.Timedelta("24h"),
        min_debias_sample_size: int = 6,
        overwrite_fill: bool = False,
        modelname: str | None = None,
        modelvariable: str | None = None,
    ) -> None:
        """
        Fill the gaps using model data corrected for the diurnal bias.

        This method fills the gap using model data corrected for its diurnal bias.
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
        modelname : str, optional
            The model name to filter by when multiple model data sources exist
            for the same observation type. If None, no filtering by model name
            is applied. The default is None.
        modelvariable : str, optional
            The model variable to filter by when multiple model variables exist
            for the same observation type and model. If None, no filtering by
            model variable is applied. The default is None.

        Returns
        -------
        None

        Notes
        -----
        A schematic description of the diurnal debiased model data gap fill:

        #. Check if the target_obstype is knonw, and if the corresponding modeldata is present.
        #. Iterate over the gaps of the target_obstype.
        #. Check the compatibility of the `ModelTimeSeries` with the `gap`.
        #. Construct a leading and trailing sample, and test if they meet the required conditions.
           The required conditions are tested by testing the samplesizes per hour, minute and second for the leading + trailing periods.
        #. A diurnal bias is computed by grouping to hour, minute and second, and averaging the biases.
        #. Fill the gap records by using raw (interpolated) modeldata that is corrected by subtracting the coresponding diurnal bias.
        #. Update the `gap` attributes with the interpolated values, labels, and details.

        Notes
        -----
        Note that a suitable `min_debias_sample_size` depends on the sizes of the
        leading- and trailing periods, and also on the time resolution gap (=time resolution of the corresponding SensorData).

        References
        ----------
        Jacobs A, et. al. (2024) `Filling gaps in urban temperature observations by debiasing ERA5 reanalysis data <https://doi.org/10.1016/j.uclim.2024.102226>`_
        """
        # special formatters
        leading_period_duration = fmt_timedelta_arg(leading_period_duration)
        trailing_period_duration = fmt_timedelta_arg(trailing_period_duration)

        # obstype check
        self._obstype_is_known_check(obstype=target_obstype)

        # Get modeltimeseries
        modeltimeseries = self.get_modeltimeseries(
            target_obstype, modelname=modelname, modelvariable=modelvariable
        )

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

    @log_entry
    def fill_gaps_with_weighted_diurnal_debiased_modeldata(
        self,
        target_obstype: str,
        leading_period_duration: Union[pd.Timedelta, str] = pd.Timedelta("24h"),
        trailing_period_duration: Union[pd.Timedelta, str] = pd.Timedelta("24h"),
        min_lead_debias_sample_size: int = 2,
        min_trail_debias_sample_size: int = 2,
        overwrite_fill=False,
        modelname: str | None = None,
        modelvariable: str | None = None,
    ):
        """
        Fill the gaps using a weighted sum of model data corrected for the diurnal bias and weights with respect to the start of the gap.

        This method fills the gaps using model data corrected for its diurnal bias.
        The diurnal bias is a bias that is estimated for each timestamp in the leading
        and trailing period (separately). For both periods separately, all biases are averaged over hour, minute and second, to
        obtain a diurnal bias (for each timestamp).

        In addition, a normalized weight is computed for each gap record indicating the distance (in time) to
        the start and end of the gap. The correction applied on the interpolated (in time) model data, is
        thus a weighted sum of corrections coming from both the leading and trailing period.

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
            The minimum number of leading samples required for bias estimation. If this condition is not met, the gap
            is not filled. The default is 2.
        min_trail_debias_sample_size : int, optional
            The minimum number of trailing samples required for bias estimation. If this condition is not met, the gap
            is not filled. The default is 2.
        overwrite_fill : bool, optional
            If True, the status of a `gap` and present gapfill info will be ignored and overwritten.
            If False, only gaps without gapfill data are filled. The default is False.
        modelname : str, optional
            The model name to filter by when multiple model data sources exist
            for the same observation type. If None, no filtering by model name
            is applied. The default is None.
        modelvariable : str, optional
            The model variable to filter by when multiple model variables exist
            for the same observation type and model. If None, no filtering by
            model variable is applied. The default is None.

        Returns
        -------
        None

        Notes
        -----
        A schematic description of the weighted diurnal debiased model data gap fill:

        #. Check if the target_obstype is knonw, and if the corresponding modeldata is present.
        #. Iterate over the gaps of the target_obstype.
        #. Check the compatibility of the `ModelTimeSeries` with the `gap`.
        #. Construct a leading and trailing sample, and test if they meet the required conditions. The required conditions are tested by testing the samplesizes per hour, minute and second for the leading and trailing periods (seperately).
        #. A leading and trailing set of diurnal biases are computed by grouping to hour, minute and second, and averaging the biases.
        #. A weight is computed for each gap record, that is the normalized distance to the start and end of the gap.
        #. Fill the gap records by using raw (interpolated) modeldata is corrected by a weighted sum the coresponding diurnal bias for the lead and trail periods.
        #. Update the `gap` attributes with the interpolated values, labels, and details.

        Notes
        -----
        Note that a suitable `min_debias_sample_size` depends on the sizes of the
        leading- and trailing periods, and also on the time resolution gap (=time resolution of the corresponding SensorData).

        References
        ----------
        Jacobs A, et. al. (2024) `Filling gaps in urban temperature observations by debiasing ERA5 reanalysis data <https://doi.org/10.1016/j.uclim.2024.102226>`_
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
        modeltimeseries = self.get_modeltimeseries(
            target_obstype, modelname=modelname, modelvariable=modelvariable
        )

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

    @log_entry
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
        """
        Fill the gap(s) using interpolation of SensorData.

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
            set accordingly to the method (higher order interpolation techniques require more leading and trailing anchors). The default is "time".
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
            Extra arguments that are passed to pandas.DataFrame.interpolate() structured in a dict. Defaults to {}.

        Returns
        -------
        None

        Notes
        -----
        A schematic description:

        #. Iterate over all gaps related to the target obstype.
        #. Get the leading and trailing periods of the gap.
        #. Check if the leading and trailing periods are valid.
        #. Create a combined DataFrame with the leading, trailing, and gap data.
        #. Interpolate the missing records using the specified method.
        #. Update the gap attributes with the interpolated values, labels, and details.

        Note
        -----
        The impact of `max_consec_fill` is highly dependent on the resolution
        of your records.

        Note
        -----
        If you want to use a higher-order method of interpolation, make sure to
        increase the `n_leading_anchors` and `n_trailing_anchors` accordingly.
        For example, for a cubic interpolation, you need at least 2 leading and 2 trailing anchors.
        """
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
