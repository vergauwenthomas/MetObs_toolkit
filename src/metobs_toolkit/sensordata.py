from __future__ import annotations
import logging
from typing import Literal, Union, TYPE_CHECKING


import copy
import numpy as np
import pandas as pd


from metobs_toolkit.backend_collection.dev_collection import copy_doc
from metobs_toolkit.backend_collection.datetime_collection import (
    timestamps_to_datetimeindex,
    convert_timezone,
    to_timedelta,
)

from metobs_toolkit.backend_collection.df_helpers import (
    save_concat,
    convert_to_numeric_series,
)
from metobs_toolkit.gf_collection.overview_df_constructors import (
    sensordata_gap_status_overview_df,
)
from metobs_toolkit.qc_collection.overview_df_constructor import (
    sensordata_qc_overview_df,
)
from metobs_toolkit.settings_collection import Settings
from metobs_toolkit.xrconversions import sensordata_to_xr
from metobs_toolkit.timestampmatcher import TimestampMatcher
from metobs_toolkit.obstypes import Obstype
from metobs_toolkit.gap import Gap
import metobs_toolkit.qc_collection as qc
from metobs_toolkit.backend_collection.errorclasses import (
    MetObsQualityControlError,
    MetObsAdditionError,
    MetObsInternalError,
)
from metobs_toolkit.qcresult import (
    QCresult,
    pass_cond,
    flagged_cond,
    unmet_cond,
    saved_cond,
    unchecked_cond
)
from metobs_toolkit.plot_collection import sensordata_simple_pd_plot
import metobs_toolkit.backend_collection.printing_collection as printing

from metobs_toolkit.backend_collection.decorators import log_entry
from metobs_toolkit.backend_collection.dataframe_constructors import sensordata_df

# add all imports only for type checking, but not for runtime
if TYPE_CHECKING:
    from dateutil.tz import tzfile
    from datetime import tzinfo
    from matplotlib.pyplot import Axes
    from xarray import Dataset as xrDataset
    from metobs_toolkit.modeltimeseries import ModelTimeSeries


logger = logging.getLogger("<metobs_toolkit>")


class SensorData:
    """
    Holds data for one station for one sensor.

    Parameters
    ----------
    stationname : str
        Name of the station.
    datarecords : np.ndarray
        Array of data records.
    timestamps : np.ndarray
        Array of timestamps.
    obstype : Obstype
        Observation type.
    datadtype : type, optional
        Data type of the records, by default np.float32.
    timezone : str or tzinfo, optional
        Timezone identifier for the timestamps (e.g., 'UTC', 'Europe/Amsterdam'), by default 'UTC'.
    freq_estimation_method : {'highest', 'median'}, optional
        Method to estimate frequency, by default 'median'.
    freq_estimation_simplify_tolerance : pandas.Timedelta or str, optional
        Tolerance for frequency estimation simplification, by default pandas.Timedelta('1min').
    origin_simplify_tolerance : pandas.Timedelta or str, optional
        Tolerance for origin simplification, by default pandas.Timedelta('1min').
    timestamp_tolerance : pandas.Timedelta or str, optional
        Tolerance for timestamp matching, by default pandas.Timedelta('4min').
    """

    # Class variable for internal timezone storage
    _target_tz: str = Settings.get("store_tz")

    def __init__(
        self,
        stationname: str,
        datarecords: np.ndarray,
        timestamps: np.ndarray,
        obstype: Obstype,
        datadtype: type = np.float32,
        timestamps_tz: Union[tzfile | tzinfo | str] = "UTC",
        **setupkwargs,
    ):
        """Initialize SensorData for a single station and observation type.

        Parameters
        ----------
        stationname : str
            Name of the station this sensor belongs to.
        datarecords : numpy.ndarray
            Raw observation values.  Will be stored after unit conversion
            during setup.
        timestamps : numpy.ndarray
            Timestamps corresponding to *datarecords*.
        obstype : Obstype
            Observation type describing the measured variable and its units.
        datadtype : type, optional
            Numeric dtype for the stored series.  Default is
            :data:`numpy.float32`.
        timestamps_tz : str or tzinfo, optional
            Timezone of the provided *timestamps*.  Default is ``'UTC'``.
        **setupkwargs
            Additional keyword arguments forwarded to :meth:`_setup` (e.g.
            ``freq_estimation_method``, ``timestamp_tolerance``).
        """
        self._stationname = stationname
        self.obstype = obstype
        data = pd.Series(
            data=convert_to_numeric_series(datarecords, datadtype=datadtype).to_numpy(),
            index=_format_timestamp_index(
                timestamps, timestamps_tz
            ),  # Transformed as UTC
            name=obstype.name,
        )

        # data.index.name = "datetime"
        self.series = data  # datetime as index

        # outliers
        self.outliers = []  # List of QCresult
        self.outliers_values_bin = pd.Series(dtype=datadtype)  # Series of outlier values

        # gaps
        self.gaps = []  # list of Gap's

        # Setup the SensorData --> apply qc control on import, find gaps, unit conversions etc
        self._setup(**setupkwargs)
        logger.info("SensorData initialized successfully.")

    def _id(self) -> str:
        """A physical unique id.

        In the __add__ methods, if the id of two instances differs, adding is
        a regular concatenation.
        """
        return f"{self.stationname}_{self.obstype._id()}"

    def __eq__(self, other):
        """Check equality with another SensorData object."""
        if not isinstance(other, SensorData):
            return False
        return (
            self.stationname == other.stationname
            and self.df.equals(other.df)  # the df contains outliers and gaps as well
            and self.obstype == other.obstype
        )

    def __repr__(self):
        """Return a string representation for debugging."""
        return f"{type(self).__name__}(id={self._id()})"

    def __str__(self) -> str:
        """Return a string representation of the SensorData object."""
        return f"{self.obstype.name} data of station {self.stationname}."

    #TODO: update this method to handle QCresult outliers + outliers_values_bin
    def __add__(self, other: "SensorData") -> "SensorData":
        """
        Combine two SensorData objects for the same station and obstype.


        !!! The result contains all unique records, with preference to non-NaN values from 'other'.
        This makes the addition not strict associative.
        Outliers and gaps are concatenated.
        """
        if not isinstance(other, SensorData):
            raise MetObsAdditionError("Can only add SensorData to SensorData.")
        if self._id() != other._id():
            raise MetObsAdditionError(
                f"Cannot add SensorData for different IDs ({self._id()} != {other._id()})."
            )

        # NOTE! combining outliers and gaps is NOT TRIVIAL !! Frequency is not guaranteed equal
        # and a additional gap (inbetween) can occur. Outliers and Gaps are recomputed !!!

        # Think of what will happen if you merge two sensordata series, where
        # The first is QC'ed and at a hourly resolution. The second overlaps the
        # first but at a resolution of 5 minutes. This messes up the outliers and
        # gaps.
        if not self.outliersdf.empty:
            # rolback the outliervalues
            selfrecords = self.series.combine_first(self.outliersdf["value"])
        else:
            selfrecords = self.series

        if not other.outliersdf.empty:
            # rolback the outliervalues
            otherrecords = other.series.combine_first(other.outliersdf["value"])
        else:
            otherrecords = other.series

        # Align timezones if different
        otherrecords = otherrecords.tz_convert(self.tz)

        # Combine the series, preferring non-NaN from 'other'
        combined_series = selfrecords.combine_first(otherrecords)

        # NOTE! combining outliers and gaps is NOT TRIVIAL !! Frequency is not guaranteed equal
        # and a additional gap (inbetween) can occur. Outliers and Gaps are recomputed !!!

        # Think of what will happen if you merge two sensordata series, where
        # The first is QC'ed and at a hourly resolution. The second overlaps the
        # first but at a resolution of 5 minutes. This messes up the outliers and
        # gaps.

        if (
            bool(self.gaps)
            | bool(self.outliers)
            | bool(other.gaps)
            | bool(other.outliers)
        ):
            logger.warning(
                f"All stored outliers and gap info of {self} will not be present in the combined."
            )

        # Create new SensorData instance
        combined = SensorData(
            stationname=self.stationname,
            datarecords=combined_series.values,
            timestamps=combined_series.index.values,
            obstype=self.obstype + other.obstype,
            datadtype=combined_series.dtype,
            timestamps_tz=self.tz,
            apply_unit_conv=False,
        )

        return combined

    @log_entry
    def copy(self, deep: bool = True) -> "SensorData":
        """
        Return a copy of the sensordata.

        Parameters
        ----------
        deep : bool, optional
            If True, perform a deep copy. Default is True.

        Returns
        -------
        SensorData
            The copied Sensordata.
        """

        if deep:
            return copy.deepcopy(self)
        return copy.copy(self)

    @copy_doc(sensordata_df)
    @property
    def df(self) -> pd.DataFrame:
        return sensordata_df(self)

    @copy_doc(sensordata_to_xr)
    @log_entry
    def to_xr(self) -> xrDataset:
        return sensordata_to_xr(self, fmt_datetime_coordinate=True)

    #TODO: update this method to handle QCresult outliers
    @property
    def outliersdf(self) -> pd.DataFrame:
        """Return a DataFrame of the outlier records."""
        logger.debug("Creating outliers DataFrame for %s", self.stationname)
        to_concat = []
        for qcresult in self.outliers:
            checkdf = qcresult.create_outliersdf(subset_to_outliers=True)
            to_concat.append(checkdf)

        totaldf = save_concat(to_concat)
        #add the values column (values not stored in qcresult, only labels and details)
        totaldf["value"] = self.outliers_values_bin.loc[totaldf.index]

        if totaldf.empty:
            # return empty dataframe
            totaldf = pd.DataFrame(
                columns=["value", "label", "details"],
                index=pd.DatetimeIndex([], name="datetime"),
            )
        
        totaldf.sort_index(inplace=True)

        return totaldf[['value', 'label', 'details']] #fixed column order       

        
    @property
    def gapsdf(self) -> pd.DataFrame:
        """Return a DataFrame of the gap records."""
        to_concat = []
        if bool(self.gaps):
            for gap in self.gaps:
                to_concat.append(gap.df)
            return save_concat((to_concat)).sort_index()
        else:
            return pd.DataFrame(
                columns=["value", "label", "details"],
                index=pd.DatetimeIndex([], name="datetime"),
            )

    @copy_doc(sensordata_gap_status_overview_df)
    def gap_overview_df(self) -> pd.DataFrame:
        return sensordata_gap_status_overview_df(self)

    @copy_doc(sensordata_qc_overview_df)
    def qc_overview_df(self) -> pd.DataFrame:
        return sensordata_qc_overview_df(self)
    
    @property
    def stationname(self) -> str:
        """Return the name of the station this SensorData belongs to."""
        return self._stationname

    @property
    def tz(self):
        """Return the timezone of the stored timestamps."""
        # Should always be SensorData._target_tz
        cur_tz = self.series.index.tz
        if str(cur_tz) != SensorData._target_tz:
            raise MetObsInternalError(
                f"Internal timezone mismatch: {cur_tz.zone} != {SensorData._target_tz}"
            )

        return cur_tz

    @property
    def start_datetime(self) -> pd.Timestamp:
        """Return the start datetime of the series."""
        return self.series.index.min()

    @property
    def end_datetime(self) -> pd.Timestamp:
        """Return the end datetime of the series."""
        return self.series.index.max()

    @property
    def freq(self) -> pd.Timedelta:
        """Return the frequency of the series."""
        freq = pd.infer_freq(self.series.index)
        if freq is None:
            raise ValueError("Frequency could not be computed.")
        return to_timedelta(freq)

    def _setup(
        self,
        freq_estimation_method: Literal["highest", "median"] = "median",
        freq_estimation_simplify_tolerance: Union[pd.Timedelta, str] = pd.Timedelta(
            "1min"
        ),
        origin_simplify_tolerance: Union[pd.Timedelta, str] = pd.Timedelta("1min"),
        timestamp_tolerance: Union[pd.Timedelta, str] = pd.Timedelta("4min"),
        apply_invalid_check: bool = True,
        apply_dupl_check: bool = True,
        apply_unit_conv: bool = True,
        force_origin=None,
        force_freq=None,
        force_closing=None,
    ) -> None:
        """
        Set up the SensorData object.

        This includes:

        #. Find the duplicates (remove them from observations and add them to outliers).
        #. Invalid check (records that could not be typecast to numeric) are interpreted as gaps.
        #. Convert the values to standard units and update the observation types.
        #. Find gaps in the records (duplicates are excluded from the gaps).
        #. Get a frequency estimate per station.
        #. Initiate the gaps (find missing records).
        #. Add the missing records to the dataframe.

        Parameters
        ----------
        freq_estimation_method : str
            Method to estimate frequency.
        freq_estimation_simplify_tolerance : pandas.Timedelta or str
            Tolerance for frequency estimation simplification.
        origin_simplify_tolerance : pandas.Timedelta or str
            Tolerance for origin simplification.
        timestamp_tolerance : pandas.Timedelta or str
            Tolerance for timestamp matching.
        apply_invalid_check : bool, optional
            Whether to apply invalid value check, by default True.
        apply_dupl_check : bool, optional
            Whether to apply duplicate timestamp check, by default True.
        apply_unit_conv : bool, optional
            Whether to apply unit conversion, by default True.
        force_origin : optional
            Force a specific origin.
        force_freq : optional
            Force a specific frequency.
        force_closing : optional
            Force closing parameter.
        """

        if apply_dupl_check:
            # remove duplicated timestamps
            self.duplicated_timestamp_check()
            

        if apply_invalid_check:
            # get the records that are flagged by the
            if self.outliers and self.outliers[0].checkname == 'duplicated_timestamp':
                dup_outl_ti = self.outliers[0].get_outlier_timestamps()
            else:
                dup_outl_ti = pd.DatetimeIndex([])
            # invalid check (no qcresult, these timestamps are removed, and catched by gapcheck)
            valid_records = qc.drop_invalid_values(records=self.series,
                                                   skip_records=dup_outl_ti)
            self.series = valid_records
        
        

        if apply_unit_conv:
            # convert units to standard units
            self.convert_to_standard_units()

        # format to perfect time records
        timestamp_matcher = TimestampMatcher(orig_records=self.series)
        timestamp_matcher.make_equispaced_timestamps_mapper(
            freq_estimation_method=freq_estimation_method,
            freq_estimation_simplify_tolerance=freq_estimation_simplify_tolerance,
            origin_simplify_tolerance=origin_simplify_tolerance,
            timestamp_tolerance=timestamp_tolerance,
            force_closing=force_closing,
            force_origin=force_origin,
            force_freq=force_freq,
        )

        # update all the attributes holding data
        self.series = timestamp_matcher.target_records

        # update the outliers (replace the raw timestamps with the new)
        raw_datetime_map = timestamp_matcher.get_raw_map()
        for qcresult in self.outliers:
            qcresult.remap_timestamps(mapping=raw_datetime_map)

        # remap the outliers_values_bin timestamps to match the new equispaced timestamps
        if not self.outliers_values_bin.empty:
            # Drop outlier values whose timestamps are not in the new time grid
            self.outliers_values_bin = self.outliers_values_bin[
                self.outliers_values_bin.index.isin(raw_datetime_map.keys())
            ]
            self.outliers_values_bin.index = self.outliers_values_bin.index.map(
                lambda ts: raw_datetime_map[ts]
            )

        # create gaps

        if bool(self.gaps):
            logger.warning(
                "The present gaps are removed, new gaps are constructed for %s.", self
            )
            self.gaps = []

        # Construct gaps
        self.gaps = self._find_gaps(
            missingrecords=timestamp_matcher.gap_records,
            target_freq=pd.to_timedelta(timestamp_matcher.target_freq),
        )
        
    def _update_outliers(self,
                        qcresult: QCresult,
                        ) -> None:
        """Record a QCresult and mask the detected outlier timestamps as NaN.

        If a QCresult with the same ``checkname`` already exists, the new
        results are merged into the existing one so that repeated application
        of the same check does not create duplicate entries.  Specifically:

        * Flags in the existing result are updated for timestamps where the
          new result detected an outlier (``flagged``) that was previously
          ``passed`` or ``unchecked``.
        * Details are updated for newly flagged timestamps.
        * Only genuinely new outlier timestamps (not already present in
          :attr:`outliers_values_bin`) are added.

        Parameters
        ----------
        qcresult : QCresult
            Result object from a QC check containing flags and metadata.
        """
        new_outlier_ts = qcresult.get_outlier_timestamps()

        # Check if a QCresult with the same checkname already exists
        existing = None
        for qc in self.outliers:
            if qc.checkname == qcresult.checkname:
                existing = qc
                break

        if existing is not None:
            # Merge: update the existing QCresult with newly flagged timestamps
            newly_flagged = new_outlier_ts.difference(existing.get_outlier_timestamps())
            if not newly_flagged.empty:
                existing.flags.loc[newly_flagged] = qcresult.flags.loc[newly_flagged]
                existing.details.loc[newly_flagged] = qcresult.details.loc[newly_flagged]
        else:
            self.outliers.append(qcresult)

        # Only add values for outlier timestamps not already in the bin
        truly_new_ts = new_outlier_ts.difference(self.outliers_values_bin.index)
        if not truly_new_ts.empty:
            self.outliers_values_bin = pd.concat([self.outliers_values_bin,
                                                    self.series.loc[truly_new_ts]])

        # Convert the outlier timestamps to NaN in the series
        self.series.loc[new_outlier_ts] = np.nan
        
        
    def _find_gaps(self, missingrecords: pd.Series, target_freq: pd.Timedelta) -> list:
        """
        Identify gaps in the missing records based on the target frequency.

        Parameters
        ----------
        missingrecords : pd.Series
            A pandas Series containing the missing records with datetime index.
        target_freq : pandas.Timedelta
            The target frequency to identify gaps.

        Returns
        -------
        list
            A list of Gap objects representing the identified gaps.

        Raises
        ------
        TypeError
            If input types are incorrect.
        """

        missing = missingrecords.sort_index().to_frame()
        missing["diff"] = missing.index.to_series().diff()
        missing["gap_group"] = (missing["diff"] != target_freq).cumsum()
        gaps = []
        for _idx, gapgroup in missing.groupby("gap_group"):
            gap = Gap(
                gaprecords=pd.date_range(
                    gapgroup.index.min(),
                    gapgroup.index.max(),
                    freq=pd.to_timedelta(target_freq),
                ),
                obstype=self.obstype,
                stationname=self.stationname,
            )
            gaps.append(gap)
        return gaps

    def _rename(self, trgname: str) -> None:
        """Rename the station and update gaps accordingly."""
        self._stationname = str(trgname)
        for gap in self.gaps:
            gap.name = str(trgname)

    @log_entry
    #TODO: update this method to handle QCresult outliers
    def convert_outliers_to_gaps(self) -> None:
        """
        Convert all outliers to gaps.

        This method will convert all outliers to gaps. Doing so new gaps are constructed.

        Returns
        -------
        None.

        Warning
        -------
        All progress on present gaps is erased, since new gaps are constructed.
        Information on the value and QC flag of the outliers will be lost.
        """

        cur_freq = self.freq
        # Create holes for all the outliers timestamps
        self.series.loc[self.outliersdf.index] = np.nan

        # Flush the outliers
        logger.warning(f"Outliers are flushed for {self}!")
        self.outliers = []

        # Flush the gaps
        if bool(self.gaps):
            logger.warning(f"Flushing current gaps for {self}")
        self.gaps = []

        # Finding new gaps
        self.gaps = self._find_gaps(
            missingrecords=self.series[self.series.isnull()],
            target_freq=cur_freq,
        )

    @log_entry
    def resample(
        self,
        target_freq: Union[str, pd.Timedelta],
        shift_tolerance: pd.Timedelta = pd.Timedelta("4min"),
        origin=None,
        origin_simplify_tolerance: pd.Timedelta = pd.Timedelta("4min"),
    ) -> None:
        """
        Resample to a new time resolution.

        All observational records, outliers, and gaps are resampled to a new
        target frequency. Each present timestamp is mapped to a target timestamp,
        present at the timeseries of target_freq, respecting a maximum shift
        set by the shift_tolerance.

        A new origin (start timestamp) can be set by the argument, or it can be
        deduced from the current present origin.

        Parameters
        ----------
        target_freq : str or pandas.Timedelta
            The target frequency to coarsen all records to.
        shift_tolerance : pandas.Timedelta, optional
            The maximum translation (in time) to map a timestamp to a target timestamp.
        origin : datetime.datetime, optional
            Define the origin (first timestamp) for the observations.
        origin_simplify_tolerance : pandas.Timedelta, optional
            Tolerance for origin simplification.

        Warning
        -------
        Since the gaps depend on the record's frequency and origin, all gaps are
        removed and re-located. All progress in gap filling will be lost.

        Note
        ----
        It is technically possible to increase the time resolution. This will
        not result in an information increase; more gaps are created instead.
        """

        target_freq = pd.to_timedelta(target_freq)
        # Create a timestampmatcher
        timestampmatcher = TimestampMatcher(orig_records=self.series)
        timestampmatcher.make_equispaced_timestamps_mapper(
            freq_estimation_method="highest",  # irrelevant
            freq_estimation_simplify_tolerance=pd.Timedelta(0),  # irrelevant
            origin_simplify_tolerance=origin_simplify_tolerance,
            timestamp_tolerance=shift_tolerance,
            force_freq=target_freq,
            force_origin=origin,
        )
        # update all the attributes holding data
        self.series = timestampmatcher.target_records

        # update the outliers (replace the raw timestamps with the new)
        raw_datetime_map = timestampmatcher.get_raw_map()  # mapping from original timestamps to new timestamps
        for outlinfo in self.outliers:
            outlinfo.remap_timestamps(mapping=raw_datetime_map)

        # remap the outliers_values_bin timestamps to match the new equispaced timestamps
        if not self.outliers_values_bin.empty:
            # Drop outlier values whose timestamps are not in the new time grid
            self.outliers_values_bin = self.outliers_values_bin[
                self.outliers_values_bin.index.isin(raw_datetime_map.keys())
            ]
            self.outliers_values_bin.index = self.outliers_values_bin.index.map(
                lambda ts: raw_datetime_map[ts]
            )

        # create gaps
        orig_gapsdf = self.gapsdf
        if bool(self.gaps):
            logger.warning(
                "The present gaps are removed, new gaps are constructed for %s.", self
            )
            self.gaps = []

        # new created-by-resampling missing timestamps
        new_missing = timestampmatcher.gap_records
        # the original gaps timestamp
        orig_missing = orig_gapsdf["value"]
        orig_missing = orig_missing[
            orig_missing.index.isin(self.series.index)
        ]  # drop the records belonging to previous freq that do not exist anymore

        # combine both sets and construct new gaps
        all_missing = pd.concat([new_missing, orig_missing]).sort_index()

        # Construct gaps
        self.gaps = self._find_gaps(
            missingrecords=all_missing,
            target_freq=pd.to_timedelta(timestampmatcher.target_freq),
        )

    @log_entry
    def get_info(self, printout: bool = True) -> Union[str, None]:
        """
        Retrieve and optionally print basic information about the sensor data.

        Parameters
        ----------
        printout : bool, optional
            If True, the information will be printed to the console. If False,
            the information will be returned as a string. Default is True.

        Returns
        -------
        str or None
            If `printout` is False, returns a string containing the information
            about the sensor data. If `printout` is True, returns None.

        """

        infostr = ""
        infostr += printing.print_fmt_title("General info of SensorData")
        infostr += printing.print_fmt_line(
            f"{self.obstype.name} records of {self.stationname}:", 0
        )
        infostr += self._get_info_core(nident_root=1)

        if printout:
            print(infostr)
        else:
            return infostr

    def _get_info_core(self, nident_root=1) -> str:
        """Build a formatted info string with core sensor-data attributes.

        Parameters
        ----------
        nident_root : int, optional
            Base indentation level for printed lines.  Default is 1.

        Returns
        -------
        str
            Formatted string listing observation type, time range, resolution,
            and outlier/gap counts.
        """
        infostr = ""
        infostr += printing.print_fmt_line(
            f"{self.obstype.name} observations in {self.obstype.std_unit}", nident_root
        )
        infostr += printing.print_fmt_line(
            f" from {self.start_datetime} -> {self.end_datetime}", nident_root
        )
        infostr += printing.print_fmt_line(
            f" At a resolution of {self.freq}", nident_root
        )

        # outliers info:
        if self.outliersdf.empty:
            infostr += printing.print_fmt_line("No outliers present.", nident_root)
        else:
            infostr += printing.print_fmt_line(
                f"A total of {self.outliersdf.shape[0]} flagged observations (QC outliers).",
                nident_root,
            )
            infostr += printing.print_fmt_line("label counts: ", nident_root + 1)
            infostr += printing.print_fmt_dict(
                self.outliersdf["label"].value_counts().to_dict(), nident_root + 2
            )

        # gaps info:
        if not self.gaps:
            infostr += printing.print_fmt_line("No gaps present.", nident_root)
        else:
            infostr += printing.print_fmt_line(
                f"{len(self.gaps)} gaps present, a total of {self.gapsdf.shape[0]} missing timestamps.",
                nident_root,
            )
            infostr += printing.print_fmt_line("label counts: ", nident_root + 1)
            infostr += printing.print_fmt_dict(
                self.gapsdf["label"].value_counts().to_dict(), nident_root + 2
            )

        return infostr

    # ------------------------------------------
    #    Specials
    # ------------------------------------------

    @log_entry
    def convert_to_standard_units(self) -> None:
        """
        Convert the data records to the standard units defined in the observation type.
        """

        self.series = self.obstype.convert_to_standard_units(
            input_data=self.series, input_unit=self.obstype.original_unit
        )

    # ------------------------------------------
    #    plots
    # ------------------------------------------

    # plots are defined on station and dataset level

    @copy_doc(sensordata_simple_pd_plot)
    def pd_plot(self, show_labels: list = ["ok"], **pdplotkwargs) -> Axes:
        return sensordata_simple_pd_plot(self, show_labels=show_labels, **pdplotkwargs)

    # ------------------------------------------
    #    Quality Control (technical qc + value-based qc)
    # ------------------------------------------

    

    @log_entry
    def duplicated_timestamp_check(self) -> None:
        """
        Check for duplicated timestamps in the series.

        Raises
        ------
        MetobsQualityControlError
            If the check is already applied.
        """
        qcresult = qc.duplicated_timestamp_check(records=self.series)

        # Drop duplicates from the series, keep only the first occurrence.
        # After this, self.series has a unique index so _update_outliers can
        # safely store values and set them to NaN.
        self.series = self.series[~self.series.index.duplicated(keep="first")]

        # Store QCresult and move flagged values to outliers_values_bin.
        # The timestamps are still raw (not yet remapped to perfect spacing);
        # they will be remapped later in _setup via qcresult.remap_timestamps.
        self._update_outliers(qcresult=qcresult)

    @log_entry
    def gross_value_check(self, **qckwargs) -> None:
        """
        Perform a gross value check on the series.

        Parameters
        ----------
        **qckwargs : dict
            Additional keyword arguments for the check.
        """

        qcresult = qc.gross_value_check(records=self.series, **qckwargs)
        self._update_outliers(qcresult=qcresult)

    @log_entry
    def persistence_check(self, **qckwargs) -> None:
        """
        Perform a persistence check on the series.

        Parameters
        ----------
        **qckwargs : dict
            Additional keyword arguments for the check.
        """

        qcresult = qc.persistence_check(records=self.series, **qckwargs)
        self._update_outliers(qcresult=qcresult)

    @log_entry
    def repetitions_check(self, **qckwargs) -> None:
        """
        Perform a repetitions check on the series.

        Parameters
        ----------
        **qckwargs : dict
            Additional keyword arguments for the check.
        """

        qcresult = qc.repetitions_check(records=self.series, **qckwargs)
        self._update_outliers(qcresult=qcresult)

    @log_entry
    def step_check(self, **qckwargs) -> None:
        """
        Perform a step check on the series.

        Parameters
        ----------
        **qckwargs : dict
            Additional keyword arguments for the check.
        """

        qcresult = qc.step_check(records=self.series, **qckwargs)
        self._update_outliers(qcresult=qcresult)

    @log_entry
    def window_variation_check(self, **qckwargs) -> None:
        """
        Perform a window variation check on the series.

        Parameters
        ----------
        **qckwargs : dict
            Additional keyword arguments for the check.
        """

        qcresult = qc.window_variation_check(records=self.series, **qckwargs)
        self._update_outliers(qcresult=qcresult)

    @log_entry
    def get_qc_freq_statistics(self) -> pd.DataFrame:
        """
        Generate quality control (QC) frequency statistics.

        This method calculates the frequency statistics for various QC checks
        applied, including the number of records labeled as
        'good', 'gap', and outliers for each QC check. The results are returned
        as a pandas DataFrame.

        Returns
        -------
        pandas.DataFrame
            A DataFrame containing the QC frequency statistics. The DataFrame
            has a multi-index with the station name and QC check label, and
            includes the following columns:

            * `N_all`: Total number of records in the dataset (including gaps).
            * `N_labeled`: Number of records with the specific label.
            * `N_checked`: Number of records checked for the specific QC check.
              This is not necessarily the same as `N_all`, as some records may be
              excluded from the check due to previous QC checks.

        """
        
        
        empty_flags = {
            flagged_cond: 0,
            pass_cond: 0,
            unmet_cond: 0,
            saved_cond: 0,
            unchecked_cond: 0}
        
        qcdict = {}
        for qcres in self.outliers:
            qcdict[qcres.checkname] = empty_flags | qcres.flags.value_counts().to_dict()
    
        #Convert to a pandas series with  multiindex ['checkname', 'flag'] and the name is 'counts'

        qcdf = pd.DataFrame.from_dict(qcdict, orient='index')
        qcdf.index.name = 'checkname'
        qcseries = qcdf.stack(future_stack=True)
        qcseries.index = qcseries.index.set_names('flag', level=-1)
        qcseries.name = 'counts'    
        return qcseries
        
        

    # ------------------------------------------
    #    Gaps related
    # ------------------------------------------
    @log_entry
    def fill_gap_with_modeldata(
        self,
        modeltimeseries: ModelTimeSeries,
        method: Literal[
            "raw", "debiased", "diurnal_debiased", "weighted_diurnal_debiased"
        ] = "raw",
        overwrite_fill: bool = False,
        method_kwargs: dict = {},
    ) -> None:
        """
        Fill gaps using model data.

        Parameters
        ----------
        modeltimeseries : pd.Series or similar
            Model data timeseries to use for filling.
        method : str, optional
            Gap filling method, by default "raw".
        overwrite_fill : bool, optional
            Whether to overwrite existing fills, by default False.
        method_kwargs : dict, optional
            Additional keyword arguments for the method, by default {}.

        Raises
        ------
        NotImplementedError
            If the specified method is not implemented.
        """

        for gap in self.gaps:
            if not gap.flag_can_be_filled(
                overwrite_fill
            ):  # if flag_can_be_filled returns False, Gaps won't be filled
                logger.warning(
                    f"{gap} cannot be filled (because it has a fill status {gap.fillstatus}), and overwrite fill is {overwrite_fill}."
                )
                continue
            if overwrite_fill:
                # clear previous fill info
                gap.flush_fill()

            logger.debug(f"Filling {gap} with {method} model data.")

            if method == "raw":
                gap.raw_model_gapfill(modeltimeseries=modeltimeseries, **method_kwargs)
            elif method == "debiased":
                gap.debiased_model_gapfill(
                    sensordata=self,
                    modeltimeseries=modeltimeseries,
                    **method_kwargs,
                )
            elif method == "diurnal_debiased":
                gap.diurnal_debiased_model_gapfill(
                    sensordata=self,
                    modeltimeseries=modeltimeseries,
                    **method_kwargs,
                )
            elif method == "weighted_diurnal_debiased":
                gap.weighted_diurnal_debiased_model_gapfill(
                    sensordata=self,
                    modeltimeseries=modeltimeseries,
                    **method_kwargs,
                )
            else:
                raise NotImplementedError(
                    f"Model data gapfill method: {method} is not implemented!"
                )

    @log_entry
    def interpolate_gaps(
        self,
        method: str = "time",
        max_gap_duration_to_fill: Union[str, pd.Timedelta] = pd.Timedelta("3h"),
        n_leading_anchors: int = 1,
        n_trailing_anchors: int = 1,
        max_lead_to_gap_distance: Union[pd.Timedelta, str, None] = None,
        max_trail_to_gap_distance: Union[pd.Timedelta, str, None] = None,
        method_kwargs: dict = {},
        overwrite_fill: bool = False,
    ) -> None:
        """
        Interpolate gaps in the data.

        Parameters
        ----------

        method : str, optional
            Interpolation method, by default "time".
        max_gap_duration_to_fill : str or pandas.Timedelta, optional
            The maximum gap duration of to fill with interpolation. The result is
            independent on the time-resolution of the gap. Defaults to 3 hours.
        n_leading_anchors : int, optional
            Number of leading anchors, by default 1.
        n_trailing_anchors : int, optional
            Number of trailing anchors, by default 1.
        max_lead_to_gap_distance : optional
            Maximum distance from leading anchor to gap.
        max_trail_to_gap_distance : optional
            Maximum distance from trailing anchor to gap.
        method_kwargs : dict, optional
            Additional keyword arguments for the interpolation method, by default {}.
        overwrite_fill : bool, optional
            Whether to overwrite existing fills, by default False.
        """

        max_gap_duration_to_fill = to_timedelta(max_gap_duration_to_fill)

        for gap in self.gaps:
            if not gap.flag_can_be_filled(overwrite_fill):
                logger.warning(
                    f"{gap} cannot be filled (because it has a fill status {gap.fillstatus}), and overwrite fill is {overwrite_fill}."
                )
                continue

            # clear previous fill info
            gap.flush_fill()

            logger.debug(f"Filling {gap} with {method} interpolation.")
            gap.interpolate(
                sensordata=self,
                method=method,
                max_gap_duration_to_fill=max_gap_duration_to_fill,
                n_leading_anchors=n_leading_anchors,
                n_trailing_anchors=n_trailing_anchors,
                max_lead_to_gap_distance=max_lead_to_gap_distance,
                max_trail_to_gap_distance=max_trail_to_gap_distance,
                method_kwargs=method_kwargs,
            )


def _format_timestamp_index(
    timestamps: np.ndarray,
    input_tz: Union[tzfile | tzinfo | str],
) -> pd.DatetimeIndex:
    """Convert a raw timestamp array to a timezone-aware :class:`pandas.DatetimeIndex`.

    Localises *timestamps* to *input_tz* and then converts to the internal
    storage timezone (``SensorData._target_tz``, default ``UTC``).

    Parameters
    ----------
    timestamps : numpy.ndarray
        Array of timestamp values (any format accepted by
        :func:`timestamps_to_datetimeindex`).
    input_tz : str or tzinfo
        Timezone of the input *timestamps*.

    Returns
    -------
    pandas.DatetimeIndex
        Timezone-aware DatetimeIndex in the internal storage timezone.
    """

    # Localize the timestamps with the input timezone (to timezone aware)
    dt_index = timestamps_to_datetimeindex(
        timestamps=timestamps,
        current_tz=input_tz,
        name="datetime",
    )
    # Convert to the tz to store in (default is UTC)
    return convert_timezone(dt_index, target_tz=SensorData._target_tz)
