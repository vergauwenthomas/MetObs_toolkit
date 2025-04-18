import logging
from typing import Literal, Union

import numpy as np
import pandas as pd
from matplotlib.pyplot import Axes

from metobs_toolkit.backend_collection.df_helpers import save_concat, to_timedelta
from metobs_toolkit.settings_collection import label_def
from metobs_toolkit.timestampmatcher import TimestampMatcher
from metobs_toolkit.obstypes import Obstype
from metobs_toolkit.newgap import Gap
import metobs_toolkit.qc_collection as qc
import metobs_toolkit.plot_collection as plotting
from metobs_toolkit.backend_collection.errorclasses import *

logger = logging.getLogger(__file__)


class SensorData:
    """Holds data for one station for one sensor."""

    def __init__(
        self,
        stationname: str,
        datarecords: np.ndarray,
        timestamps: np.ndarray,
        obstype: Obstype,
        datadtype=np.float32,
        timezone: Union[str, pd.Timedelta] = "UTC",
        freq_estimation_method: Literal["highest", "median"] = "median",
        freq_estimation_simplify_tolerance: Union[pd.Timedelta, str] = pd.Timedelta(
            "1min"
        ),
        origin_simplify_tolerance: Union[pd.Timedelta, str] = pd.Timedelta("1min"),
        timestamp_tolerance: Union[pd.Timedelta, str] = pd.Timedelta("4min"),
    ):
        """
        Initialize SensorData.

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
        timezone : str or pd.Timedelta, optional
            Timezone of the timestamps, by default 'UTC'.
        freq_estimation_method : {'highest', 'median'}, optional
            Method to estimate frequency, by default 'median'.
        freq_estimation_simplify_tolerance : pd.Timedelta or str, optional
            Tolerance for frequency estimation simplification, by default pd.Timedelta('1min').
        origin_simplify_tolerance : pd.Timedelta or str, optional
            Tolerance for origin simplification, by default pd.Timedelta('1min').
        timestamp_tolerance : pd.Timedelta or str, optional
            Tolerance for timestamp matching, by default pd.Timedelta('4min').
        """
        logger.info("Initializing SensorData with obstype: %s", obstype.name)

        if not isinstance(stationname, str):
            raise TypeError("stationname must be a string")
        if not isinstance(datarecords, np.ndarray):
            raise TypeError("datarecords must be a numpy array")
        if not isinstance(timestamps, np.ndarray):
            raise TypeError("timestamps must be a numpy array")
        if not isinstance(obstype, Obstype):
            raise TypeError("obstype must be an instance of Obstype")

        # Set data
        self._stationname = stationname
        self.obstype = obstype
        data = pd.Series(
            data=pd.to_numeric(datarecords, errors="coerce").astype(datadtype),
            index=self._format_timestamp_index(timestamps, timezone),
            name=obstype.name,
        )

        data.index.name = "datetime"
        self.series = data  # datetime as index

        # outliers
        self.outliers = []  # List of {'checkname': ..., 'df': ....., 'settings': }

        # gaps
        self.gaps = []  # list of Gap's

        # Setup the SensorData --> apply qc control on import, find gaps, unit conversions etc
        self._setup(
            freq_estimation_method=freq_estimation_method,
            freq_estimation_simplify_tolerance=freq_estimation_simplify_tolerance,
            origin_simplify_tolerance=origin_simplify_tolerance,
            timestamp_tolerance=timestamp_tolerance,
            apply_invalid_check=True,
            apply_dupl_check=True,
            apply_unit_conv=True,
        )

        logger.info("SensorData initialized successfully.")

    def __eq__(self, other):
        if not isinstance(other, SensorData):
            return False
        return (
            self.stationname == other.stationname
            and self.df.equals(other.df)  # the df contains outliers and gaps as well
            and self.obstype == other.obstype
        )

    def __repr__(self):
        return f"Sensordata instance of {self.obstype.name} -> {self.stationname}"

    def __str__(self) -> str:
        """Return a string representation of the SensorData object."""
        return f"{self.obstype.name} data of station {self.stationname}."

    def _setup(
        self,
        freq_estimation_method,
        freq_estimation_simplify_tolerance,
        origin_simplify_tolerance,
        timestamp_tolerance,
        apply_invalid_check=True,
        apply_dupl_check=True,
        apply_unit_conv=True,
        force_origin=None,
        force_freq=None,
        force_closing=None,
    ) -> None:
        """
        Setup the SensorData object.

        This includes:

        1. Find the duplicates (remove them from observations + add them to outliers)
        2. Invalid check (records that could not be typecast to numeric) --> are interpreted as gaps.
        3. Convert the values to standard units + update the observation types
        4. Find gaps in the records (duplicates are excluded from the gaps)
        5. Get a frequency estimate per station
        6. Initiate the gaps (find missing records)
        7. Add the missing records to the dataframe
        """
        logger.debug("Setting up SensorData for %s", self.stationname)

        if apply_dupl_check:
            # remove duplicated timestamps
            self.duplicated_timestamp_check()

        if apply_invalid_check:
            # invalid check
            self.invalid_value_check(
                skip_records=self.outliers[0]["df"].index
            )  # skip the records already labeld as duplicates

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
        outl_datetime_map = timestamp_matcher.get_outlier_map()
        for outlinfo in self.outliers:
            outlinfo["df"]["new_datetime"] = outlinfo["df"].index.map(outl_datetime_map)
            outlinfo["df"] = (
                outlinfo["df"]
                .reset_index()
                .rename(
                    columns={"datetime": "raw_timestamp", "new_datetime": "datetime"}
                )
                .set_index("datetime")
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

    def convert_outliers_to_gaps(self):
        """Convert all outliers to gaps.

        This method will convert all outliers to gaps. Doing so new gaps are constructed.


        Returns
        -------
        None.

        Warning
        ------
        All progress on present gaps is ereased, since new gaps are constructed.
        Information on the value and QC flag of the outliers will be lost.
        """

        cur_freq = self.freq
        # Create holes for all the outliers timestamps
        self.series.loc[self.outliersdf.index] = np.nan

        # Flush the outliers
        logger.warning(f"Outliers are flushed of {self}!")
        self.outliers = []

        # Flush the gaps
        if bool(self.gaps):
            logger.warning(f"Flushing current gaps of {self}")
        self.gaps = []

        # Finding new gaps
        self.gaps = self._find_gaps(
            missingrecords=self.series[self.series.isnull()],
            target_freq=cur_freq,
        )

    def resample(
        self,
        target_freq,
        shift_tolerance=pd.Timedelta("4min"),
        origin=None,
        origin_simplify_tolerance=pd.Timedelta("4min"),
    ):
        """Resample to a new time resolution.

        All observational records, outliers and gaps are resampled to a new
        target frequency. Each present timestamp is mapped to a target timestamp,
        present at the timeseries of target_freq, respecting a maximum shift
        set by the shift_tolerance.

        A new origin (=start timestamp) can be set by the argument, or it can be
        deduced from the current present origin.



        Parameters
        ----------
        origin : datetime.datetime, optional
            Define the origin (first timestamp) for the observations. The origin
            is timezone naive and is assumed to have the same timezone as the
            obervations. If None, the earliest occurring timestamp is used as
            origin. The default is None.
        freq : DateOffset, Timedelta or str, optional
            The target frequency to coarsen all records to.
            Ex: '15min' is 15 minutes, '1h', is one hour. The default is '60min'.
        direction : 'backward', 'forward', or 'nearest'
            Whether to search for prior, subsequent, or closest matches for
            mapping to ideal timestamps. The default is 'nearest'.
        timestamp_shift_tolerance : Timedelta or str
            The tolerance string or object representing the maximum translation
            (in time) to map a timestamp to a target timestamp.
            Ex: '5min' is 5 minutes. The default is '4min'.


        Returns
        -------
        None.


        Warning
        ---------
        Since the gaps depend on the record's frequency and origin, all gaps are
        removed and re-located. All progress in gap(filling) will be lost.

        Note
        -------
        It is technically possible to increase the time resolution. This will
        however not result in an information increase, more gaps are
        created instead.

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
        outl_datetime_map = timestampmatcher.get_outlier_map()
        for outlinfo in self.outliers:
            # add mapped timestamps
            outlinfo["df"]["new_datetime"] = outlinfo["df"].index.map(outl_datetime_map)
            # reformat the dataframe
            outlinfo["df"] = (
                outlinfo["df"]
                .reset_index()
                .rename(
                    columns={"datetime": "raw_timestamp", "new_datetime": "datetime"}
                )
                .set_index("datetime")
            )
            # Drop references to NaT datetimes (when qc is applied before resampling)
            outlinfo["df"] = outlinfo["df"].loc[outlinfo["df"].index.notnull()]

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
        ]  # drop the records belonging to previous freq that does not exists anymore

        # combine both sets and construct new gaps
        all_missing = pd.concat([new_missing, orig_missing]).sort_index()

        # Construct gaps
        self.gaps = self._find_gaps(
            missingrecords=all_missing,
            target_freq=pd.to_timedelta(timestampmatcher.target_freq),
        )

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

        Notes
        -----
        The information includes:
        - Observation type and station name.
        - Start and end datetime of the records.
        - Assumed frequency of the data.
        - Number of records and the count of outliers.
        - Number of gaps in the data.
        """

        infostr = ""
        infostr += f"{self.obstype.name} records of {self.stationname}:\n"
        infostr += f"  * from {self.start_datetime} --> {self.end_datetime}\n"
        infostr += f"  * assumed frequency: {self.freq}\n"
        infostr += f"  * Number of records: {self.series.shape[0]}\n"
        infostr += f"    - of which outliers: {self.outliersdf.shape[0]}\n"
        infostr += f"  * Number of gaps: {len(self.gaps)}\n"

        if printout:
            print(infostr)
        else:
            return infostr

    @property
    def df(self) -> pd.DataFrame:
        """
        Create a DataFrame of the sensor records.

        Returns
        -------
        pd.DataFrame
            DataFrame with ['datetime', 'obstype'] as index and ['value', 'label'] as columns.
        """
        logger.debug(
            "Creating DataFrame from SensorData series for %s", self.stationname
        )
        # get all records
        df = (
            self.series.to_frame()
            .rename(columns={self.obstype.name: "value", self.stationname: "value"})
            .assign(label=label_def["goodrecord"]["label"])
        )

        outliersdf = self.outliersdf[["value", "label"]]

        gapsdf = self.gapsdf[["value", "label"]]

        # concat all together (do not change order)
        to_concat = [df]
        if not outliersdf.empty:
            to_concat.append(outliersdf)
        if not gapsdf.empty:
            to_concat.append(gapsdf)
        df = save_concat((to_concat))
        # remove duplicates
        df = df[~df.index.duplicated(keep="last")].sort_index()

        # add 'obstype' as index
        df = (
            df.assign(obstype=self.obstype.name)
            .reset_index()
            .set_index(["datetime", "obstype"])
        )
        return df

    @property
    def outliersdf(self) -> pd.DataFrame:
        """
        Create a DataFrame of the outlier-records.

        Returns
        -------
        pd.DataFrame
            DataFrame containing all outliers and their labels.
        """
        logger.debug("Creating outliers DataFrame for %s", self.stationname)
        to_concat = []
        for outlierinfo in self.outliers:
            checkname = outlierinfo["checkname"]
            checkdf = outlierinfo["df"]
            checkdf["label"] = label_def[checkname]["label"]
            to_concat.append(checkdf)

        totaldf = save_concat(to_concat)

        if totaldf.empty:
            # return empty dataframe
            totaldf = pd.DataFrame(
                columns=["value", "label"], index=pd.DatetimeIndex([], name="datetime")
            )
        else:
            totaldf.sort_index(inplace=True)

        logger.debug("Outliers DataFrame created successfully for %s", self.stationname)
        return totaldf

    @property
    def gapsdf(self) -> pd.DataFrame:
        """
        Create a DataFrame of the gap-records. Gap-records are all the
        timestamps that make up a gap.

        Returns
        -------
        pd.DataFrame
            DataFrame containing all gap-records and their labels.
        """

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

    @property
    def stationname(self) -> str:
        """
        Return the name of the station this SensorData belongs to.

        Returns
        -------
        str
            station name
        """
        return self._stationname

    @property
    def tz(self):
        """
        Return the (up-to-date) timezone of the stored timestamps.

        Returns
        -------
        str
            Timezone of the stored timestamps.
        """

        return self.series.index.tz

    @property
    def start_datetime(self) -> pd.Timestamp:
        """
        Return the start datetime of the series.

        Returns
        -------
        pd.Timestamp
            Start datetime of the series.
        """
        return self.series.index.min()

    @property
    def end_datetime(self) -> pd.Timestamp:
        """
        Return the end datetime of the series.

        Returns
        -------
        pd.Timestamp
            End datetime of the series.
        """
        return self.series.index.max()

    @property
    def freq(self) -> pd.Timedelta:
        """
        Return the frequency of the series.

        Returns
        -------
        pd.Timedelta
            Frequency of the series.
        """
        freq = pd.infer_freq(self.series.index)
        if freq is None:
            raise ValueError("Frequency could not be computed.")
        return to_timedelta(freq)

    def _format_timestamp_index(
        self, timestamps: np.ndarray, tz: Union[str, pd.Timedelta]
    ) -> pd.DatetimeIndex:
        """
        Format the timestamp index.

        Parameters
        ----------
        timestamps : np.ndarray
            Array of timestamps.
        tz : str or pd.Timedelta
            Timezone of the timestamps.

        Returns
        -------
        pd.DatetimeIndex
            Formatted timestamp index.
        """
        logger.debug("Formatting timestamp index for %s", self.stationname)
        if not isinstance(timestamps, np.ndarray):
            raise TypeError("timestamps must be a numpy array")
        if not isinstance(tz, (str, pd.Timedelta)):
            raise TypeError("tz must be a string or pandas Timedelta")
        return pd.DatetimeIndex(data=timestamps, tz=tz)

    def _update_outliers(
        self,
        qccheckname: str,
        outliertimestamps: pd.DatetimeIndex,
        check_kwargs: dict,
        extra_columns: dict = {},
        overwrite: bool = False,
    ) -> None:
        """
        Update the outliers attribute.

        Parameters
        ----------
        qccheckname : str
            Name of the quality control check.
        outliertimestamps : pd.DatetimeIndex
            Datetime index of the outliers.
        check_kwargs : dict
            Additional arguments for the check.
        extra_columns : dict, optional
            Extra columns to add to the outliers DataFrame, by default {}.
        overwrite : bool, optional
            Whether to overwrite existing outliers, by default False.
        """
        logger.debug(
            "Updating outliers for %s with check %s", self.stationname, qccheckname
        )

        if not isinstance(qccheckname, str):
            raise TypeError("qccheckname must be a string")
        if not isinstance(outliertimestamps, pd.DatetimeIndex):
            raise TypeError("outliertimestamps must be a pandas DatetimeIndex")
        if not isinstance(check_kwargs, dict):
            raise TypeError("check_kwargs must be a dictionary")
        if not isinstance(extra_columns, dict):
            raise TypeError("extra_columns must be a dictionary")
        if not isinstance(overwrite, bool):
            raise TypeError("overwrite must be a boolean")

        for applied_qc_info in self.outliers:
            if qccheckname == applied_qc_info.keys():
                if overwrite:
                    self.outliers.remove(applied_qc_info)
                else:
                    raise MetobsQualityControlError(
                        f"The {qccheckname} is already applied on {self}. Fix error or set overwrite=True"
                    )

        outlier_values = self.series.loc[outliertimestamps]
        outlier_values = outlier_values[~outlier_values.index.duplicated(keep="first")]

        datadict = {"value": outlier_values.to_numpy()}
        datadict.update(extra_columns)
        df = pd.DataFrame(data=datadict, index=outlier_values.index)

        self.outliers.append(
            {"checkname": qccheckname, "df": df, "settings": check_kwargs}
        )

        self.series.loc[outliertimestamps] = np.nan

    def _find_gaps(self, missingrecords: pd.Series, target_freq: pd.Timedelta) -> list:
        """
        Identify gaps in the missing records based on the target frequency.

        Parameters
        ----------
        missingrecords : pd.Series
            A pandas Series containing the missing records with datetime index.
        target_freq : pd.Timedelta
            The target frequency to identify gaps.

        Returns
        -------
        list
            A list of Gap objects representing the identified gaps.
        """
        logger.debug("Finding gaps for %s", self.stationname)

        if not isinstance(missingrecords, pd.Series):
            raise TypeError("missingrecords must be a pandas Series")
        if not isinstance(target_freq, pd.Timedelta):
            raise TypeError("target_freq must be a pandas Timedelta")

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

    def _rename(self, trgname: str):
        self._stationname = str(trgname)
        for gap in self.gaps:
            gap.name = str(trgname)

    # ------------------------------------------
    #    Specials
    # ------------------------------------------

    def convert_to_standard_units(self) -> None:
        """
        Convert the data records to the standard units defined in the observation type.
        """
        logger.info(
            "Converting data records to standard units for %s", self.stationname
        )

        self.series = self.obstype.convert_to_standard_units(
            input_data=self.series, input_unit=self.obstype.original_unit
        )

    # ------------------------------------------
    #    plots
    # ------------------------------------------

    # plots are defined on station and dataset level

    # ------------------------------------------
    #    Quality Control (techincal qc + valuebased qc)
    # ------------------------------------------

    def invalid_value_check(self, skip_records: pd.DatetimeIndex) -> None:
        """
        Check for invalid values in the series.

        Invalid values are those that could not be cast to numeric.

        Raises
        ------
        MetobsQualityControlError
            If the check is already applied.
        """

        logger.info("Performing invalid value check for %s", self.stationname)

        skipped_data = self.series.loc[skip_records]
        targets = self.series.drop(skip_records)
        # Option 1: Create a outlier label for these invalid inputs,
        # and treath them as outliers
        # outlier_timestamps = targets[~targets.notnull()]

        # self._update_outliers(
        #     qccheckname="invalid_input",
        #     outliertimestamps=outlier_timestamps.index,
        #     check_kwargs={},
        #     extra_columns={},
        #     overwrite=False,
        # )

        # Option 2: Since there is not numeric value present, these timestamps are
        # interpreted as gaps --> remove the timestamp, so that it is captured by the
        # gap finder.

        # Note: do not treat the first/last timestamps differently. That is
        # a philosiphycal choice.

        self.series = targets[targets.notnull()]  # subset to numerical casted values
        # add the skipped records back
        self.series = pd.concat([self.series, skipped_data]).sort_index()

    def duplicated_timestamp_check(self) -> None:
        """
        Check for duplicated timestamps in the series.

        Raises
        ------
        MetobsQualityControlError
            If the check is already applied.
        """
        logger.info("Performing duplicated timestamp check for %s", self.stationname)

        duplicates = pd.Series(
            data=self.series.index.duplicated(keep=False), index=self.series.index
        )
        duplicates = duplicates.loc[duplicates]
        duplicates = duplicates[duplicates.index.duplicated(keep="first")]

        self._update_outliers(
            qccheckname="duplicated_timestamp",
            outliertimestamps=duplicates.index,
            check_kwargs={},
            extra_columns={},
            overwrite=False,
        )

        self.series = self.series[~self.series.index.duplicated(keep="first")]

    def gross_value_check(self, **qckwargs):
        outlier_timestamps = qc.gross_value_check(records=self.series, **qckwargs)
        self._update_outliers(
            qccheckname="gross_value",
            outliertimestamps=outlier_timestamps,
            check_kwargs={**qckwargs},
            extra_columns={},
            overwrite=False,
        )

    def persistence_check(self, **qckwargs):

        outlier_timestamps = qc.persistence_check(records=self.series, **qckwargs)

        self._update_outliers(
            qccheckname="persistence",
            outliertimestamps=outlier_timestamps,
            check_kwargs={**qckwargs},
            extra_columns={},
            overwrite=False,
        )

    def repetitions_check(self, **qckwargs):
        outlier_timestamps = qc.repetitions_check(records=self.series, **qckwargs)

        self._update_outliers(
            qccheckname="repetitions",
            outliertimestamps=outlier_timestamps,
            check_kwargs={**qckwargs},
            extra_columns={},
            overwrite=False,
        )

    def step_check(self, **qckwargs):
        outlier_timestamps = qc.step_check(records=self.series, **qckwargs)

        self._update_outliers(
            qccheckname="step",
            outliertimestamps=outlier_timestamps,
            check_kwargs={**qckwargs},
            extra_columns={},
            overwrite=False,
        )

    def window_variation_check(self, **qckwargs):
        outlier_timestamps = qc.window_variation_check(records=self.series, **qckwargs)

        self._update_outliers(
            qccheckname="window_variation",
            outliertimestamps=outlier_timestamps,
            check_kwargs={**qckwargs},
            extra_columns={},
            overwrite=False,
        )

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
            - `N_all`: Total number of records in the dataset (including gaps).
            - `N_labeled`: Number of records with the specific label.
            - `N_checked`: Number of records checked for the specific QC check.
              This is not nesesarily the same as `N_all`, as some records may be
              excluded from the check due to previous QC checks.

        """
        infodict = {}  # checkname : details
        ntotal = self.series.shape[0]  # gaps included !!
        already_rejected = self.gapsdf.shape[0]  # initial gap records
        # add the 'ok' labels
        infodict[label_def["goodrecord"]["label"]] = {
            "N_all": ntotal,
            "N_labeled": self.series[self.series.notnull()].shape[0],
        }
        # add the 'gap' labels

        infodict[label_def["regular_gap"]["label"]] = {
            "N_all": ntotal,
            "N_labeled": already_rejected,
        }

        # add the qc check labels
        for check in self.outliers:
            n_outliers = check["df"].shape[0]
            n_checked = ntotal - already_rejected
            outlierlabel = label_def[check["checkname"]]["label"]
            infodict[outlierlabel] = {
                "N_labeled": n_outliers,
                "N_checked": n_checked,
                "N_all": ntotal,
            }

            # remove the outliers of the previous check
            already_rejected = already_rejected + n_outliers

        # Convert to a dataframe
        checkdf = pd.DataFrame(infodict).transpose()
        checkdf.index.name = "qc_check"
        checkdf["name"] = self.stationname
        checkdf = checkdf.reset_index().set_index(["name", "qc_check"])

        return checkdf

    # ------------------------------------------
    #    Gaps related
    # ------------------------------------------
    def fill_gap_with_modeldata(
        self, modeltimeseries, method="raw", overwrite_fill=False, method_kwargs={}
    ):
        for gap in self.gaps:
            if not gap.flag_can_be_filled(
                overwrite_fill
            ):  # if flag_can_be_filled returns False, Gaps won't be filled
                logger.warning(
                    f"{gap} cannot be filled because it already contains filled values, and overwrite fill is {overwrite_fill}."
                )
                continue
            if overwrite_fill:
                # clear previous fill info
                gap.flush_fill()

            logger.debug(f"filling {gap} with {method} modeldata.")

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
                    f"modeldata gapfill method: {method} is not implemented!"
                )

    def interpolate_gaps(
        self,
        overwrite_fill=False,
        method="time",
        max_consec_fill=10,
        n_leading_anchors=1,
        n_trailing_anchors=1,
        max_lead_to_gap_distance=None,
        max_trail_to_gap_distance=None,
        method_kwargs={},
    ):

        for gap in self.gaps:
            if not gap.flag_can_be_filled(overwrite_fill):
                logger.warning(
                    f"{gap} cannot be filled because it already contains filled values, and overwrite fill is {overwrite_fill}."
                )
                continue

            # clear previous fill info
            gap.flush_fill()

            logger.debug(f"filling {gap} with {method} interpolation.")
            gap.interpolate(
                sensordata=self,
                method=method,
                max_consec_fill=max_consec_fill,
                n_leading_anchors=n_leading_anchors,
                n_trailing_anchors=n_trailing_anchors,
                max_lead_to_gap_distance=max_lead_to_gap_distance,
                max_trail_to_gap_distance=max_trail_to_gap_distance,
                method_kwargs=method_kwargs,
            )
