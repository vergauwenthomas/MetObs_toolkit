import logging
from typing import Literal

import numpy as np
import pandas as pd
from matplotlib.pyplot import Axes

import metobs_toolkit.settings_files.default_formats_settings as defaults
from metobs_toolkit.timestampmatcher import TimestampMatcher
from metobs_toolkit.obstypes import Obstype
from metobs_toolkit.gap import Gap
import metobs_toolkit.qc_collection as qc
import metobs_toolkit.backend_collection.timeseries_plotting as plotting
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
        timezone: str | pd.Timedelta = "UTC",
        freq_estimation_method: Literal["highest", "median"] = "median",
        freq_estimation_simplify_tolerance: pd.Timedelta | str = pd.Timedelta("1min"),
        origin_simplify_tolerance: pd.Timedelta | str = pd.Timedelta("1min"),
        timestamp_tolerance: pd.Timedelta | str = pd.Timedelta("4min"),
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

        self._freqsettings = {
            "freq_estimation_method": freq_estimation_method,
            "freq_estimation_simplify_tolerance": freq_estimation_simplify_tolerance,
            "orig_simplify_tolerance": origin_simplify_tolerance,
            "timestamp_tolerance": timestamp_tolerance,
        }

        # Setup the SensorData --> apply qc control on import, find gaps, unit conversions etc
        self._setup()
        logger.info("SensorData initialized successfully.")

    def __str__(self) -> str:
        """Return a string representation of the SensorData object."""
        return f"{self.obstype.name} data of station {self.stationname}."

    def _setup(self) -> None:
        """
        Setup the SensorData object.

        This includes:
        1. Invalid check (records that could not be typecast to numeric)
        2. Find the duplicates (remove them from observations + add them to outliers)
        3. Convert the values to standard units + update the observation types
        4. Find gaps in the records (duplicates are excluded from the gaps)
        5. Get a frequency estimate per station
        6. Initiate the gaps (find missing records)
        7. Add the missing records to the dataframe
        """
        logger.debug("Setting up SensorData for %s", self.stationname)

        # invalid check
        self.invalid_value_check()  # must be applied first!

        # remove duplicated timestamps
        self.duplicated_timestamp_check()

        # convert units to standard units
        self.convert_to_standard_units()

        # format to perfect time records
        timestamp_matcher = TimestampMatcher(orig_records=self.series)
        timestamp_matcher._make_equispaced_timestamps_mapper(
            freq_estimation_method=self._freqsettings["freq_estimation_method"],
            freq_estimation_simplify_tolerance=self._freqsettings[
                "freq_estimation_simplify_tolerance"
            ],
            origin_simplify_tolerance=self._freqsettings["orig_simplify_tolerance"],
            timestamp_tolerance=self._freqsettings["timestamp_tolerance"],
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

        self.gaps = self._find_gaps(
            missingrecords=timestamp_matcher.gap_records,
            target_freq=pd.to_timedelta(timestamp_matcher.target_freq),
        )

    def resample(
        self,
        target_freq,
        shift_tolerance=pd.Timedelta("4min"),
        origin=None,
        direction="nearest",
    ):

        target_freq = pd.to_timedelta(target_freq)
        # Create a timestampmatcher
        timestampmatcher = TimestampMatcher(orig_records=self.series)
        timestampmatcher._reindex_from_perfect_to_perfect(
            target_freq=target_freq,
            shift_tolerance=shift_tolerance,
            origin=origin,
            direction=direction,
        )
        # update all the attributes holding data
        self.series = timestampmatcher.target_records

        # update the outliers (replace the raw timestamps with the new)
        outl_datetime_map = timestampmatcher.get_outlier_map()
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

        self.gaps = self._find_gaps(
            missingrecords=timestampmatcher.gap_records,
            target_freq=pd.to_timedelta(timestampmatcher.target_freq),
        )

    def get_info(self, printout: bool = True) -> str | None:

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
        Create a DataFrame from SensorData series.

        Returns
        -------
        pd.DataFrame
            DataFrame with ['datetime', 'obstype'] as index and ['value'] as column.
        """
        logger.debug(
            "Creating DataFrame from SensorData series for %s", self.stationname
        )
        # get all records
        df = (
            self.series.to_frame()
            .rename(columns={self.obstype.name: "value", self.stationname: "value"})
            .assign(label=defaults.label_def["goodrecord"]["label"])
        )

        outliersdf = self.outliersdf[["value", "label"]]
        gapsdf = self.gapsdf[["value", "label"]]

        # concat all together (do not change order)
        to_concat = [df]
        if not outliersdf.empty:
            to_concat.append(outliersdf)
        if not gapsdf.empty:
            to_concat.append(gapsdf)
        df = pd.concat(to_concat)
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
        Format all outliers and their labels in one pandas.DataFrame.

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
            checkdf["label"] = defaults.label_def[checkname]["label"]
            to_concat.append(checkdf)

        totaldf = pd.concat(to_concat)
        totaldf.sort_index(inplace=True)
        logger.debug("Outliers DataFrame created successfully for %s", self.stationname)
        return totaldf

    @property
    def gapsdf(self) -> pd.DataFrame:
        to_concat = []
        if bool(self.gaps):
            for gap in self.gaps:
                gapdf = (
                    gap.gapdf.reset_index()
                    .set_index("datetime")
                    .drop(columns=["name"])
                    .rename(columns={self.obstype.name: "value"})
                    .assign(label="gap")
                )

                to_concat.append(gapdf)
            return pd.concat(to_concat).sort_index()
        else:
            return pd.DataFrame(
                columns=["value", "label"], index=pd.DatetimeIndex([], name="datetime")
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
        # note: sometimes 'h' is returned, and this gives issues, so add a 1 in front
        if not freq[0].isdigit():
            freq = "1" + freq

        return pd.Timedelta(freq)

    def _format_timestamp_index(
        self, timestamps: np.ndarray, tz: str | pd.Timedelta
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
                name=self.stationname,
                startdt=gapgroup.index.min(),
                enddt=gapgroup.index.max(),
                obstype=self.obstype,
                records_freq=pd.to_timedelta(target_freq),
            )
            gaps.append(gap)
        return gaps

    # ------------------------------------------
    #    plots
    # ------------------------------------------

    def make_plot(
        self,
        colorby: Literal["station", "label"] = "label",
        linecolor=None,
        show_outliers=True,
        show_gaps=True,
        ax=None,
        figkwargs: dict = {},
        title: str | None = None,
    ) -> Axes:
        # define figure
        if ax is None:
            ax = plotting.create_axes(**figkwargs)

        if colorby == "station":
            # Define a color
            if linecolor is None:
                # create a new color
                color = plotting.create_station_color_map(["dummy"])["dummy"]
            else:
                color = linecolor

            ax = plotting.plot_timeseries_as_one_color(
                sensordata=self,
                color=color,
                ax=ax,
                show_gaps=show_gaps,
                show_outliers=show_outliers,
            )
        elif colorby == "label":
            if linecolor is not None:
                logger.warning(
                    f"linecolor (={linecolor} is ignore when colorby={colorby})"
                )
            # TODO: forward plotkwargs
            ax = plotting.plot_timeseries_color_by_label(
                sensordata=self, show_gaps=show_gaps, show_outliers=show_outliers, ax=ax
            )

        else:
            raise Exception(
                f'colorby is either "label" or "station", but not {colorby}'
            )

        # Add Styling attributes
        # Set title:
        if title is None:
            plotting.set_title(
                ax, f"{self.obstype.name} data for station {self.stationname}"
            )
        else:
            plotting.set_title(ax, title)
        # Set ylabel
        plotting.set_ylabel(ax, self.obstype._get_plot_y_label())

        # Set xlabel
        plotting.set_xlabel(ax, f"Timestamps (in {self.tz})")

        # Add legend
        plotting.set_legend(ax)
        return ax

    # ------------------------------------------
    #    Quality Control (techincal qc + valuebased qc)
    # ------------------------------------------

    def invalid_value_check(self) -> None:
        """
        Check for invalid values in the series.

        Invalid values are those that could not be cast to numeric.

        Raises
        ------
        MetobsQualityControlError
            If the check is already applied.
        """
        logger.info("Performing invalid value check for %s", self.stationname)

        outlier_timestamps = self.series[~self.series.notnull()]

        self._update_outliers(
            qccheckname="invalid_input",
            outliertimestamps=outlier_timestamps.index,
            check_kwargs={},
            extra_columns={},
            overwrite=False,
        )

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


# ------------------------------------------
#    Helpers
# ------------------------------------------
