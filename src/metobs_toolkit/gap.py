import logging
from typing import Union
import numpy as np
import pandas as pd

from metobs_toolkit.obstypes import Obstype
from metobs_toolkit.modeltimeseries import ModelTimeSeries
from metobs_toolkit.backend_collection.argumentcheckers import (
    fmt_timedelta_arg,
)

from metobs_toolkit.backend_collection.df_helpers import convert_to_numeric_series
import metobs_toolkit.backend_collection.printing_collection as printing
from metobs_toolkit.settings_collection import label_def
import metobs_toolkit.gf_collection.gf_common_methods as gf_methods
from metobs_toolkit.gf_collection.debias_gapfill import fill_regular_debias
from metobs_toolkit.gf_collection.diurnal_debias_gapfill import (
    fill_with_diurnal_debias,
    fill_with_weighted_diurnal_debias,
)

from metobs_toolkit.backend_collection.loggingmodule import log_entry
from metobs_toolkit.backend_collection.dev_collection import copy_doc
from metobs_toolkit.backend_collection.dataframe_constructors import gap_df

logger = logging.getLogger("<metobs_toolkit>")

_unfilled_label = "unfilled"
_failed_label = "failed gapfill"
_successful_label = "successful gapfill"


class Gap:
    """
    Represents a gap in observational data for a specific station and observation type.

    Parameters
    ----------
    gaprecords : pd.DatetimeIndex
        The datetime index representing the gap records.
    obstype : Obstype
        The type of observation (e.g., temperature, humidity).
    stationname : str
        The name of the station where the gap occurred.
    """

    def __init__(
        self,
        gaprecords: pd.DatetimeIndex,
        obstype: Obstype,
        stationname: str,
    ):
        """Initialize a Gap object."""
        gaprecords.name = "datetime"
        self._records = pd.Series(data=np.nan, index=gaprecords, name="value")
        self._labels = pd.Series(
            data=label_def["regular_gap"]["label"], index=gaprecords, name="label"
        )
        self._extra_info = pd.Series(
            data="no details", index=gaprecords, name="details"
        )
        self._fillkwargs = {}
        self._obstype = obstype
        self._stationname = stationname

    def __repr__(self):
        """Instance representation."""
        return f"{type(self).__name__}(station={self.stationname}, obstype={self.obstype.name}, start={self.start_datetime}, end={self.end_datetime}, status={self.fillstatus})"

    @property
    def records(self) -> pd.Series:
        """Return the records of the gap."""

        return convert_to_numeric_series(self._records, datadtype=np.float32)

    @property
    def obstype(self) -> Obstype:
        """Return the observation type."""
        return self._obstype

    @property
    def stationname(self) -> str:
        """Return the station name."""
        return self._stationname

    @property
    def fillsettings(self) -> dict:
        """
        Return the settings used for filling the gap.

        The settings are the kwargs (keyword arguments) used in the gapfill methods.

        Returns
        -------
        dict
            A dictionary containing the settings used for filling the gap.
        """
        return self._fillkwargs

    @property
    def fillstatus(self) -> str:
        """
        Returns the fill status of the gap.

        Returns
        -------
        str
            The fill status, which can be one of the following:

            * 'unfilled'
            * 'failed gapfill'
            * 'successful gapfill'

        """
        if self.records.isna().all() and not bool(self._fillkwargs):
            return _unfilled_label
        elif self.records.isna().all() and bool(self._fillkwargs):
            return _failed_label
        elif not self.records.isna().any() and bool(self._fillkwargs):
            return _successful_label
        else:
            raise NotImplementedError(
                "This situation is unforeseen! Please notify developers."
            )

    @property
    def start_datetime(self) -> pd.Timestamp:
        """Return the start datetime of the gap."""
        return min(self.records.index)

    @property
    def end_datetime(self) -> pd.Timestamp:
        """Return the end datetime of the gap."""
        return max(self.records.index)

    @copy_doc(gap_df)
    @property
    def df(self) -> pd.DataFrame:
        return gap_df(self)

    # ------------------------------------------
    #    Get info methods
    # ------------------------------------------

    @log_entry
    def flag_can_be_filled(self, overwrite: bool = False) -> bool:
        """
        Determine if the gap can be filled.

        By default, a gap can be filled if it is not already filled or if the previous gapfill method failed for the gap.
        A gap that is already filled can only be updated if the overwrite flag is set to True.

        Parameters
        ----------
        overwrite : bool, optional
            If True, allows filling regardless of the current fill status. Default is False.

        Returns
        -------
        bool
            True if the gap can be filled, False otherwise.
        """
        logger.debug(
            f"Entering flag_can_be_filled for {self} with overwrite={overwrite}"
        )
        if not isinstance(overwrite, bool):
            raise TypeError("Argument 'overwrite' must be of type bool.")

        if overwrite:
            return True
        if self.fillstatus == _successful_label:
            return False
        else:
            return True

    @log_entry
    def get_info(self, printout: bool = True) -> Union[str, None]:
        """
        Print or return detailed information about the Gap.

        Parameters
        ----------
        printout : bool, optional
            If True, prints the information. If False, returns the information as a string. Default is True.

        Returns
        -------
        str or None
            The gap information as a string if printout is False, otherwise None.
        """

        infostr = ""
        infostr += printing.print_fmt_title("General info of Gap")
        infostr += printing.print_fmt_section("Gap details")

        infostr += printing.print_fmt_line(
            f"Gap of {self.obstype.name} for station: {self.stationname}", 0
        )
        infostr += printing.print_fmt_line(
            f"From {self.start_datetime} -> {self.end_datetime}", 1
        )
        infostr += printing.print_fmt_line(
            f"Duration gap: {self.end_datetime - self.start_datetime}", 1
        )

        infostr += printing.print_fmt_section("Gap filling details")
        infostr += printing.print_fmt_line(f"Gap status: {self.fillstatus}")
        infostr += printing.print_fmt_line("Gapfill settings used:")
        infostr += printing.print_fmt_dict(d=self.fillsettings, identlvl=2)

        if printout:
            print(infostr)
        else:
            return infostr

    @log_entry
    def debiased_model_gapfill(
        self,
        sensordata: "SensorData",  # type: ignore #noqa: F821
        modeltimeseries: ModelTimeSeries,
        leading_period_duration: Union[str, pd.Timedelta],
        min_leading_records_total: int,
        trailing_period_duration: Union[str, pd.Timedelta],
        min_trailing_records_total: int,
    ) -> None:
        """
        Fill the gaps using model data corrected for the bias.

        This method fills the gap using model data corrected for bias. The bias is estimated using a leading (before the gap)
        and trailing (after the gap) period. The bias is computed by combining the leading and trailing period, and comparing
        the model with the observations (not labeled as outliers). The model data is then interpolated to the missing
        records, and corrected with the estimated bias.

        Parameters
        ----------
        sensordata : SensorData
            The corresponding SensorData used in the computation of the bias. Only
            the observations that are not labeled as outliers are used to compute the bias.
        modeltimeseries : ModelTimeSeries
            The model time series used to fill the gap records. The model data
            must be compatible (equivalent obstype and related to the same Station as the gap.)
        leading_period_duration : str or pd.Timedelta
            The duration of the leading period.
        min_leading_records_total : int
            The minimum number of records required in the leading period.
        trailing_period_duration : str or pd.Timedelta
            The duration of the trailing period.
        min_trailing_records_total : int
            The minimum number of records required in the trailing period.

        Returns
        ----------
        None.

        Notes
        -----
        A schematic description of the debiased modeldata gap fill:

        1. Check the compatibility of the `ModelTimeSeries` with the `gap`.
        2. Construct a leading and trailing sample, and test if they meet the required conditions.
        3. Compute the bias of the modeldata (combine leading and trailing samples).
        4. Fill the gap records by using raw (interpolated) modeldata that is corrected by subtracting the bias.
        5. Update the `gap` attributes with the interpolated values, labels, and details.

        """

        leading_period_duration = fmt_timedelta_arg(leading_period_duration)
        trailing_period_duration = fmt_timedelta_arg(trailing_period_duration)

        self._fillkwargs = {
            "applied_gapfill_method": "debias_model_gapfill",
            "leading_period_duration": leading_period_duration,
            "min_leading_records_total": min_leading_records_total,
            "trailing_period_duration": trailing_period_duration,
            "min_trailing_records_total": min_trailing_records_total,
        }

        # 1. Check validity of modeltimeseries
        is_compat, err_msg = gf_methods.check_if_modeltimeseries_is_compatible(
            gap=self,
            modeltimeseries=modeltimeseries,
            lp_duration=leading_period_duration,
            tp_duration=trailing_period_duration,
        )
        if not is_compat:
            self._labels[:] = label_def["failed_debias_modeldata_fill"]["label"]
            self._extra_info[:] = err_msg
            logger.warning(
                f"Incompatible modeldata for debias_model_gapfill: \n{err_msg}"
            )
            return

        # 2. Construct and validity-test leading and trailing periods
        (
            lead_period,
            trail_period,
            continueflag,
        ) = self._setup_lead_and_trail_for_debias_gapfill(
            sensordata=sensordata,
            fail_label=label_def["failed_debias_modeldata_fill"]["label"],
            leading_period_duration=leading_period_duration,
            min_leading_records_total=min_leading_records_total,
            trailing_period_duration=trailing_period_duration,
            min_trailing_records_total=min_trailing_records_total,
        )
        if not continueflag:
            # warnings and gap attributes are already updated
            return

        # 3. Fill the gap
        combdf = gf_methods.create_a_combined_df(
            leadseries=lead_period, trailseries=trail_period, gap=self
        )
        # add modeldata to combdf
        combdf = gf_methods.add_modeldata_to_combdf(
            combineddf=combdf, modeltimeseries=modeltimeseries
        )
        # Fill the missing records
        filleddf = fill_regular_debias(df=combdf)
        filleddf = filleddf.loc[self.records.index]  # subset to gap records

        # 4. Update attributes
        self._records = filleddf["fillvalue"].rename(
            "value"
        )  # set the new filled records

        # set labels
        self._labels.loc[self.records.notna()] = label_def["debias_modeldata_fill"][
            "label"
        ]
        self._labels.loc[self.records.isna()] = label_def[
            "failed_debias_modeldata_fill"
        ]["label"]

        # update details
        self._extra_info = filleddf["msg"].rename("details")

    @log_entry
    def diurnal_debiased_model_gapfill(
        self,
        sensordata: "SensorData",  # type: ignore #noqa: F821
        modeltimeseries: ModelTimeSeries,
        leading_period_duration: pd.Timedelta,
        trailing_period_duration: pd.Timedelta,
        min_debias_sample_size: int,
    ) -> None:
        """
        Fill the gaps using model data corrected for the diurnal bias.

        This method fills the gap using model data corrected for its diurnal bias.
        The diurnal bias is a bias that is estimated for each timestamp in the leading
        and trailing period. All biases are averaged over hour, minute and second, to
        obtain a diurnal bias (for each timestamp).

        Parameters
        ----------
        sensordata : SensorData
            The corresponding SensorData used in the computation of the bias. Only
            the observations that are not labeled as outliers are used to compute the bias.
        modeltimeseries : ModelTimeSeries
            The model time series used to fill the gap records. The model data
            must be compatible (equivalent obstype and related to the same Station as the gap.)
        leading_period_duration : pd.Timedelta
            The duration of the leading period. That is the period before the gap, used
            for bias estimation.
        trailing_period_duration : pd.Timedelta
            The duration of the trailing period. That is the period after the gap, used
            for bias estimation.
        min_debias_sample_size : int
            The minimum number of samples required for bias estimation. If this condition is not met, the gap
            is not filled.

        Returns
        ---------
        None.

        Notes
        -----
        A schematic description of the diurnal debiased modeldata gap fill:

        #. Check the compatibility of the `ModelTimeSeries` with the `gap`.
        #. Construct a leading and trailing sample, and test if they meet the required conditions.
           The required conditions are tested by testing the samplesizes per hour, minute and second for the leading + trailing periods.
        #. A diurnal bias is computed by grouping to hour, minute and second, and averaging the biases.
        #. Fill the gap records by using raw (interpolated) modeldata that is corrected by subtracting the coresponding diurnal bias.
        #. Update the `gap` attributes with the interpolated values, labels, and details.

        Notes
        --------
        Note that a suitable `min_debias_sample_size` depends on the sizes of the
        leading- and trailing periods, and also on the time resolution gap (=time resolution of the corresponding SensorData).

        References
        -----------
        Jacobs .A, et. al. (2024) `Filling gaps in urban temperature observations by debiasing ERA5 reanalysis data <https://doi.org/10.1016/j.uclim.2024.102226>`_

        """
        self._fillkwargs = {
            "applied_gapfill_method": "diurnal_debias_model_gapfill",
            "leading_period_duration": leading_period_duration,
            "trailing_period_duration": trailing_period_duration,
            "min_debias_sample_size": min_debias_sample_size,
        }

        # 1. Check validity of modeltimeseries
        is_compat, err_msg = gf_methods.check_if_modeltimeseries_is_compatible(
            gap=self,
            modeltimeseries=modeltimeseries,
            lp_duration=leading_period_duration,
            tp_duration=trailing_period_duration,
        )
        if not is_compat:
            self._labels[:] = label_def["failed_diurnal_debias_modeldata_fill"]["label"]
            self._extra_info[:] = err_msg
            logger.warning(
                f"Incompatible modeldata for diurnal_debias_model_gapfill: \n{err_msg}"
            )
            return

        # 2. Construct and validity-test leading and trailing periods
        (
            lead_period,
            trail_period,
            continueflag,
        ) = self._setup_lead_and_trail_for_debias_gapfill(
            sensordata=sensordata,
            fail_label=label_def["failed_diurnal_debias_modeldata_fill"]["label"],
            leading_period_duration=leading_period_duration,
            min_leading_records_total=min_debias_sample_size,
            trailing_period_duration=trailing_period_duration,
            min_trailing_records_total=min_debias_sample_size,
        )
        if not continueflag:
            # warnings and gap attributes are already been updated
            return

        # 3. Fill the gap
        combdf = gf_methods.create_a_combined_df(
            leadseries=lead_period, trailseries=trail_period, gap=self
        )
        # add modeldata to combdf
        combdf = gf_methods.add_modeldata_to_combdf(
            combineddf=combdf, modeltimeseries=modeltimeseries
        )
        # Fill the missing records
        filleddf = fill_with_diurnal_debias(
            df=combdf, min_sample_size=int(min_debias_sample_size)
        )
        filleddf = filleddf.loc[self.records.index]  # subset to gap records

        # 4. Update attributes
        self._records = filleddf["fillvalue"].rename(
            "value"
        )  # set the new filled records

        # set labels
        self._labels.loc[self.records.notna()] = label_def[
            "diurnal_debias_modeldata_fill"
        ]["label"]
        self._labels.loc[self.records.isna()] = label_def[
            "failed_diurnal_debias_modeldata_fill"
        ]["label"]

        # update details
        self._extra_info = filleddf["msg"].rename("details")

    @log_entry
    def weighted_diurnal_debiased_model_gapfill(
        self,
        sensordata: "SensorData",  # type: ignore #noqa: F821
        modeltimeseries: ModelTimeSeries,
        leading_period_duration: pd.Timedelta,
        min_lead_debias_sample_size: int,
        trailing_period_duration: pd.Timedelta,
        min_trail_debias_sample_size: int,
    ) -> None:
        """
        Fill the gaps using a weighted sum of model data corrected for the diurnal bias and weights with respect to the start of the gap.

        This method fills the gap using model data corrected for its diurnal bias.
        The diurnal bias is a bias that is estimated for each timestamp in the leading
        and trailing period (separately). For both periods separately, all biases are averaged over hour, minute and second, to
        obtain a diurnal bias (for each timestamp).

        In addition, a normalized weight is computed for each gap record indicating the distance (in time) to
        the start and end of the gap. The correction applied on the interpolated (in time) model data is
        thus a weighted sum of corrections coming from both the leading and trailing period.

        Parameters
        ----------
        sensordata : SensorData
            The corresponding SensorData used in the computation of the bias. Only
            the observations that are not labeled as outliers are used to compute the bias.
        modeltimeseries : ModelTimeSeries
            The model time series used to fill the gap records. The model data
            must be compatible (equivalent obstype and related to the same Station as the gap.)
        leading_period_duration : pd.Timedelta
            The duration of the leading period. That is the period before the gap, used
            for bias estimation.
        min_lead_debias_sample_size : int
            The minimum number of leading samples required for bias estimation. If this condition is not met, the gap
            is not filled.
        trailing_period_duration : pd.Timedelta
            The duration of the trailing period. That is the period after the gap, used
            for bias estimation.
        min_trail_debias_sample_size : int
            The minimum number of trailing samples required for bias estimation. If this condition is not met, the gap
            is not filled.

        Returns
        --------
        None.


        Notes
        -----
        A schematic description of the weighted diurnal debiased modeldata gap fill:

        #. Check the compatibility of the `ModelTimeSeries` with the `gap`.
        #. Construct a leading and trailing sample, and test if they meet the required conditions.
           The required conditions are tested by testing the samplesizes per hour, minute and second for the leading and trailing periods (seperatly).
        #. A leading and trailing set of diurnal biases are computed by grouping to hour, minute and second, and averaging the biases.
        #. A weight is computed for each gap record, that is the normalized distance to the start and end of the gap.
        #. Fill the gap records by using raw (interpolated) modeldata is corrected by a weighted sum the coresponding diurnal bias for the lead and trail periods.
        #. Update the `gap` attributes with the interpolated values, labels, and details.

        Notes
        --------
        Note that a suitable `min_debias_sample_size` depends on the sizes of the
        leading- and trailing periods, and also on the time resolution gap (=time resolution of the corresponding SensorData).

        References
        -----------
        Jacobs .A, et. al. (2024) `Filling gaps in urban temperature observations by debiasing ERA5 reanalysis data <https://doi.org/10.1016/j.uclim.2024.102226>`_

        """

        self._fillkwargs = {
            "applied_gapfill_method": "weighted_diurnal_debias_model_gapfill",
            "leading_period_duration": leading_period_duration,
            "trailing_period_duration": trailing_period_duration,
            "min_lead_debias_sample_size": min_lead_debias_sample_size,
            "min_trail_debias_sample_size": min_trail_debias_sample_size,
        }

        # 1. Check validity of modeltimeseries
        is_compat, err_msg = gf_methods.check_if_modeltimeseries_is_compatible(
            gap=self,
            modeltimeseries=modeltimeseries,
            lp_duration=leading_period_duration,
            tp_duration=trailing_period_duration,
        )
        if not is_compat:
            self._labels[:] = label_def[
                "failed_weighted_diurnal_debias_modeldata_fill"
            ]["label"]
            self._extra_info[:] = err_msg
            logger.warning(
                f"Incompatible modeldata for weighted_diurnal_debias_model_gapfill: \n{err_msg}"
            )
            return

        # 2. Construct and validity-test leading and trailing periods
        (
            lead_period,
            trail_period,
            continueflag,
        ) = self._setup_lead_and_trail_for_debias_gapfill(
            sensordata=sensordata,
            fail_label=label_def["failed_weighted_diurnal_debias_modeldata_fill"][
                "label"
            ],
            leading_period_duration=leading_period_duration,
            min_leading_records_total=min_lead_debias_sample_size,
            trailing_period_duration=trailing_period_duration,
            min_trailing_records_total=min_trail_debias_sample_size,
        )
        if not continueflag:
            # warnings and gap attributes are already been updated
            return

        # 3. Fill the gap
        combdf = gf_methods.create_a_combined_df(
            leadseries=lead_period, trailseries=trail_period, gap=self
        )
        # add modeldata to combdf
        combdf = gf_methods.add_modeldata_to_combdf(
            combineddf=combdf, modeltimeseries=modeltimeseries
        )
        # Fill the missing records
        filleddf = fill_with_weighted_diurnal_debias(
            df=combdf,
            min_lead_sample_size=min_lead_debias_sample_size,
            min_trail_sample_size=min_trail_debias_sample_size,
        )

        filleddf = filleddf.loc[self.records.index]  # subset to gap records

        # 4. Update attributes
        self._records = filleddf["fillvalue"].rename(
            "value"
        )  # set the new filled records

        # set labels
        self._labels.loc[self.records.notna()] = label_def[
            "weighted_diurnal_debias_modeldata_fill"
        ]["label"]
        self._labels.loc[self.records.isna()] = label_def[
            "failed_weighted_diurnal_debias_modeldata_fill"
        ]["label"]

        # update details
        self._extra_info = filleddf["msg"].rename("details")

    @log_entry
    def raw_model_gapfill(self, modeltimeseries: ModelTimeSeries) -> None:
        """
        Fill the gap using model data without correction.

        This method fills the gap by directly interpolating
        the model data to the missing records.

        Parameters
        ----------
        modeltimeseries : ModelTimeSeries
            The model time series used to fill the gap records. The model data
            must be compatible (equivalent obstype and related to the same Station as the gap.)

        Returns
        -------
        None

        Notes
        -----
        A schematic description of the raw model data gap fill:

        #. Check the compatibility of the `ModelTimeSeries` with the `gap`.
        #. Ensure both the `ModelTimeSeries` and `gap` have the same timezone.
        #. Interpolate the model data to match the missing records in the gap.
        #. Update the `gap` attributes with the interpolated values, labels, and details.

        """
        self._fillkwargs = {"applied_gapfill_method": "raw_model_gapfill"}

        # 1. Check validity of modeltimeseries
        is_compat, err_msg = gf_methods.check_if_modeltimeseries_is_compatible(
            gap=self,
            modeltimeseries=modeltimeseries,
            lp_duration=pd.Timedelta(0),
            tp_duration=pd.Timedelta(0),
        )
        if not is_compat:
            self._labels.loc[self.records.isna()] = label_def[
                "failed_raw_modeldata_fill"
            ]["label"]
            self._extra_info.loc[self.records.isna()] = err_msg
            logger.warning(f"Incompatible modeldata for raw_model_gapfill: \n{err_msg}")
            return

        modelseries = modeltimeseries.series
        gapseries = self.records
        # 2. Ensure both series have the same timezone
        if modelseries.index.tz != gapseries.index.tz:
            modelseries = modelseries.tz_convert(gapseries.index.tz)

        # 3. Fill the gap
        # 3. Reindex modelseries to match gapseries, interpolating if necessary
        modelseries_reindexed = (
            pd.concat([modelseries, gapseries])
            .sort_index()
            .interpolate(method="time", limit_area="inside")
        )
        # duplicates are introduced when timestamps are both in modelseries and gapseries
        modelseries_reindexed = modelseries_reindexed[
            ~modelseries_reindexed.index.duplicated(keep="first")
        ]

        # 4. Update attributes
        self._records = modelseries_reindexed.loc[
            self.records.index
        ]  # (save) set the new filled records

        # set labels
        self._labels.loc[self.records.notna()] = label_def["raw_modeldata_fill"][
            "label"
        ]
        self._labels.loc[self.records.isna()] = label_def["failed_raw_modeldata_fill"][
            "label"
        ]

        # update details
        self._extra_info.loc[self.records.notna()] = (
            f"Successful raw modeldata fill using {modeltimeseries.modelvariable} (but converted to {self.obstype.std_unit}) of {modeltimeseries.modelname}"
        )
        self._extra_info.loc[self.records.isna()] = "Unsuccessful raw modeldata fill."

    @log_entry
    def interpolate(
        self,
        sensordata: "SensorData",  # type: ignore #noqa: F821
        method: str = "time",
        max_consec_fill: int = 10,
        n_leading_anchors: int = 1,
        n_trailing_anchors: int = 1,
        max_lead_to_gap_distance: Union[pd.Timedelta, None] = None,
        max_trail_to_gap_distance: Union[pd.Timedelta, None] = None,
        method_kwargs: dict = {},
    ) -> None:
        """
        Fill the gap using interpolation of SensorData.

        The gap is interpolated using the leading and trailing periods of the gap. One can select different
        interpolation methods. By using restrictions on the leading and trailing periods, one can
        ensure that the interpolation is only done when there are enough leading and trailing data available.

        Parameters
        ----------
        sensordata : SensorData
            The corresponding SensorData used to interpolate the gap.
        method : str, optional
            Interpolation technique to use. See pandas.DataFrame.interpolate
            'method' argument for possible values. Make sure that
            `n_leading_anchors`, `n_trailing_anchors` and `method_kwargs` are
            set accordingly to the method (higher order interpolation techniques require more leading and trailing anchors). The default is "time".
        max_consec_fill : int, optional
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
        max_lead_to_gap_distance : pd.Timedelta or None, optional
            The maximum time difference between the start of the gap and a
            leading anchor(s). If None, no time restriction is applied on the leading anchors. The default is None.
        max_trail_to_gap_distance : pd.Timedelta or None, optional
            The maximum time difference between the end of the gap and a
            trailing anchor(s). If None, no time restriction is applied on the trailing anchors. Defaults to None.
        method_kwargs : dict, optional
            Extra arguments that are passed to pandas.DataFrame.interpolate() structured in a dict. Defaults to {}.

        Notes
        -----
        A schematic description:

        #. Get the leading and trailing periods of the gap.
        #. Check if the leading and trailing periods are valid.
        #. Create a combined DataFrame with the leading, trailing, and gap data.
        #. Interpolate the missing records using the specified method.
        #. Update the gap attributes with the interpolated values, labels, and details.

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
        # store fill settings
        self._fillkwargs = {
            "applied_gapfill_method": "interpolation",
            "method": method,
            "max_consec_fill": max_consec_fill,
            "n_leading_anchors": n_leading_anchors,
            "n_trailing_anchors": n_trailing_anchors,
            "max_lead_to_gap_distance": max_lead_to_gap_distance,
            "max_trail_to_gap_distance": max_trail_to_gap_distance,
            **method_kwargs,
        }

        # 1. Get leading period
        lead_period, continueflag, err_msg = gf_methods.get_leading_period(
            gap=self,
            sensordata=sensordata,
            n_records=n_leading_anchors,
            duration=max_lead_to_gap_distance,
            fixed_by_records=True,
            fixed_by_duration=False,
        )

        if not continueflag:
            # Interpolation failed due to failing leading period
            self._labels[:] = label_def["failed_interpolation_gap"]["label"]
            self._extra_info[:] = err_msg
            logger.warning(
                f"Cannot interpolate {self} because no valid leading period can be found."
            )
            return

        # 2. Get trailing period
        trail_period, continueflag, err_msg = gf_methods.get_trailing_period(
            gap=self,
            sensordata=sensordata,
            n_records=n_trailing_anchors,
            duration=max_trail_to_gap_distance,
            fixed_by_records=True,
            fixed_by_duration=False,
        )

        if not continueflag:
            # Interpolation failed due to failing trailing period
            self._labels[:] = label_def["failed_interpolation_gap"]["label"]
            self._extra_info[:] = err_msg
            logger.warning(
                f"Cannot interpolate {self} because no valid trailing period can be found."
            )
            return

        # 3. Check if the gap records do not exceed the max_consec_fill
        if self.records.shape[0] > max_consec_fill:
            self._labels[:] = label_def["failed_interpolation_gap"]["label"]
            self._extra_info[:] = (
                f"Gap is too large ({self.records.shape[0]} records) to be filled with interpolation (and max_consec_fill={max_consec_fill})."
            )
            logger.warning(
                f"Cannot interpolate {self} because the gap is too large ({self.records.shape[0]} records) to be filled with interpolation (and max_consec_fill={max_consec_fill}). Increase the max_consec_fill or use another gapfill method."
            )
            return

        # 4. Combine the anchors with the observations
        combdf = gf_methods.create_a_combined_df(
            leadseries=lead_period, trailseries=trail_period, gap=self
        )
        tofill_series = combdf["value"]

        # 4. Replace the NaN's (GAPFILLING)
        # Interpolate series
        tofill_series = tofill_series.interpolate(
            method=method,
            limit=max_consec_fill,
            limit_area="inside",
            **method_kwargs,
        )

        # Update attributes
        self._records = tofill_series.loc[
            self.records.index
        ]  # set the new filled records

        # set labels
        self._labels.loc[self.records.notna()] = label_def["interpolated_gap"]["label"]
        self._labels.loc[self.records.isna()] = label_def["failed_interpolation_gap"][
            "label"
        ]

        # update details
        self._extra_info.loc[self.records.notna()] = "Successful interpolation"
        self._extra_info.loc[self.records.isna()] = (
            "Unsuccessful interpolation, likely due to an error when calling pandas.Series.interpolate. See the error logs for further details."
        )

        return

    # ------------------------------------------
    #    Helping methods
    # ------------------------------------------

    @log_entry
    def flush_fill(self) -> None:
        """
        Clear all fill information for this gap.

        This method resets the records, labels, extra info, and fillkwargs to their initial state.
        """
        logger.debug(f"Flushing fill values of {self}")
        # 1. set nan for all records
        self._records.loc[:] = np.nan

        # 2. Convert all labels to 'gap'
        self._labels.loc[:] = label_def["regular_gap"]["label"]

        # 3. Clears the extra data (per record)
        self._extra_info.loc[:] = "no details"

        # 4. Empty the fillkwargs
        self._fillkwargs = {}

    def _setup_lead_and_trail_for_debias_gapfill(
        self,
        sensordata: "SensorData",  # type: ignore #noqa: F821
        fail_label: str,
        leading_period_duration: pd.Timedelta,
        min_leading_records_total: int,
        trailing_period_duration: pd.Timedelta,
        min_trailing_records_total: int,
    ) -> tuple[Union[pd.Series, None], Union[pd.Series, None], bool]:
        """
        Construct leading and trailing periods for debias gapfill.

        This method is shared by multiple gap-filling methods to construct and validate
        the leading and trailing periods required for bias estimation.

        Parameters
        ----------
        sensordata : SensorData
            The corresponding SensorData used to compute the bias.
        fail_label : str
            The label to assign to the gap records if the setup fails.
        leading_period_duration : pd.Timedelta
            The duration of the leading period.
        min_leading_records_total : int
            The minimum number of records required in the leading period.
        trailing_period_duration : pd.Timedelta
            The duration of the trailing period.
        min_trailing_records_total : int
            The minimum number of records required in the trailing period.

        Returns
        -------
        tuple of (pd.Series or None, pd.Series or None, bool)
            A tuple containing the leading period, trailing period, and a flag indicating
            whether the setup was successful. If unsuccessful, the leading and trailing
            periods will be None, and the flag will be False.
        """

        # Validate argument types

        # 1. Get leading period
        lead_period, continueflag, err_msg = gf_methods.get_leading_period(
            gap=self,
            sensordata=sensordata,
            n_records=min_leading_records_total,
            duration=leading_period_duration,
            fixed_by_records=False,
            fixed_by_duration=True,
        )

        if not continueflag:
            # Setup failed due to an invalid leading period
            self._labels[:] = fail_label
            self._extra_info[:] = err_msg
            logger.warning(
                f"Cannot fill {self} because no valid leading period can be found."
            )
            return None, None, False

        # 2. Get trailing period
        trail_period, continueflag, err_msg = gf_methods.get_trailing_period(
            gap=self,
            sensordata=sensordata,
            n_records=min_trailing_records_total,
            duration=trailing_period_duration,
            fixed_by_records=False,
            fixed_by_duration=True,
        )

        if not continueflag:
            # Setup failed due to an invalid trailing period
            self._labels[:] = fail_label
            self._extra_info[:] = err_msg
            logger.warning(
                f"Cannot fill {self} because no valid trailing period can be found."
            )
            return None, None, False

        logger.debug(f"Exiting _setup_lead_and_trail_for_debias_gapfill for {self}")
        return lead_period, trail_period, True
