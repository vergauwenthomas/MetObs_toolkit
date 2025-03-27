import logging
import pandas as pd
import numpy as np

from metobs_toolkit.backend_collection.errorclasses import *
from metobs_toolkit.settings_collection import label_def

import metobs_toolkit.gf_collection.gf_common_methods as gf_methods
from metobs_toolkit.gf_collection.debias_gapfill import fill_regular_debias
from metobs_toolkit.gf_collection.diurnal_debias_gapfill import (
    fill_with_diurnal_debias,
    fill_with_weighted_diurnal_debias,
)

logger = logging.getLogger(__file__)


class Gap:
    def __init__(self, gaprecords: pd.DatetimeIndex, obstype, stationname: str):

        gaprecords.name = "datetime"  # set dtindex name
        # holds the values of the gap (nan when not filled)
        self._records = pd.Series(
            data=np.nan, index=gaprecords, name="value"
        )  # perfect freq records with nan values (initial)

        # holds the labels for each gap record
        self._labels = pd.Series(
            data=label_def["regular_gap"]["label"], index=gaprecords, name="label"
        )

        # holds extra details per record of the gapfill method
        self._extra_info = pd.Series(
            data="no details", index=gaprecords, name="details"
        )

        self._fillkwargs = {}  # store the settings applied to fill method

        self._obstype = obstype
        self._stationname = stationname

    @property
    def records(self):
        return self._records

    @property
    def obstype(self):
        return self._obstype

    @property
    def stationname(self):
        return self._stationname

    @property
    def fillstatus(self):
        if self.records.isna().all() and not bool(self._fillkwargs):
            return "unfilled"

        elif self.records.isna().all() and bool(self._fillkwargs):
            return "failed gapfill"

        elif self.records.isna().any() and bool(self._fillkwargs):
            return "partially succesfull gapfill"

        elif not self.records.isna().any() and bool(self._fillkwargs):
            return "succesfull gapfill"
        else:
            raise NotImplementedError(
                "This situation is unforseen! Please notify developers."
            )

    @property
    def start_datetime(self):
        return min(self.records.index)

    @property
    def end_datetime(self):
        return max(self.records.index)

    @property
    def df(self):
        return pd.DataFrame(
            {"value": self.records, "label": self._labels, "details": self._extra_info}
        )

    # ------------------------------------------
    #    Get info methods
    # ------------------------------------------

    def flag_can_be_filled(self, overwrite: bool = False):
        if overwrite:
            return True
        if self.fillstatus == "unfilled":
            return True
        else:
            return False

    def get_info(self, printout: bool = True):
        """Print out detailed info about the Gap.

        Returns
        -------
        None.

        """

        retstr = ""
        retstr += "---- Gap info -----\n"

        retstr += f"  * Gap of {self.obstype.name} for station: {self.stationname} \n"
        retstr += f"  * Start gap: {self.start_datetime}\n"
        retstr += f"  * End gap: {self.end_datetime}\n"
        retstr += f"  * Duration gap: {self.end_datetime - self.start_datetime}\n"
        retstr += f"  * Gap status: {self.fillstatus}"

        if printout:
            print(retstr)
        else:
            return retstr

    def debiased_model_gapfill(
        self,
        sensordata,
        modeltimeseries,
        leading_period_duration: pd.Timedelta,
        min_leading_records_total: int,
        trailing_period_duration: pd.Timedelta,
        min_trailing_records_total: int,
    ):

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
        lead_period, trail_period, continueflag = (
            self._setup_lead_and_trail_for_debias_gapfill(
                sensordata=sensordata,
                fail_label=label_def["failed_debias_modeldata_fill"]["label"],
                leading_period_duration=leading_period_duration,
                min_leading_records_total=min_leading_records_total,
                trailing_period_duration=trailing_period_duration,
                min_trailing_records_total=min_trailing_records_total,
            )
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
        filleddf = fill_regular_debias(df=combdf)
        filleddf = filleddf.loc[self.records.index]  # subset to gap records

        # 4. Update attributes
        # Update attributes
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

    def diurnal_debiased_model_gapfill(
        self,
        sensordata,
        modeltimeseries,
        leading_period_duration: pd.Timedelta,
        trailing_period_duration: pd.Timedelta,
        min_debias_sample_size: int,
    ):

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
        lead_period, trail_period, continueflag = (
            self._setup_lead_and_trail_for_debias_gapfill(
                sensordata=sensordata,
                fail_label=label_def["failed_diurnal_debias_modeldata_fill"]["label"],
                leading_period_duration=leading_period_duration,
                min_leading_records_total=min_debias_sample_size,  # at least!
                trailing_period_duration=trailing_period_duration,
                min_trailing_records_total=min_debias_sample_size,  # at least!
            )
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
        # Update attributes
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

    def weighted_diurnal_debiased_model_gapfill(
        self,
        sensordata,
        modeltimeseries,
        leading_period_duration: pd.Timedelta,
        min_lead_debias_sample_size: int,
        trailing_period_duration: pd.Timedelta,
        min_trail_debias_sample_size: int,
    ):

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
        lead_period, trail_period, continueflag = (
            self._setup_lead_and_trail_for_debias_gapfill(
                sensordata=sensordata,
                fail_label=label_def["failed_weighted_diurnal_debias_modeldata_fill"][
                    "label"
                ],
                leading_period_duration=leading_period_duration,
                min_leading_records_total=min_lead_debias_sample_size,  # at least!
                trailing_period_duration=trailing_period_duration,
                min_trailing_records_total=min_trail_debias_sample_size,  # at least!
            )
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
        # Update attributes
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

    def raw_model_gapfill(self, modeltimeseries):
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
        # dulicates are introduced when timestamps are both in modelseries and gapseries
        modelseries_reindexed = modelseries_reindexed[
            ~modelseries_reindexed.index.duplicated(keep="first")
        ]

        # 4. Update attributes
        # Update attributes
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
            f"Succesfull raw modeldata fill using {modeltimeseries.modelvariable} (but converted to {self.obstype.std_unit}) of {modeltimeseries.modelname}"
        )
        self._extra_info.loc[self.records.isna()] = f"Unsuccesfull raw modeldata fill."

    def interpolate(
        self,
        sensordata,
        method="time",
        max_consec_fill=10,
        n_leading_anchors=1,
        n_trailing_anchors=1,
        max_lead_to_gap_distance: pd.Timedelta | None = None,
        max_trail_to_gap_distance: pd.Timedelta | None = None,
        method_kwargs={},
    ):
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
                f"Cannot interpolate {self} because leading no valid leading period can be found."
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
                f"Cannot interpolate {self} because leading no valid trailing period can be found."
            )
            return

        # 3. Combine the anchors with the observations
        combdf = gf_methods.create_a_combined_df(
            leadseries=lead_period, trailseries=trail_period, gap=self
        )
        tofill_series = combdf["value"]

        # 4. Replace the NaN's (GAPFILLING)
        # Interpolate series
        tofill_series = tofill_series.interpolate(
            method=method,
            limit=max_consec_fill,  # Maximum number of consecutive NaNs to fill. Must be greater than 0.
            limit_area="inside",
            **method_kwargs,
        )

        # Update attributes
        self._records = tofill_series.loc[
            self.records.index
        ]  # set the new filled records

        # set labels
        self._labels.loc[self.records.notna()] = label_def["interpolated_gap"][
            "label"
        ]  # interpolated
        self._labels.loc[self.records.isna()] = label_def["failed_interpolation_gap"][
            "label"
        ]  # failed

        # update details
        self._extra_info.loc[self.records.notna()] = "Succesfull interpolation"
        self._extra_info.loc[self.records.isna()] = (
            f"Unsuccesfull interpolation, likely due to exceeding maximumn number of consec fill (={max_consec_fill})."
        )

        return

    # ------------------------------------------
    #    Helping methods
    # ------------------------------------------

    def flush_fill(self):
        """Clears all fill info"""

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
        sensordata,
        fail_label,
        leading_period_duration: pd.Timedelta,
        min_leading_records_total: int,
        trailing_period_duration: pd.Timedelta,
        min_trailing_records_total: int,
    ):
        """In a seperate method because it is shared by multiple methods"""

        # 2. Get leading period
        lead_period, continueflag, err_msg = gf_methods.get_leading_period(
            gap=self,
            sensordata=sensordata,
            n_records=min_leading_records_total,
            duration=leading_period_duration,
            fixed_by_records=False,
            fixed_by_duration=True,
        )

        if not continueflag:
            # Interpolation failed due to failing leading period
            self._labels[:] = fail_label
            self._extra_info[:] = err_msg
            logger.warning(
                f"Cannot fill {self} because leading no valid leading period can be found."
            )
            return None, None, False

        # 3. Get trailing period
        trail_period, continueflag, err_msg = gf_methods.get_trailing_period(
            gap=self,
            sensordata=sensordata,
            n_records=min_trailing_records_total,
            duration=trailing_period_duration,
            fixed_by_records=False,
            fixed_by_duration=True,
        )

        if not continueflag:
            # Interpolation failed due to failing trailing period
            self._labels[:] = fail_label
            self._extra_info[:] = err_msg
            logger.warning(
                f"Cannot fill {self} because leading no valid trailing period can be found."
            )
            return None, None, False

        return lead_period, trail_period, True
