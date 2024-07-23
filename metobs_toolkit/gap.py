#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module contains the Gap class and all its methods.

A Gap contains all information and methods of a data-gap.
"""
import numpy as np
import pandas as pd

import logging
import datetime
import math


import metobs_toolkit.gap_filling as gap_filling

from metobs_toolkit.df_helpers import (
    format_outliersdf_to_doubleidx,
    init_multiindexdf,
    concat_save,
    get_likely_frequency,
    _find_closes_occuring_date,
)

from metobs_toolkit.obstypes import Obstype as Obstype_class
from metobs_toolkit.df_helpers import init_multiindex, xs_save

from metobs_toolkit.settings_files.default_formats_settings import (
    label_def,
)  # to assign labels

logger = logging.getLogger(__name__)


# =============================================================================
# Gap class

# a gap is a sequence of repeting missing obs
# =============================================================================


class Gap:
    """Gap class holds all gap information and methods for gaps."""

    def __init__(self, name, startdt, enddt, obstype, records_freq):
        """
        Initiate Gap object.


        Parameters
        ----------
        name : String
            Station name where the gap occures.
        startdt : datetime.datetime
            Start datetime of the gap (included).
        enddt : datetime.datetime
            End datetime of the gap (included)
        obstype : metobs_toolkit.Obstype
           The corresponding observation type of the gap.
        records_freq : pd.Timedelta()
           The assumded frequency of records for this observationtype.

        Note
        -----
        The records_freq is required for the gap to expand into missing records.

        Returns
        -------
        None

        """
        # Set key attributes
        self.name = str(name)

        assert isinstance(
            startdt, datetime.datetime
        ), f"{startdt} not a datetime instance."
        assert isinstance(enddt, datetime.datetime), f"{enddt} not a datetime instance."
        assert enddt >= startdt, f"Enddt: {enddt} is not after startdt: {startdt}."

        self.startdt = startdt
        self.enddt = enddt
        self.duration = enddt - startdt

        assert isinstance(
            obstype, Obstype_class
        ), f"{obstype} is not a metobs_toolkit.Obstype."
        self.obstype = obstype

        assert isinstance(
            records_freq, pd.Timedelta
        ), f"{records_freq} is not a pd.Timedelta."
        self.records_freq = records_freq

        # gapdf holds all the records in the gap (depending on the freq of the dataset!)
        self.gapdf = (
            init_multiindexdf()
        )  # columns: obstype.name, obstype.name"_fill', "fill_method", "msg"
        # anchordf holds all the records outside the gap, used by fill methods.
        self.anchordf = (
            init_multiindexdf()
        )  # columns: obstype.name, obstype.name"_fill', "fill_method", "msg"

        # Construct derived attributes and initiate them
        self._setup_of_gapdf()

    def __str__(self):
        return f"{self.obstype.name}-gap of {self.name} for {self.startdt} --> {self.enddt}, duration: {self.duration} or {self.gapdf.shape[0]} records."

    def __repr__(self):
        return self.__str__()

    # =============================================================================
    # Setup (methods automatically applied after creation)
    # =============================================================================
    def _setup_of_gapdf(self):
        """Setup of the gapdf attribute from the core attributes."""
        # Compute missing records
        missing_records = pd.date_range(
            start=self.startdt,
            end=self.enddt,
            freq=self.records_freq,
            tz=self.startdt.tz,
        )

        # ----  Construct gapdf ------
        gapdf = init_multiindexdf()
        # Create index
        gapdf.index = pd.MultiIndex.from_arrays(
            arrays=[
                [self.name] * len(missing_records),
                missing_records,
            ],
            names=["name", "datetime"],
        )

        # Set columns with default values
        gapdf[f"{self.obstype.name}"] = np.nan  # fill with nan values
        gapdf[f"{self.obstype.name}_fill"] = np.nan  # fill with nan values
        gapdf["fill_method"] = "not filled"  # on gap creation, will be overwritten
        gapdf["msg"] = None  # Nonedataset

        self._set_gapdf(gapdf)

    def _set_gapdf(self, gapdf):
        # TODO: checks
        if list(gapdf.index.names) != ["name", "datetime"]:
            raise MetobsGapError(f"The gapdf does not have a correct index: {gapdf}")

        try:
            if gapdf.index.get_level_values("datetime").tz is None:
                raise MetobsGapError(
                    f"The gapdf-datetimeindex is tz-naive. {gapdf.index}"
                )
            pass
        except:
            raise MetobsGapError(
                f"The gapdf-datetimeindex is not a datetime dtype. {gapdf.index}"
            )
        if list(gapdf.columns) != [
            f"{self.obstype.name}",
            f"{self.obstype.name}_fill",
            "fill_method",
            "msg",
        ]:
            raise MetobsGapError(f"The gapdf does not have a correct columns: {gapdf}")

        self.gapdf = gapdf

    def _set_anchordf(self, anchordf):
        # TODO : Checks
        if list(anchordf.index.names) != ["name", "datetime"]:
            raise MetobsGapError(
                f"The anchordf does not have a correct index: {anchordf}"
            )
        try:
            if anchordf.index.get_level_values("datetime").tz is None:
                raise MetobsGapError(
                    f"The gapdf-datetimeindex is tz-naive. {anchordf.index}"
                )
            pass
        except:
            raise MetobsGapError(
                f"The gapdf-datetimeindex is not a datetime dtype. {anchordf.index}"
            )

        if list(anchordf.columns) != [self.obstype.name, "fill_method", "msg"]:
            raise MetobsGapError(
                f"The anchordf does not have a correct columns: {anchordf}"
            )

        self.anchordf = anchordf

    def _can_be_filled(self, overwrite_fill):
        """Returns True if no fill values are present."""

        if self.gapdf[f"{self.obstype.name}_fill"].isnull().all():
            return True
        elif overwrite_fill:
            logger.warning(f"The filling of {self}, will be overwritten!")
            self._setup_of_gapdf()  # reset gapfill
            return True

        else:
            return False

    def _get_gapfill_status(self):
        # Warning! do not change these labels (or replace all hardcode occurences)
        # if not self.does_gap_holds_missing_records():
        #     return "Gap does not exist in observation space"
        if self.gapdf[f"{self.obstype.name}_fill"].isnull().all():
            return "Unfilled gap"
        if self.gapdf[f"{self.obstype.name}_fill"].isnull().any():
            return "Partially filled gap"
        else:
            return "Filled gap"

    def get_info(self):
        """Print detailed information of a gap."""

        print("---- Gap info -----")
        print(
            "(Note: gaps start and end are defined on the frequency estimation of the native dataset.)"
        )
        print(f"  * Gap for station: {self.name}")
        print(f"  * Start gap: {self.startdt}")
        print(f"  * End gap: {self.enddt}")
        print(f"  * Duration gap: {self.duration}")
        print(f"  * For {self.obstype.name}")
        print(f"  * Gapfill status >>>> {self._get_gapfill_status()}")
        # if not self.does_gap_holds_missing_records():
        #     print(
        #         "The gap is unreal (=no missing records at the current frequency resolution)."
        #     )

        if not self.anchordf.empty:
            print("---- Anchor Data Frame ----- ")
            print(self.anchordf)

        print("---- Gap Data Frame -----")
        print(self.gapdf)

    # def get_info(self):
    #     """Print detailed information of a gap."""
    #     print(f"Gap for {self.name} with:")
    #     print("---- Gap info -----")
    #     print(
    #         "(Note: gaps are defined on the frequency estimation of the native dataset.)"
    #     )
    #     print(f"  * Start gap: {self.startgap}")
    #     print(f"  * End gap: {self.endgap}")
    #     print(f"  * Duration gap: {self.duration}")
    #     print("---- Gap fill info -----")
    #     obstypes = self.gapfill_df.columns.to_list()
    #     obstypes = [obst for obst in obstypes if not obst.endswith("_final_label")]
    #     if self.gapfill_df.empty:
    #         print("(No gapfill applied)")
    #     elif self.gapfill_technique == "gap_interpolation":
    #         for obstype in obstypes:
    #             print(f"  * On observation type: {obstype}")
    #             print(f"  * Technique: {self.gapfill_technique}")
    #             if bool(self.leading_val):
    #                 leading_val = self.leading_val[obstype]
    #             else:
    #                 leading_val = "No leading observation value"
    #             print(
    #                 f"  * Leading timestamp: {self.leading_timestamp} with  {obstype} = {leading_val}"
    #             )
    #             if bool(self.trailing_val):
    #                 trailing_val = self.trailing_val[obstype]
    #             else:
    #                 trailing_val = "No trailing observation value"
    #             print(
    #                 f"  * Trailing timestamp: {self.trailing_timestamp} with  {obstype} = {trailing_val}"
    #             )
    #             print(f"  * Filled values: {self.gapfill_df[obstype]}")
    #             if obstype in self.gapfill_errormessage:
    #                 print(f"  * Gapfill message: {self.gapfill_errormessage[obstype]}")
    #             if self.gapfill_info is not None:
    #                 print(f"  * Gapfill info: {self.gapfill_info.head()}")
    #                 print(
    #                     "    (Extract the gapfill info dataframe by using the .gapfill_info attribute)"
    #                 )

    #     elif self.gapfill_technique == "gap_debiased_era5":
    #         for obstype in obstypes:
    #             print(f"  * On observation type: {obstype}")
    #             print(f"  * Technique: {self.gapfill_technique}")
    #             # print(f'  * Leading timestamp: {self.leading_timestamp} with  {obstype} = {self.leading_val[obstype]}\n')
    #             # print(f'  * Trailing timestamp: {self.trailing_timestamp} with  {obstype} = {self.trailing_val[obstype]}\n')
    #             print(f"  * Filled values: {self.gapfill_df[obstype]}")
    #             if obstype in self.gapfill_errormessage:
    #                 print(f"  * Gapfill message: {self.gapfill_errormessage[obstype]}")
    #             if self.gapfill_info is not None:
    #                 print(f"  * Gapfill info: {self.gapfill_info.head()}")
    #                 print(
    #                     "    (Extract the gapfill info dataframe by using the .gapfill_info attribute)"
    #                 )

    #     else:
    #         print("technique not implemented in yet in show")

    # def to_df(self):
    #     """Convert a Gap object to a dataframe (with one row).

    #     The station name is the index and two colums ('start_gap', 'end_gap')
    #     are constructed.

    #     Returns
    #     -------
    #     pandas.DataFrame()
    #         Gap in dataframe format.

    #     """
    #     returndf = pd.DataFrame(
    #         index=[self.name],
    #         data={
    #             "start_gap": self.startgap,
    #             "end_gap": self.endgap,
    #             "duration": self.duration,
    #         },
    #     )
    #     returndf.index.name = "name"
    #     return returndf

    # def update_leading_trailing_obs(self, obsdf, outliersdf, obs_only=False):
    #     """Update leading and trailing periods in the attributes.

    #     Add the leading (last obs before gap) and trailing (first obs after gap)
    #     as extra columns to the self.df.

    #     One can specify to look for leading and trailing in the obsdf or in both
    #     the obsdf and outliersdf.

    #     The gap leading and trailing timestamps and value attributes are updated.

    #     If no leading/trailing timestamp is found, it is set to the gaps startdt/enddt.

    #     Parameters
    #     ----------
    #     obsdf : pandas.DataFrame
    #         Dataset.df
    #     outliersdf : pandas.DataFrame
    #         Dataset.outliersdf
    #     obs_only: bool, optional
    #         If True, only the obsdf will be used to search for leading and trailing.

    #     Returns
    #     -------
    #     None.

    #     """
    #     sta_obs = xs_save(obsdf, self.name, level="name").index
    #     if obs_only:
    #         sta_comb = sta_obs
    #     else:

    #         outliersdf = format_outliersdf_to_doubleidx(outliersdf)

    #         # combine timestamps of observations and outliers
    #         sta_outl = xs_save(outliersdf, self.name, level="name").index
    #         if sta_outl.empty:
    #             sta_comb = sta_obs
    #         else:
    #             sta_comb = sta_obs.append(sta_outl)

    #     # find minimium timediff before
    #     before_diff = _find_closes_occuring_date(
    #         refdt=self.startgap, series_of_dt=sta_comb, where="before"
    #     )

    #     # if no timestamps are before gap, assume gap at the start of the observations
    #     if math.isnan(before_diff):
    #         before_diff = 0.0

    #     # find minimum timediff after gap
    #     after_diff = _find_closes_occuring_date(
    #         refdt=self.endgap, series_of_dt=sta_comb, where="after"
    #     )
    #     # if no timestamps are after gap, assume gap at the end of the observations
    #     if math.isnan(after_diff):
    #         after_diff = 0.0

    #     # get before and after timestamps
    #     self.leading_timestamp = self.startgap - datetime.timedelta(seconds=before_diff)
    #     self.trailing_timestamp = self.endgap + datetime.timedelta(seconds=after_diff)

    #     # get the values
    #     try:
    #         self.leading_val = obsdf.loc[(self.name, self.leading_timestamp)].to_dict()
    #     except KeyError:
    #         logger.warning("Leading value not found in the observations")
    #         self.leading_val = {}
    #     try:
    #         self.trailing_val = obsdf.loc[
    #             (self.name, self.trailing_timestamp)
    #         ].to_dict()
    #     except KeyError:
    #         logger.warning("Trailing value not found in the observations")
    #         self.trailing_val = {}

    # def update_gaps_indx_in_obs_space(self, obsdf, outliersdf, dataset_res):
    #     """Get the gap records in observation-space.

    #     Explode the gap, to the dataset resolution and format to a multiindex
    #     with name -- datetime.

    #     In addition the last observation before the gap (leading), and first
    #     observation (after) the gap are computed and stored in the df attribute.
    #     (the outliers are used to look for leading and trailing observations.)

    #     Parameters
    #     ----------
    #     obsdf : Dataset.df
    #         The Dataset.df attribute. (Needed to extract trailing/leading
    #                                    observations.)
    #     outliersdf : Dataset.outliersdf
    #         The Dataset.outliersdf attribute.(Needed to extract trailing/leading
    #                                           observations.))
    #     resolutionseries : Datetime.timedelta
    #         Resolution of the station observations in the dataset.

    #     Returns
    #     -------
    #     None

    #     """
    #     outliersdf = format_outliersdf_to_doubleidx(outliersdf)
    #     self.update_leading_trailing_obs(obsdf, outliersdf)

    #     gaprange = pd.date_range(
    #         start=self.leading_timestamp,
    #         end=self.trailing_timestamp,
    #         freq=dataset_res,
    #         inclusive="neither",
    #     )

    #     self.exp_gap_idx = pd.MultiIndex.from_arrays(
    #         arrays=[[self.name] * len(gaprange), gaprange], names=["name", "datetime"]
    #     )

    # =============================================================================
    #         Gapfill
    # =============================================================================
    def interpolate_gap(
        self,
        records,
        method="time",
        max_consec_fill=10,
        max_lead_to_gap_distance=None,
        max_trail_to_gap_distance=None,
    ):
        gapdf = self.gapdf

        obsname = self.obstype.name
        sta_obs_series = xs_save(records, self.name, "name", drop_level=True)
        sta_obs_series = xs_save(sta_obs_series, obsname, "obstype", drop_level=True)
        sta_obs_series = sta_obs_series["value"]

        # 1. Get leading and trailing info
        # get leading record, check validity and add to the gapfilldf
        (_, lead_dt, lead_val, lead_msg) = self.get_leading_record(
            observations_series=sta_obs_series,
            max_lead_to_gap_distance=max_lead_to_gap_distance,
        )
        (_, trail_dt, trail_val, trail_msg) = self.get_trailing_record(
            observations_series=sta_obs_series,
            max_trail_to_gap_distance=max_trail_to_gap_distance,
        )

        # 2. Update the anchordf
        idx = pd.MultiIndex.from_arrays(
            arrays=[[self.name, self.name], [lead_dt, trail_dt]],
            names=["name", "datetime"],
        )

        anchor_df = pd.DataFrame(
            data={
                obsname: [lead_val, trail_val],
                "fill_method": ["leading", "trailing"],
                "msg": [lead_msg, trail_msg],
            },
            index=idx,
        )
        # Update attributes
        self._set_anchordf(anchor_df)

        gapdf["fill_method"] = label_def["interpolated_gap"]["label"]

        # check if gapfill can proceed
        if (lead_msg != "ok") & (trail_msg != "ok"):
            logger.warning(
                f"Cannot fill {self}, because leading and trailing records are not valid."
            )
            print(
                f"Warning! Cannot fill {self}, because leading and trailing periods are not valid."
            )

            gapdf[f"{obsname}_fill"] = np.nan
            gapdf["msg"] = f"{lead_msg} and {trail_msg}"

            self._set_gapdf(gapdf)
            return
        elif (lead_msg != "ok") & (trail_msg == "ok"):
            logger.warning(f"Cannot fill {self}, because leading record is not valid.")
            print(f"Warning! Cannot fill {self}, because leading record is not valid.")

            gapdf[f"{obsname}_fill"] = np.nan
            gapdf["msg"] = f"{lead_msg}"
            self._set_gapdf(gapdf)
            return
        elif (lead_msg == "ok") & (trail_msg != "ok"):
            logger.warning(f"Cannot fill {self}, because trailing record is not valid.")
            print(f"Warning! Cannot fill {self}, because trailing record is not valid.")

            gapdf[f"{obsname}_fill"] = np.nan
            gapdf["msg"] = f"{trail_msg}"
            self._set_gapdf(gapdf)
            return
        else:
            pass

        # 3. Combine the anchors with the observations
        tofilldf = pd.concat([gapdf[[obsname]], self.anchordf[[obsname]]])

        # 4. Replace the NaN's (GAPFILLING)
        # make shure only datetimes are in the index, and sorted by datetimes
        tofilldf = tofilldf.reset_index().set_index("datetime").sort_index()

        # Interpolate series
        tofilldf[obsname] = tofilldf[obsname].interpolate(
            method=method,
            limit=max_consec_fill,  # Maximum number of consecutive NaNs to fill. Must be greater than 0.
            limit_area="inside",
        )
        # revert to multiindex
        tofilldf = tofilldf.reset_index().set_index(["name", "datetime"]).sort_index()

        # 5. Get messages for all filled records
        # filter to gap records
        tofilldf = tofilldf.loc[gapdf.index]

        # messager
        def _msg_interp(row, obstypename):
            if np.isnan(row[obstypename]):
                return "Permitted_by_max_consecutive_fill"
            else:
                return "ok"

        if (lead_msg == "ok") & (
            trail_msg == "ok"
        ):  # labels are result of interpolation
            tofilldf["msg"] = tofilldf.apply(
                _msg_interp, axis="columns", args=([obsname])
            )
        elif (lead_msg == "ok") & (
            trail_msg != "ok"
        ):  # labels are result of invalid trailing
            tofilldf["msg"] = "Invalid trailing anchor"
        elif (lead_msg != "ok") & (
            trail_msg == "ok"
        ):  # labels are result of invalid leading
            tofilldf["msg"] = "Invalid leading anchor"
        elif (lead_msg != "ok") & (
            trail_msg != "ok"
        ):  # labels are result of invalid trailing
            tofilldf["msg"] = "Invalid leading and trailing anchor"
        else:
            raise MetobsGapError(
                "unforseen situation of lead and trail messages in interpolation method."
            )

        # 6. Update the gapdf
        gapdf[f"{obsname}_fill"] = tofilldf[obsname]  # update values
        gapdf["msg"] = tofilldf["msg"]  # give the same msg to all

        # 7. Add FAIL labels if needed
        gapdf.loc[gapdf[f"{obsname}_fill"].isna(), "fill_method"] = label_def[
            "failed_interpolation_gap"
        ]["label"]

        self._set_gapdf(gapdf)

        return

    # =============================================================================
    # Gapfill using modeldata
    # =============================================================================

    def raw_model_gapfill(self, Dataset, Modeldata):
        obsname = self.obstype.name

        # 1. Get leading and trailing info
        # 2. Update the anchordf
        # 3. Combine the anchors with the observations
        label = label_def["raw_modeldata_fill"]["label"]
        fail_label = label_def["failed_raw_modeldata_fill"]["label"]

        # add the gap period
        filldf = Modeldata.interpolate_modeldata(self.gapdf.index)
        filldf = filldf[[obsname]]  # get only relevant obstype
        filldf = filldf.rename(columns={obsname: "modelvalues"})

        # combine the modeldata and gapdata in one df
        filldf = pd.merge(
            left=self.gapdf, right=filldf, how="left", left_index=True, right_index=True
        )

        # 4. Fill the gap
        fill_values, msg = gap_filling._raw_modeldata_fill(filldf)

        # 5. create a gapdf
        gapdf = filldf
        gapdf = gapdf[[f"{obsname}", f"{obsname}_fill", "fill_method", "msg"]]

        gapdf[f"{obsname}_fill"] = fill_values
        gapdf["msg"] = msg

        gapdf["fill_method"] = label_def["raw_modeldata_fill"]["label"]
        gapdf.loc[gapdf[f"{obsname}_fill"].isna(), "fill_method"] = label_def[
            "failed_raw_modeldata_fill"
        ]["label"]

        self._set_gapdf(gapdf)

    def debias_model_gapfill(
        self,
        Dataset,
        Modeldata,
        leading_period_duration="24H",
        min_leading_records_total=60,
        trailing_period_duration="24H",
        min_trailing_records_total=60,
    ):

        # 1. Get leading and trailing info
        obsname = self.obstype.name

        # Create anchor periods
        anchordf = gap_filling._create_anchor_df_for_leading_trailing_periods(
            Gap=self,
            Dataset=Dataset,
            leading_period_duration=leading_period_duration,
            trailing_period_duration=trailing_period_duration,
        )
        anchordf["msg"] = "_init_"
        # Check validity of anchors
        # 1. Check leading period total sample size and add msg
        if (
            anchordf[anchordf["fill_method"] == "leading period"].shape[0]
            < min_leading_records_total
        ):
            msg_lead = f"minimum records ({min_leading_records_total}) for leading period not met."
        else:
            msg_lead = "ok"

        anchordf.loc[anchordf["fill_method"] == "leading period", "msg"] = msg_lead

        # 1B. Check trailing period total sample size and add msg
        if (
            anchordf[anchordf["fill_method"] == "trailing period"].shape[0]
            < min_trailing_records_total
        ):
            msg_trail = f"minimum records ({min_trailing_records_total}) for trailing period not met."
        else:
            msg_trail = "ok"

        anchordf.loc[anchordf["fill_method"] == "trailing period", "msg"] = msg_lead

        # 2. Update the Anchor attribute
        self._set_anchordf(anchordf)

        # 3.Test if anchor period is valid, if not write default attribute
        if (msg_lead != "ok") & (msg_trail != "ok"):
            logmsg = f"Cannot fill {self}, because leading and trailing periods are not valid."
            gapmsg = f"{msg_lead} and {msg_trail}"
            err = True
        elif (msg_lead != "ok") & (msg_trail == "ok"):
            logmsg = f"Cannot fill {self}, because leading period is not valid."
            gapmsg = f"{msg_lead}"
            err = True
        elif (msg_lead == "ok") & (msg_trail != "ok"):
            logmsg = f"Cannot fill {self}, because trailing period is not valid."
            gapmsg = f"{msg_trail}"
            err = True
        else:
            err = False

        gapdf = self.gapdf
        if err:
            logger.warning(logmsg)
            print("Warning! ", logmsg)
            gapdf[f"{obsname}_fill"] = np.nan
            gapdf["msg"] = gapmsg
            gapdf["fill_method"] = label_def["failed_debias_modeldata_fill"]["label"]
        else:

            # 4. combine learning and gap period
            filldf = gap_filling._combine_learning_and_gap_to_one_df(
                Modeldata=Modeldata,
                anchordf=self.anchordf,
                gapdf=gapdf,
                obsname=self.obstype.name,
            )

            # 5. Fill the gap
            fill_values, msg = gap_filling._debias_modeldata_fill(filldf)

            # 6. Update gapdf attribute
            gapdf[f"{obsname}_fill"] = fill_values
            gapdf["msg"] = msg
            gapdf["fill_method"] = label_def["debias_modeldata_fill"]["label"]
            gapdf.loc[gapdf[f"{obsname}_fill"].isna(), "fill_method"] = label_def[
                "failed_debias_modeldata_fill"
            ]["label"]

        # 7. Update gapdf attribute
        self._set_gapdf(gapdf)

    def diurnal_debias_model_gapfill(
        self,
        Dataset,
        Modeldata,
        leading_period_duration="24H",
        min_debias_sample_size=6,
        trailing_period_duration="24H",
    ):
        # 1. Get leading and trailing info

        # Create anchor periods
        anchordf = gap_filling._create_anchor_df_for_leading_trailing_periods(
            Gap=self,
            Dataset=Dataset,
            leading_period_duration=leading_period_duration,
            trailing_period_duration=trailing_period_duration,
        )

        # Check validity of anchors and label them
        anchordf = gap_filling._label_anchors_for_diurnal_gapfilling(
            anchordf=anchordf,
            gapdf=self.gapdf,
            obsname=self.obstype.name,
            min_anchors_per_diurnal_timestamp=min_debias_sample_size,
        )

        self._set_anchordf(anchordf)

        # 4. combine learning and gap period
        filldf = gap_filling._combine_learning_and_gap_to_one_df(
            Modeldata=Modeldata,
            anchordf=self.anchordf,
            gapdf=self.gapdf,
            obsname=self.obstype.name,
        )

        # 5. apply blacklisting (remove anchors that do not fulfill criteria)
        blacklist = self.anchordf[self.anchordf["msg"] != "ok"].index
        filldf = filldf[~filldf.index.isin(blacklist)]

        # 6. Fill the gap
        fill_values, msg = gap_filling._diurnal_debias_modeldata_fill(filldf)

        # 7. Update gapdf attribute
        gapdf = self.gapdf
        gapdf[f"{self.obstype.name}_fill"] = fill_values
        gapdf["msg"] = msg
        gapdf["fill_method"] = label_def["diurnal_debias_modeldata_fill"]["label"]
        gapdf.loc[gapdf[f"{self.obstype.name}_fill"].isna(), "fill_method"] = label_def[
            "failed_diurnal_debias_modeldata_fill"
        ]["label"]

        self._set_gapdf(gapdf)

    def weighted_diurnal_debias_model_gapfill(
        self,
        Dataset,
        Modeldata,
        leading_period_duration="48h",
        min_lead_debias_sample_size=3,
        trailing_period_duration="48h",
        min_trail_debias_sample_size=3,
    ):

        obsname = self.obstype.name
        # Create anchor periods
        # TODO: make shure datetime index is datetime dtype
        anchordf = gap_filling._create_anchor_df_for_leading_trailing_periods(
            Gap=self,
            Dataset=Dataset,
            leading_period_duration=leading_period_duration,
            trailing_period_duration=trailing_period_duration,
        )

        # Check validity of anchors and label them
        lead_anchorsdf = anchordf[anchordf["fill_method"] == "leading period"]
        lead_anchorsdf = gap_filling._label_anchors_for_diurnal_gapfilling(
            anchordf=lead_anchorsdf,
            gapdf=self.gapdf,
            obsname=self.obstype.name,
            min_anchors_per_diurnal_timestamp=min_lead_debias_sample_size,
        )

        trail_anchorsdf = anchordf[anchordf["fill_method"] == "trailing period"]
        trail_anchorsdf = gap_filling._label_anchors_for_diurnal_gapfilling(
            anchordf=trail_anchorsdf,
            gapdf=self.gapdf,
            obsname=self.obstype.name,
            min_anchors_per_diurnal_timestamp=min_trail_debias_sample_size,
        )

        anchordf = pd.concat([lead_anchorsdf, trail_anchorsdf])
        self._set_anchordf(anchordf)

        # 4. combine learning and gap period
        filldf = gap_filling._combine_learning_and_gap_to_one_df(
            Modeldata=Modeldata,
            anchordf=self.anchordf,
            gapdf=self.gapdf,
            obsname=self.obstype.name,
        )

        # 5. apply blacklisting (remove anchors that do not fulfill criteria)
        blacklist = self.anchordf[self.anchordf["msg"] != "ok"].index
        filldf = filldf[~filldf.index.isin(blacklist)]

        # 6. Fill the gap
        fill_values, msg = gap_filling._weighted_diurnal_debias_modeldata(filldf)

        # 7. Update gapdf attribute
        gapdf = self.gapdf
        gapdf[f"{obsname}_fill"] = fill_values
        gapdf["msg"] = msg
        gapdf["fill_method"] = label_def["weighted_diurnal_debias_modeldata_fill"][
            "label"
        ]
        gapdf.loc[gapdf[f"{self.obstype.name}_fill"].isna(), "fill_method"] = label_def[
            "failed_weighted_diurnal_debias_modeldata_fill"
        ]["label"]

        self._set_gapdf(gapdf)

    # =============================================================================
    #  Gapfill anchor methods
    # =============================================================================
    def get_leading_period(
        self,
        observations_series,
        leading_period_duration="24H",
    ):
        """TODO update docstring


        Find the leading observation (last ok-record before the gap).

        Get the last observation before the gap, for which an observation exists.

        It is also possible to use the max_lead_to_gap_distance to specify the
        maximum distantance (in time) between the start of the gap and the leading
        record.

        Parameters
        ----------
        observation_series : pandas.Series
            The observation values for the specific station and obstype. This
            is thus a pandas.Series with datetime as index. Outliers can be
            represented by nan's.
        max_lead_to_gap_distance : Timedelta or str, optional
            The max distance, in time, between the gap start and leading record.

        Returns
        -------
        tuple : (name, timestamp, obs-value)
            All information on the leading record stored in a tuple.

        """
        # sta_obs = xs_save(Dataset.df, self.name, level="name")
        # sta_obs = sta_obs[[self.obstype.name]]  # subset to the gap variable
        sta_obs = observations_series.dropna()  # remove the Nans (from qc)
        sta_obs_records = sta_obs.index

        period_start = self.startdt - pd.Timedelta(leading_period_duration)
        leading_period = sta_obs_records[
            ((sta_obs_records < self.startdt) & (sta_obs_records >= period_start))
        ]

        if leading_period.empty:  # No records found
            return tuple(
                [
                    self.name,
                    pd.DatetimeIndex([]),
                    np.nan,
                    # "no leading record candidates could be found.",
                ]
            )

        else:
            return tuple(
                [
                    self.name,
                    leading_period,
                    sta_obs.loc[leading_period].to_list(),
                    # "ok",
                ]
            )

    def get_trailing_period(self, observations_series, trailing_period_duration="24H"):
        """TODO update docstring


        Find the leading observation (last ok-record before the gap).

        Get the last observation before the gap, for which an observation exists.

        It is also possible to use the max_lead_to_gap_distance to specify the
        maximum distantance (in time) between the start of the gap and the leading
        record.

        Parameters
        ----------
        observation_series : pandas.Series
            The observation values for the specific station and obstype. This
            is thus a pandas.Series with datetime as index. Outliers can be
            represented by nan's.
        max_lead_to_gap_distance : Timedelta or str, optional
            The max distance, in time, between the gap start and leading record.

        Returns
        -------
        tuple : (name, timestamp, obs-value)
            All information on the leading record stored in a tuple.

        """
        # sta_obs = xs_save(Dataset.df, self.name, level="name")
        # sta_obs = sta_obs[[self.obstype.name]]  # subset to the gap variable
        sta_obs = observations_series.dropna()  # remove the Nans (from qc)
        sta_obs_records = sta_obs.index

        period_end = self.enddt + pd.Timedelta(trailing_period_duration)
        trailing_period = sta_obs_records[
            ((sta_obs_records > self.enddt) & (sta_obs_records <= period_end))
        ]

        if trailing_period.empty:  # No records found
            return tuple(
                [
                    self.name,
                    pd.DatetimeIndex([]),
                    [],
                    # "no trailing record candidates could be found.",
                ]
            )

        else:
            return tuple(
                [
                    self.name,
                    trailing_period,
                    sta_obs.loc[trailing_period].to_list(),
                    # "ok",
                ]
            )

    def get_leading_record(self, observations_series, max_lead_to_gap_distance=None):
        """Find the leading observation (last ok-record before the gap).

        Get the last observation before the gap, for which an observation exists.

        It is also possible to use the max_lead_to_gap_distance to specify the
        maximum distantance (in time) between the start of the gap and the leading
        record.

        Parameters
        ----------
        observation_series : pandas.Series
            The observation values for the specific station and obstype. This
            is thus a pandas.Series with datetime as index. Outliers can be
            represented by nan's.
        max_lead_to_gap_distance : Timedelta or str, optional
            The max distance, in time, between the gap start and leading record.

        Returns
        -------
        tuple : (name, timestamp, obs-value, status-message)
            All information on the leading record stored in a tuple.

        """
        # sta_obs = xs_save(Dataset.df, self.name, level="name")
        # sta_obs = sta_obs[[self.obstype.name]]  # subset to the gap variable
        sta_obs = observations_series.dropna()  # remove the Nans (from qc)
        sta_obs_records = sta_obs.index
        leading_timestamp = sta_obs_records[sta_obs_records < self.startdt].max()
        if pd.isnull(leading_timestamp):
            return tuple(
                [
                    self.name,
                    pd.Timestamp("NaT").to_pydatetime(),
                    np.nan,
                    "no leading record candidate could be found.",
                ]
            )

        if max_lead_to_gap_distance is not None:
            delta_t = pd.Timedelta(max_lead_to_gap_distance)
            if self.startdt - leading_timestamp > delta_t:
                return tuple(
                    [
                        self.name,
                        pd.Timestamp("NaT").to_pydatetime(),
                        np.nan,
                        "leading record distance to gap is to large.",
                    ]
                )

        leading_value = sta_obs.loc[leading_timestamp]
        return tuple([self.name, leading_timestamp, leading_value, "ok"])

    def get_trailing_record(self, observations_series, max_trail_to_gap_distance=None):
        """Find the trailing observation (first ok-record after gap).

        Get the first observation after the gap, for which an observation exists.

        It is also possible to use the max_trail_to_gap_distance to specify the
        maximum distantance (in time) between the end of the gap and the trailing
        record.

        Parameters
        ----------
        observation_series : pandas.Series
            The observation values for the specific station and obstype. This
            is thus a pandas.Series with datetime as index. Outliers can be
            represented by nan's.
        max_trail_to_gap_distance : Timedelta or str, optional
            The max distance, in time, between the gap end and trailing record.

        Returns
        -------
        tuple : (name, timestamp, obs-value, status-message)
            All information on the trailing record stored in a tuple.

        """
        # sta_obs = xs_save(Dataset.df, self.name, level="name")
        # sta_obs = sta_obs[[self.obstype.name]]  # subset to the gap variable
        sta_obs = observations_series.dropna()  # remove the Nans (from qc)
        sta_obs_records = sta_obs.index
        trailing_timestamp = sta_obs_records[sta_obs_records > self.enddt].min()
        if pd.isnull(trailing_timestamp):
            return tuple(
                [
                    self.name,
                    pd.Timestamp("NaT").to_pydatetime(),
                    np.nan,
                    "no trailing record candidate could be found.",
                ]
            )

        if max_trail_to_gap_distance is not None:
            delta_t = pd.Timedelta(max_trail_to_gap_distance)
            if trailing_timestamp - self.enddt > delta_t:
                return tuple(
                    [
                        self.name,
                        pd.Timestamp("NaT").to_pydatetime(),
                        np.nan,
                        "trailing record distance to gap is to large.",
                    ]
                )

        trailing_value = sta_obs.loc[trailing_timestamp]
        return tuple([self.name, trailing_timestamp, trailing_value, "ok"])


# =============================================================================
# Find gaps and missing values
# =============================================================================
def get_station_gaps(gapslist, name):
    """Extract a Gap_collection specific to one station.

    If no gaps are found for the station, an empty Gap_collection is
    returned.

    Parameters
    ----------
    name : String
        Name of the station to extract a Gaps_collection from.

    Returns
    -------
    Gap_collection
        A Gap collection specific of the specified station.

    """
    return [gap for gap in gapslist if gap.name == name]


# =============================================================================
# Helpers
# =============================================================================


# =============================================================================
# Gap finders
# =============================================================================
def find_gaps(df, metadf, outliersdf, obstypes):

    gap_list = []
    # Replace the nans in the df, which are assigned as outliers, to a default value
    df.loc[outliersdf.index, "value"] = (
        -999
    )  # WARNING! df is a pointer, so this affects the self.df of the dataset !!!
    # make sure to revert this!!

    # Subset to the missing records
    missing_records = df[df["value"].isnull()]

    # Group by station and obstype to locate and define gaps
    for groupidx, groupdf in missing_records.reset_index().groupby(["name", "obstype"]):
        staname = groupidx[0]
        obsname = groupidx[1]
        sta_res = metadf.loc[staname, "dataset_resolution"]

        # calculate differences (at level of obstype) because else
        # gaps are computed wrong !! (DO NOT TRY TO REDUCE COMPUTATIONS ON THIS STEP! )
        groupdf["diff"] = groupdf["datetime"] - groupdf["datetime"].shift()

        # find groups of consecutive missing records
        groupdf["gap_def"] = ((groupdf["diff"] != sta_res)).cumsum()
        gapgroups = groupdf.groupby("gap_def")

        for _idx, gapgroup in gapgroups:
            # Construct gap
            gap = Gap(
                name=staname,
                startdt=gapgroup["datetime"].min(),
                enddt=gapgroup["datetime"].max(),
                obstype=obstypes[obsname],
                records_freq=sta_res,
            )

            # # Compute missing records
            # missing_records = pd.date_range(
            #     start=gap.startdt, end=gap.enddt, freq=sta_res, tz=gap.startdt.tz
            # )
            # gap._construct_gapdf(missing_datetime_records=missing_records)

            # add gap to list
            gap_list.append(gap)

    # Make sure to revert all outlier values to Nan (because they are present as -999)
    df.loc[outliersdf.index, "value"] = np.nan

    return gap_list


# =============================================================================
# Errors
# =============================================================================


class MetobsGapError(Exception):
    """Exception raised for errors in the Gap() class"""

    pass
