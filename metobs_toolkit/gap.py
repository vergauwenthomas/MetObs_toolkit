#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module contains the Gap class and all its methods.

A Gap contains all information and methods of a data-gap.
"""

import pandas as pd
import numpy as np
import logging
from datetime import timedelta
import math
import sys
import datetime


import metobs_toolkit.gap_filling as gap_filling

# from metobs_toolkit.gap_filling import (
#     interpolate_gap,
#     create_leading_trailing_debias_periods,
#     make_era_bias_correction,
# )

from metobs_toolkit.df_helpers import (
    init_multiindexdf,
    get_likely_frequency,
    xs_save,
)

# from metobs_toolkit.df_helpers import (
#     format_outliersdf_to_doubleidx,
#     concat_save,
#     get_likely_frequency,
#     _find_closes_occuring_date,
# )


# from metobs_toolkit.df_helpers import init_multiindex, xs_save

# from metobs_toolkit.missingobs import Missingob_collection
from metobs_toolkit.obstypes import Obstype as Obstype_class

logger = logging.getLogger(__name__)


# =============================================================================
# Gap class

# a gap is a sequence of repeting missing obs
# =============================================================================


class Gap:
    """Gap class holds all gap information and methods for gaps."""

    def __init__(self, name, startdt, enddt, obstype):
        """
        Initiate Gap object with a name, startdt and enddt.

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

        # Construct derived attributes and initiate them

        # gapdf holds all the records in the gap (depending on the freq of the dataset!)
        self.gapdf = (
            init_multiindexdf()
        )  # columns: obstype.name, obstype.name"_fill', "fill_method", "msg"
        # anchordf holds all the records outside the gap, used by fill methods.
        self.anchordf = (
            init_multiindexdf()
        )  # columns: obstype.name, obstype.name"_fill', "fill_method", "msg"

    def __str__(self):
        return f"{self.obstype.name}-gap of {self.name} for {self.startdt} --> {self.enddt}, duration: {self.duration}"

    def __repr__(self):
        return self.__str__()

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
        if not self.does_gap_holds_missing_records():
            print(
                "The gap is unreal (=no missing records at the current frequency resolution)."
            )
            return

        if not self.anchordf.empty:
            print("---- Anchor Data Frame ----- ")
            print(self.anchordf)

        print("---- Gap Data Frame -----")
        print(self.gapdf)

    # =============================================================================
    # Helpers
    # =============================================================================
    def _get_gapfill_status(self):
        # Warning! do not change these labels (or replace all hardcode occurences)
        if not self.does_gap_holds_missing_records():
            return "Gap does not exist in observation space"
        if self.gapdf[f"{self.obstype.name}_fill"].isnull().all():
            return "Unfilled gap"
        if self.gapdf[f"{self.obstype.name}_fill"].isnull().any():
            return "Partially filled gap"
        else:
            return "Filled gap"

    def does_gap_holds_missing_records(self):
        if self.gapdf.empty:
            return False
        else:
            return True

    def _initiate_gapdf(self, Dataset):
        """Initiate the gapdf attribute (which is dep on the observation freq)"""

        logger.debug(f"initiate gapdf for {self}")

        if not self.gapdf.empty:
            logger.debug(f"The gapdf of {self} will be overwritten and initialized")
            # print(f"The gapdf of {self} will be overwritten and initialized")
            self.gapdf = init_multiindexdf()

        # Get the missing records in the observationspace-frequency
        dtrecords = self.get_missing_records(Dataset)

        # Create index
        records_multiidx = pd.MultiIndex.from_arrays(
            arrays=[[self.name] * len(dtrecords), dtrecords], names=["name", "datetime"]
        )
        self.gapdf.index = records_multiidx

        # Set columns with default values
        self.gapdf[f"{self.obstype.name}"] = np.nan  # fill with nan values
        self.gapdf[f"{self.obstype.name}_fill"] = np.nan  # fill with nan values
        self.gapdf["fill_method"] = None  # fill with None
        self.gapdf["msg"] = None  # Nonedataset

    # =============================================================================
    # Gap --- Records methods
    # =============================================================================

    def get_missing_records(self, Dataset):
        """Find the missing records of a gap in observation-space.

        Gaps are located (at import frequency) and the period is define by
        a start and end of the gap (thus not in records!!). This mehod will
        find the missing records (timestamps), for the gap given an
        an observation space (for the time-resolution)

        Parameters
        ----------
        Dataset : metobs_toolkit.Dataset
            The dataset with the observationspace to create missing records for.

        Returns
        -------
        missing_records : pd.datetimeindex.
            DESCRIPTION.

        """
        logger.debug(f"Get missing records for {self}")
        assert not Dataset.df.empty, f"{Dataset} is empty."
        assert self.name in Dataset.df.index.get_level_values(
            "name"
        ), f"{self.name} not found in the observation records of {Dataset}"

        # filter out start, end and freq of observations
        obsfreq = Dataset.metadf.loc[self.name, "dataset_resolution"]
        obsstart = Dataset.df.xs(self.name, level="name").index.min()
        obsend = Dataset.df.xs(self.name, level="name").index.max()

        # create ideal datetimerange
        obs_records_ideal = pd.date_range(start=obsstart, end=obsend, freq=obsfreq)

        # Subset to the gap-period (so only missing records are kept)
        missing_records = obs_records_ideal[
            (obs_records_ideal <= self.enddt) & ((obs_records_ideal >= self.startdt))
        ]
        return missing_records

    # =============================================================================
    #  Gapfill anchor methods
    # =============================================================================
    def get_leading_period(
        self,
        Dataset,
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
        Dataset : metobs_toolkit.Dataset
            The dataset with the observations to look for the leading record.
        max_lead_to_gap_distance : Timedelta or str, optional
            The max distance, in time, between the gap start and leading record.

        Returns
        -------
        tuple : (name, timestamp, obs-value)
            All information on the leading record stored in a tuple.

        """
        sta_obs = xs_save(Dataset.df, self.name, level="name")
        sta_obs = sta_obs[[self.obstype.name]]  # subset to the gap variable
        sta_obs = sta_obs.dropna()  # remove the Nans (from qc)
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
                    sta_obs.loc[leading_period][self.obstype.name].to_list(),
                    # "ok",
                ]
            )

    def get_trailing_period(self, Dataset, trailing_period_duration="24H"):
        """TODO update docstring


        Find the leading observation (last ok-record before the gap).

        Get the last observation before the gap, for which an observation exists.

        It is also possible to use the max_lead_to_gap_distance to specify the
        maximum distantance (in time) between the start of the gap and the leading
        record.

        Parameters
        ----------
        Dataset : metobs_toolkit.Dataset
            The dataset with the observations to look for the leading record.
        max_lead_to_gap_distance : Timedelta or str, optional
            The max distance, in time, between the gap start and leading record.

        Returns
        -------
        tuple : (name, timestamp, obs-value)
            All information on the leading record stored in a tuple.

        """
        sta_obs = xs_save(Dataset.df, self.name, level="name")
        sta_obs = sta_obs[[self.obstype.name]]  # subset to the gap variable
        sta_obs = sta_obs.dropna()  # remove the Nans (from qc)
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
                    sta_obs.loc[trailing_period][self.obstype.name].to_list(),
                    # "ok",
                ]
            )

    def get_leading_record(self, Dataset, max_lead_to_gap_distance=None):
        """Find the leading observation (last ok-record before the gap).

        Get the last observation before the gap, for which an observation exists.

        It is also possible to use the max_lead_to_gap_distance to specify the
        maximum distantance (in time) between the start of the gap and the leading
        record.

        Parameters
        ----------
        Dataset : metobs_toolkit.Dataset
            The dataset with the observations to look for the leading record.
        max_lead_to_gap_distance : Timedelta or str, optional
            The max distance, in time, between the gap start and leading record.

        Returns
        -------
        tuple : (name, timestamp, obs-value, status-message)
            All information on the leading record stored in a tuple.

        """
        sta_obs = xs_save(Dataset.df, self.name, level="name")
        sta_obs = sta_obs[[self.obstype.name]]  # subset to the gap variable
        sta_obs = sta_obs.dropna()  # remove the Nans (from qc)
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

        leading_value = sta_obs.loc[leading_timestamp, self.obstype.name]
        return tuple([self.name, leading_timestamp, leading_value, "ok"])

    def get_trailing_record(self, Dataset, max_trail_to_gap_distance=None):
        """Find the trailing observation (first ok-record after gap).

        Get the first observation after the gap, for which an observation exists.

        It is also possible to use the max_trail_to_gap_distance to specify the
        maximum distantance (in time) between the end of the gap and the trailing
        record.

        Parameters
        ----------
        Dataset : metobs_toolkit.Dataset
            The dataset with the observations to look for the leading record.
        max_trail_to_gap_distance : Timedelta or str, optional
            The max distance, in time, between the gap end and trailing record.

        Returns
        -------
        tuple : (name, timestamp, obs-value, status-message)
            All information on the trailing record stored in a tuple.

        """
        sta_obs = xs_save(Dataset.df, self.name, level="name")
        sta_obs = sta_obs[[self.obstype.name]]  # subset to the gap variable
        sta_obs = sta_obs.dropna()  # remove the Nans (from qc)
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

        trailing_value = sta_obs.loc[trailing_timestamp, self.obstype.name]
        return tuple([self.name, trailing_timestamp, trailing_value, "ok"])

    # =============================================================================
    # Fill methods
    # =============================================================================

    def raw_model_gapfill(self, Dataset, Modeldata):
        obsname = self.obstype.name

        # add the gap period
        gapdf = Modeldata.interpolate_modeldata(self.gapdf.index)
        gapdf = gapdf[[obsname]]
        gapdf = gapdf.rename(columns={obsname: "modelvalues"})
        gapdf["obsvalues"] = np.nan
        gapdf["fill_method"] = "gap"

        filldf = gapdf
        filldf = filldf.sort_index()

        # 4. Fill the gap
        fill_values, msg = gap_filling._raw_modeldata_fill(filldf)

        # 5. Update gapdf attribute
        self.gapdf[f"{obsname}_fill"] = fill_values
        self.gapdf["fill_method"] = "raw_modeldata"
        self.gapdf["msg"] = msg

    def debias_model_gapfill(
        self,
        Dataset,
        Modeldata,
        leading_period_duration="24H",
        min_leading_records_total=60,
        trailing_period_duration="24H",
        min_trailing_records_total=60,
    ):

        obsname = self.obstype.name
        self.gapdf["fill_method"] = "debias_modeldata"
        # Create anchor periods
        anchordf = self._create_anchor_df_for_leading_trailing_periods(
            Dataset, leading_period_duration, trailing_period_duration
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
        self.anchordf = anchordf

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

        if err:
            logger.warning(logmsg)
            print("Warning! ", logmsg)
            self.gapdf[f"{obsname}_fill"] = np.nan
            self.gapdf["msg"] = gapmsg
            return

        # 4. combine learning and gap period
        filldf = self._combine_learning_and_gap_to_one_df(Modeldata)

        # 5. Fill the gap
        fill_values, msg = gap_filling._debias_modeldata_fill(filldf)

        # 6. Update gapdf attribute
        self.gapdf[f"{obsname}_fill"] = fill_values
        self.gapdf["fill_method"] = "debias_modeldata"
        self.gapdf["msg"] = msg

    def diurnal_debias_model_gapfill(
        self,
        Dataset,
        Modeldata,
        leading_period_duration="24H",
        min_debias_sample_size=6,
        trailing_period_duration="24H",
    ):

        obsname = self.obstype.name
        # Create anchor periods
        anchordf = self._create_anchor_df_for_leading_trailing_periods(
            Dataset, leading_period_duration, trailing_period_duration
        )
        # Check validity of anchors and label them
        anchordf = self._label_anchors_for_diurnal_gapfilling(
            anchordf=anchordf, min_anchors_per_diurnal_timestamp=min_debias_sample_size
        )

        self.anchordf = anchordf

        # 4. combine learning and gap period
        filldf = self._combine_learning_and_gap_to_one_df(Modeldata)

        # 5. apply blacklisting (remove anchors that do not fulfill criteria)
        blacklist = self.anchordf[self.anchordf["msg"] != "ok"].index
        filldf = filldf[~filldf.index.isin(blacklist)]

        # 6. Fill the gap
        fill_values, msg = gap_filling._diurnal_debias_modeldata_fill(filldf)

        # 7. Update gapdf attribute
        self.gapdf[f"{obsname}_fill"] = fill_values
        self.gapdf["fill_method"] = "diurnal_debias_modeldata"
        self.gapdf["msg"] = msg

    def weighted_diurnal_debias_model_gapfill(
        self,
        Dataset,
        Modeldata,
        leading_period_duration="48H",
        min_lead_debias_sample_size=3,
        trailing_period_duration="48H",
        min_trail_debias_sample_size=3,
    ):

        obsname = self.obstype.name
        # Create anchor periods
        anchordf = self._create_anchor_df_for_leading_trailing_periods(
            Dataset, leading_period_duration, trailing_period_duration
        )

        # Check validity of anchors and label them
        lead_anchorsdf = anchordf[anchordf["fill_method"] == "leading period"]
        lead_anchorsdf = self._label_anchors_for_diurnal_gapfilling(
            anchordf=lead_anchorsdf,
            min_anchors_per_diurnal_timestamp=min_lead_debias_sample_size,
        )

        trail_anchorsdf = anchordf[anchordf["fill_method"] == "trailing period"]
        trail_anchorsdf = self._label_anchors_for_diurnal_gapfilling(
            anchordf=trail_anchorsdf,
            min_anchors_per_diurnal_timestamp=min_trail_debias_sample_size,
        )

        anchordf = pd.concat([lead_anchorsdf, trail_anchorsdf])

        self.anchordf = anchordf

        # 4. combine learning and gap period
        filldf = self._combine_learning_and_gap_to_one_df(Modeldata)

        # 5. apply blacklisting (remove anchors that do not fulfill criteria)
        blacklist = self.anchordf[self.anchordf["msg"] != "ok"].index
        filldf = filldf[~filldf.index.isin(blacklist)]

        # 6. Fill the gap
        fill_values, msg = gap_filling._weighted_diurnal_debias_modeldata(filldf)

        # 7. Update gapdf attribute
        self.gapdf[f"{obsname}_fill"] = fill_values
        self.gapdf["fill_method"] = "weighted_diurnal_debias_modeldata"
        self.gapdf["msg"] = msg

    def interpolate_gap(
        self,
        Dataset,
        method="time",
        max_consec_fill=10,
        max_lead_to_gap_distance=None,
        max_trail_to_gap_distance=None,
    ):

        obsname = self.obstype.name
        # 1. Get leading and trailing info
        # get leading record, check validity and add to the gapfilldf
        (_, lead_dt, lead_val, lead_msg) = self.get_leading_record(
            Dataset=Dataset, max_lead_to_gap_distance=max_lead_to_gap_distance
        )
        (_, trail_dt, trail_val, trail_msg) = self.get_trailing_record(
            Dataset=Dataset, max_trail_to_gap_distance=max_trail_to_gap_distance
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
        self.anchordf = anchor_df
        self.gapdf["fill_method"] = "interpolation"

        # check if gapfill can proceed
        if (lead_msg != "ok") & (trail_msg != "ok"):
            logger.warning(
                f"Cannot fill {self}, because leading and trailing records are not valid."
            )
            print(
                f"Warning! Cannot fill {self}, because leading and trailing periods are not valid."
            )

            self.gapdf[f"{obsname}_fill"] = np.nan
            self.gapdf["msg"] = f"{lead_msg} and {trail_msg}"
            return
        elif (lead_msg != "ok") & (trail_msg == "ok"):
            logger.warning(f"Cannot fill {self}, because leading record is not valid.")
            print(f"Warning! Cannot fill {self}, because leading record is not valid.")

            self.gapdf[f"{obsname}_fill"] = np.nan
            self.gapdf["msg"] = f"{lead_msg}"
            return
        elif (lead_msg == "ok") & (trail_msg != "ok"):
            logger.warning(f"Cannot fill {self}, because trailing record is not valid.")
            print(f"Warning! Cannot fill {self}, because trailing record is not valid.")

            self.gapdf[f"{obsname}_fill"] = np.nan
            self.gapdf["msg"] = f"{trail_msg}"
            return
        else:
            pass

        # 3. Combine the anchors with the observations
        tofilldf = pd.concat([self.gapdf[[obsname]], self.anchordf[[obsname]]])

        # 4. Replace the NaN's (GAPFILLING)
        # make shure only datetimes are in the index, and sorted by datetimes
        tofilldf = tofilldf.reset_index().set_index("datetime").sort_index()
        # Interpolate series
        tofilldf[obsname].interpolate(
            method=method,
            limit=max_consec_fill,  # Maximum number of consecutive NaNs to fill. Must be greater than 0.
            limit_area="inside",
            inplace=True,
        )
        # revert to multiindex
        tofilldf = tofilldf.reset_index().set_index(["name", "datetime"]).sort_index()

        # 5. Get messages for all filled records
        # filter to gap records
        tofilldf = tofilldf.loc[self.gapdf.index]

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
            sys.exit("*unforseen situation* Report this error!!")

        # 6. Update the gapdf
        self.gapdf[f"{obsname}_fill"] = tofilldf[obsname]  # update values

        self.gapdf["msg"] = tofilldf["msg"]

        return

    # =============================================================================
    # Fill helpers
    # =============================================================================

    def _create_anchor_df_for_leading_trailing_periods(
        self, Dataset, leading_period_duration, trailing_period_duration
    ):
        """
        Helper method to construct the anchordf (the dataframe with all anchor
        records) for a leading and trailing period.

        Parameters
        ----------
        Dataset : metobs_toolkit.Dataset
            The dataset that contains the observations for the anchor records.
        leading_period_duration : Timedelta or str
            Size (in time) of the leading period.
        trailing_period_duration : Timedelta or str
            Size (in time) of the trailing period.

        Returns
        -------
        anchor_df : pandas.dataframe
            The anchordf (multiindex name-timestamp) with a column with values
            and a "fill_method" column with labels ('leading period' or
            'trailing period').

        """
        # 1. Get leading and trailing info
        # get leading record, check validity and add to the gapfilldf
        (_, lead_period, lead_vals) = self.get_leading_period(
            Dataset=Dataset,
            leading_period_duration=leading_period_duration,
        )

        (_, trail_period, trail_vals) = self.get_trailing_period(
            Dataset=Dataset,
            trailing_period_duration=trailing_period_duration,
        )

        # 2. Create anchordf
        _leaddf_idx = pd.MultiIndex.from_arrays(
            arrays=[[self.name] * lead_period.shape[0], lead_period],
            names=["name", "datetime"],
        )

        _leaddf = pd.DataFrame(
            data={
                self.obstype.name: lead_vals,
                "fill_method": ["leading period"],
                # "msg": lead_msg,
            },
            index=_leaddf_idx,
        )

        _traildf_idx = pd.MultiIndex.from_arrays(
            arrays=[[self.name] * trail_period.shape[0], trail_period],
            names=["name", "datetime"],
        )

        _traildf = pd.DataFrame(
            data={
                self.obstype.name: trail_vals,
                "fill_method": ["trailing period"],
                # "msg": trail_msg,
            },
            index=_traildf_idx,
        )

        anchor_df = pd.concat([_leaddf, _traildf]).sort_index()
        return anchor_df

    def _combine_learning_and_gap_to_one_df(self, Modeldata):
        """
        Helper method to combine the anchorsdf and the gapdf to one dataframe.
        This combined dataframe is than used by the fill method.

        Modeldata is extracted for all these records (and interpolated to the
        to match time resolution).

        Each record is labeld by 'leading period', 'trailing period' or 'gap'

        Parameters
        ----------
        Modeldata : metobs_toolkit.Modeldata
            The modeldata to use for the gapfilling

        Returns
        -------
        filldf : pandas.DataFrame
            The dataframe structured as the Gap.gapdf, that combines the
            anchors and the gap records and the corresponding modelvalues.

        """
        anchordf = self.anchordf
        obsname = self.obstype.name

        debiasdf = Modeldata.interpolate_modeldata(anchordf.index)
        assert (
            obsname in debiasdf.columns
        ), f"{obsname} not present in the modeldata: {Modeldata}"
        debiasdf = debiasdf[[obsname]]
        debiasdf = debiasdf.rename(columns={obsname: "modelvalues"})
        debiasdf["obsvalues"] = anchordf[obsname]
        debiasdf["fill_method"] = anchordf["fill_method"]

        # add the gap period
        gapdf = Modeldata.interpolate_modeldata(self.gapdf.index)
        gapdf = gapdf[[obsname]]
        gapdf = gapdf.rename(columns={obsname: "modelvalues"})
        gapdf["obsvalues"] = np.nan
        gapdf["fill_method"] = "gap"

        filldf = pd.concat([debiasdf, gapdf])
        filldf = filldf.sort_index()
        return filldf

    def _label_anchors_for_diurnal_gapfilling(
        self, anchordf, min_anchors_per_diurnal_timestamp
    ):
        """
        Helper method to label (in the msg column) all records of the anchorsdf.
        This methods is applied when filling uses diurnal timestamps, thus the
        'training period' depends on similar timestamps. This method checks if
        all these periods have equal, or more records then set by min_anchors_per_diurnal_timestamp.

        The following labels can be written:
         * 'ok'
         * 'will not be used (diurnal timestamp not in gap)'
         * 'diurnal sampel size to small ...

        Parameters
        ----------
        anchordf : pandas.DataFrame
            The dataframe containing all the anchorrecords to check.
        min_anchors_per_diurnal_timestamp : int
            The minimum size for each timestamp group to check.

        Returns
        -------
        anchordf : pandas.DataFrame
            The anchordf with the 'msg' (string message on status -> label)
            column and 'anchors' (count of records per group in the anchorsdf)

        """

        obsname = self.obstype.name
        anchordf["msg"] = "_init_"  # these will be overwritten
        # aggregate to diurnal groups and count number of samples For all anchors
        anchordf = anchordf.dropna()
        anchordf = anchordf.reset_index().set_index("datetime")
        anchordf = gap_filling._add_diurnal_timestamps(anchordf)
        anchordf = (
            anchordf.reset_index()
        )  # else the datetimes will be lost when mergeing

        samplesizes = anchordf.groupby(["hour", "minutes", "seconds"])[obsname].count()
        samplesizes.name = "anchors"

        # aggregate the gaps records to diurnal groups, to find the diurnal timestamps that will be used
        gap_records = self.gapdf.reset_index().set_index("datetime")
        gap_records = gap_filling._add_diurnal_timestamps(gap_records)
        gap_records = gap_records.groupby(["hour", "minutes", "seconds"])[
            obsname
        ].count()

        # label the anchors
        # A. 'will not be used' label
        not_used = samplesizes[~samplesizes.index.isin(gap_records.index)]
        not_useddf = not_used.to_frame().reset_index()
        not_useddf["msg_not_used"] = "_found"  # dummy label

        anchordf = anchordf.merge(
            not_useddf[["hour", "minutes", "seconds", "msg_not_used"]],
            how="left",
            on=["hour", "minutes", "seconds"],
        )

        anchordf.loc[~anchordf["msg_not_used"].isnull(), "msg"] = (
            "will not be used (diurnal timestamp not in gap)"
        )

        # B. "sample size to small
        to_small_diurnal_stamps = samplesizes[
            (
                (~samplesizes.index.isin(not_used.index))
                & (samplesizes < min_anchors_per_diurnal_timestamp)
            )
        ]

        to_smalldf = to_small_diurnal_stamps.to_frame().reset_index()
        to_smalldf["msg_to_small"] = to_smalldf.apply(
            lambda x: f'diurnal sampel size to small ({x["anchors"]} < {min_anchors_per_diurnal_timestamp})',
            axis=1,
        )

        anchordf = anchordf.merge(
            to_smalldf[["hour", "minutes", "seconds", "msg_to_small"]],
            how="left",
            on=["hour", "minutes", "seconds"],
        )

        anchordf.loc[~anchordf["msg_to_small"].isnull(), "msg"] = anchordf.loc[
            ~anchordf["msg_to_small"].isnull(), "msg_to_small"
        ]

        # C "ok"  the remaining anchors fulfill the requirements
        anchordf["msg"] = anchordf["msg"].replace({"_init_": "ok"})

        # Subset anchordf to standard format
        anchordf = anchordf.reset_index().set_index(["name", "datetime"])
        anchordf = anchordf[[obsname, "fill_method", "msg"]]

        return anchordf


def create_gaps_overview_df(gapslist):
    """TODO: docstring"""
    gapsdf_list = []
    for gap in gapslist:
        gapdf = gap.gapdf
        gapdf["obstype"] = gap.obstype.name
        gapdf = gapdf.rename(columns={f"{gap.obstype.name}_fill": "fill_value"})
        gapsdf_list.append(gapdf[["fill_value", "obstype", "fill_method", "msg"]])

    tot_gapdf = pd.concat(gapsdf_list)
    tot_gapdf = tot_gapdf.reset_index().set_index(["name", "datetime", "obstype"])
    return tot_gapdf.sort_index()


# =============================================================================
# Gap finders
# =============================================================================


def find_gaps(df, known_obstypes, tstart=None, tend=None, freq_series=None):
    """
    #TODO update the docstring parameters
    Find gaps in the observations.

    Looking for gaps by assuming an observation frequency. The assumed frequency is the highest occuring frequency PER STATION.

    A gap in the observations is located by two conditions:

     * A timestamp is missing
     * The value of an observation is invallid (i.e. Nan or non-numerical)



    Parameters
    ----------
    df : pandas.DataFrame
        The observations dataframe of the dataset object (Dataset.df)
    known_obstypes : list
        A list of all known metobs_toolkit.Obstypes.


    Returns
    -------
    gap_list : list
        A list of gaps.

    """
    gap_list = []

    # get all obstypes that are found in the dataframe
    found_obstypes = [obst for obst in known_obstypes if obst.name in df.columns]

    # gaps per station (because some stations can have other frequencies!)
    stationnames = df.index.get_level_values(level="name").unique()
    for station in stationnames:

        stadf = xs_save(df, station, level="name")

        # Get frequency estimation per station
        if freq_series is None:
            likely_freq = get_likely_frequency(
                stadf.index, method="highest", simplify=False
            )
        else:
            likely_freq = freq_series[station]

        # Define the start and end of the timeseries
        if tstart is None:
            tstart = stadf.index.min()
        if tend is None:
            tend = stadf.index.max()

        assert likely_freq.seconds > 0, "The frequency is not positive!"
        for obstype in found_obstypes:
            # remove Nans --> so that nans are interpretted as missing records --> thus gaps
            stadfnonan = stadf.dropna(subset=[obstype.name])[[obstype.name]]
            #  Locate missing timestamps
            missing_datetimeseries = (
                pd.date_range(start=tstart, end=tend, freq=likely_freq)
                .difference(stadfnonan.index)
                .to_series()
                .diff()
            )

            if missing_datetimeseries.empty:
                continue

            # Check for gaps
            gap_defenition = ((missing_datetimeseries != likely_freq)).cumsum()
            consec_missing_groups = missing_datetimeseries.groupby(gap_defenition)
            group_sizes = consec_missing_groups.size()

            # iterate over the gabs and fill the gap_list
            for gap_idx in group_sizes.index:
                datetime_of_gap_records = consec_missing_groups.get_group(gap_idx).index
                gap_list.append(
                    Gap(
                        name=station,
                        startdt=datetime_of_gap_records.min(),
                        enddt=datetime_of_gap_records.max(),
                        obstype=obstype,
                    )
                )

    return gap_list
