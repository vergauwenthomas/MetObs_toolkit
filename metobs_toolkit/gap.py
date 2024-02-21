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
                "The gap holds no missing records at the current frequency resolution."
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
            logger.warning(f"The gapdf of {self} will be overwritten and initialized")
            print(f"The gapdf of {self} will be overwritten and initialized")
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
        self, Dataset, leading_period_duration="24H", min_leading_records_total=60
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
        tuple : (name, timestamp, obs-value, status-message)
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
                    "no leading record candidates could be found.",
                ]
            )

        elif (
            leading_period.shape[0] < min_leading_records_total
        ):  # not enouggh (total) records found
            return tuple(
                [
                    self.name,
                    leading_period,
                    sta_obs.loc[leading_period][self.obstype.name].to_list(),
                    f"minimum records ({min_leading_records_total}) for leading period not met.",
                ]
            )

        else:
            return tuple(
                [
                    self.name,
                    leading_period,
                    sta_obs.loc[leading_period][self.obstype.name].to_list(),
                    "ok",
                ]
            )

    def get_trailing_period(
        self, Dataset, trailing_period_duration="24H", min_trailing_records_total=60
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
        tuple : (name, timestamp, obs-value, status-message)
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
                    "no trailing record candidates could be found.",
                ]
            )

        elif (
            trailing_period.shape[0] < min_trailing_records_total
        ):  # not enouggh (total) records found
            return tuple(
                [
                    self.name,
                    trailing_period,
                    sta_obs.loc[trailing_period][self.obstype.name].to_list(),
                    f"minimum records ({min_trailing_records_total}) for trailing period not met.",
                ]
            )

        else:
            return tuple(
                [
                    self.name,
                    trailing_period,
                    sta_obs.loc[trailing_period][self.obstype.name].to_list(),
                    "ok",
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

    def model_gapfill(
        self,
        Dataset,
        Modeldata,
        method="debias_modeldata",
        leading_period_duration="24H",
        min_leading_records_total=60,
        trailing_period_duration="24H",
        min_trailing_records_total=60,
    ):

        # check attributes
        assert method in [
            "raw_modeldata",
            "debias_modeldata",
            "diurnal_debias_modeldata",
            "period_diurnal_debias_modeldata",
        ], f"{method} is not a valid method."

        obsname = self.obstype.name
        # 1. Get leading and trailing info
        # get leading record, check validity and add to the gapfilldf
        (_, lead_period, lead_vals, lead_msg) = self.get_leading_period(
            Dataset=Dataset,
            leading_period_duration=leading_period_duration,
            min_leading_records_total=min_leading_records_total,
        )
        (_, trail_period, trail_vals, trail_msg) = self.get_trailing_period(
            Dataset=Dataset,
            trailing_period_duration=trailing_period_duration,
            min_trailing_records_total=min_trailing_records_total,
        )

        # 2. Update anchordf
        _leaddf_idx = pd.MultiIndex.from_arrays(
            arrays=[[self.name] * lead_period.shape[0], lead_period],
            names=["name", "datetime"],
        )

        _leaddf = pd.DataFrame(
            data={
                self.obstype.name: lead_vals,
                "fill_method": ["leading period"],
                "msg": lead_msg,
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
                "msg": trail_msg,
            },
            index=_traildf_idx,
        )

        anchor_df = pd.concat([_leaddf, _traildf]).sort_index()

        # 2b. Update attribute
        self.anchordf = anchor_df
        self.gapdf["fill_method"] = method

        # check if gapfill can proceed
        if method == "raw_modeldata":
            # leading and trailing are not used for fill
            pass
        else:
            if (lead_msg != "ok") & (trail_msg != "ok"):
                logger.warning(
                    f"Cannot fill {self}, because leading and trailing periods are not valid."
                )
                print(
                    f"Warning! Cannot fill {self}, because leading and trailing periods are not valid."
                )

                self.gapdf[f"{obsname}_fill"] = np.nan
                self.gapdf["msg"] = f"{lead_msg} and {trail_msg}"
                return
            elif (lead_msg != "ok") & (trail_msg == "ok"):
                logger.warning(
                    f"Cannot fill {self}, because leading period is not valid."
                )
                print(
                    f"Warning! Cannot fill {self}, because leading period is not valid."
                )

                self.gapdf[f"{obsname}_fill"] = np.nan
                self.gapdf["msg"] = f"{lead_msg}"
                return
            elif (lead_msg == "ok") & (trail_msg != "ok"):
                logger.warning(
                    f"Cannot fill {self}, because trailing period is not valid."
                )
                print(
                    f"Warning! Cannot fill {self}, because trailing period is not valid."
                )

                self.gapdf[f"{obsname}_fill"] = np.nan
                self.gapdf["msg"] = f"{trail_msg}"
                return
            else:
                pass

        # 3. extract modeldata for leading, trailing and gap period
        debiasdf = Modeldata.interpolate_modeldata(anchor_df.index)
        assert (
            obsname in debiasdf.columns
        ), f"{obsname} not present in the modeldata: {Modeldata}"
        debiasdf = debiasdf[[obsname]]
        debiasdf = debiasdf.rename(columns={obsname: "modelvalues"})
        debiasdf["obsvalues"] = anchor_df[obsname]
        debiasdf["fill_method"] = anchor_df["fill_method"]

        # add the gap period
        gapdf = Modeldata.interpolate_modeldata(self.gapdf.index)
        gapdf = gapdf[[obsname]]
        gapdf = gapdf.rename(columns={obsname: "modelvalues"})
        gapdf["obsvalues"] = np.nan
        gapdf["fill_method"] = "gap"

        filldf = pd.concat([debiasdf, gapdf])
        filldf = filldf.sort_index()

        # 4. Fill the gap
        if method == "raw_modeldata":
            fill_values, msg = gap_filling._raw_modeldata_fill(filldf)

        elif method == "debias_modeldata":
            fill_values, msg = gap_filling._debias_modeldata_fill(filldf)

        elif method == "diurnal_debias_modeldata":
            fill_values, msg = gap_filling._diurnal_debias_modeldata_fill(filldf)

        elif method == "period_diurnal_debias_modeldata":
            fill_values, msg = gap_filling._weighted_diurnal_debias_modeldata(filldf)

        else:
            sys.exit(f"Unknown method: {method}")

        # 5. Update gapdf attribute
        self.gapdf[f"{obsname}_fill"] = fill_values
        self.gapdf["fill_method"] = method
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
        self.gapdf["fill_method"] = f"{method}-interpolation"

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

        # if obs_only:
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
        #     self.leading_timestamp = self.startgap - timedelta(seconds=before_diff)
        #     self.trailing_timestamp = self.endgap + timedelta(seconds=after_diff)

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


# class Gap:
#     """Gap class holds all gap information and methods for gaps."""

#     def __init__(self, name, startdt, enddt):
#         """
#         Initiate Gap object with a name, startdt and enddt.

#         Parameters
#         ----------
#         name : String
#             Station name where the gap occures.
#         startdt : datetime.datetime
#             Start datetime of the gap (included).
#         enddt : datetime.datetime
#             End datetime of the gap (included)

#         Returns
#         -------
#         None

#         """
#         # init attributes
#         self.name = name
#         self.startgap = startdt  # in IO space
#         self.endgap = enddt  # in IO space
#         self.duration = enddt - startdt

#         # computed attributes
#         self.leading_timestamp = None  # last ob_dt before gap in datset space
#         self.leading_val = {}  # keys are obstypes
#         self.trailing_timestamp = None  # first ob_dt after gap in dataset space
#         self.trailing_val = {}  # keys are obstypes

#         self.exp_gap_idx = None

#         # gap fill (only for conventional saving)
#         self.gapfill_df = (
#             pd.DataFrame()
#         )  # index: datetime, columns: obstypes, values: fill_values
#         self.gapfill_technique = None  # will become a string
#         self.gapfill_info = (
#             None  # detailed infomation on the gapfill technique (only for the user)
#         )
#         self.gapfill_errormessage = {}  # keys are obstypes

#

#

#     def to_df(self):
#         """Convert a Gap object to a dataframe (with one row).

#         The station name is the index and two colums ('start_gap', 'end_gap')
#         are constructed.

#         Returns
#         -------
#         pandas.DataFrame()
#             Gap in dataframe format.

#         """
#         returndf = pd.DataFrame(
#             index=[self.name],
#             data={
#                 "start_gap": self.startgap,
#                 "end_gap": self.endgap,
#                 "duration": self.duration,
#             },
#         )
#         returndf.index.name = "name"
#         return returndf

#

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
#                                     observations.)
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

#     # =============================================================================
#     #         Gapfill
#     # =============================================================================

#

# # =============================================================================
# # Find gaps and missing values
# # =========================================================


# def get_gaps_indx_in_obs_space(gapslist, obsdf, outliersdf, resolutionseries):
#     """Get all gaps in obsspace.

#     Explode the gaps, to the dataset resolution and format to a multiindex
#     with name -- datetime.

#     In addition the last observation before the gap (leading), and first
#     observation (after) the gap are computed and stored in the df attribute.
#     (the outliers are used to look for leading and trailing observations.)


#     Parameters
#     ----------
#     obsdf : pandas.DataFrame
#         Dataframe containing all the observations.
#     outliersdf : pandas.DataFrame
#         Dataframe containing all outliers
#     resolutionseries : pandas.Series
#         The resolution of each station in a Series with the stationname as an index.

#     Returns
#     -------
#     expanded_gabsidx_obsspace : pandas.index
#         Multiindex with name and datetime of gaps in obsspace.

#     """
#     outliersdf = format_outliersdf_to_doubleidx(outliersdf)

#     expanded_gabsidx_obsspace = init_multiindex()

#     for gap in gapslist:
#         gap.update_gaps_indx_in_obs_space(
#             obsdf, outliersdf, resolutionseries.loc[gap.name]
#         )
#         expanded_gabsidx_obsspace = expanded_gabsidx_obsspace.append(gap.exp_gap_idx)

#     return expanded_gabsidx_obsspace


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


# def gaps_to_df(gapslist):
#     """Combine all gaps into a dataframe as an overview.

#     Parameters
#     ----------
#     gapslist : list
#         List of gaps.

#     Returns
#     -------
#     pandas.DataFrame
#         A DataFrame with stationnames as index, and the start, end and duretion
#         of the gaps as columns.

#     """
#     gapdflist = []
#     for gap in gapslist:
#         gapdflist.append(gap.to_df())

#     if not bool(gapdflist):
#         # when no gaps, make default return
#         default_df = pd.DataFrame(data={"start_gap": [], "end_gap": [], "duration": []})
#         default_df.index.name = "name"
#         return default_df

#     return concat_save(gapdflist)


# def remove_gaps_from_obs(gaplist, obsdf):
#     """
#     Remove station - datetime records that are in the gaps from the obsdf.

#     (Usefull when filling timestamps to a df, and if you whant to remove the
#       gaps.)

#     Parameters
#     ----------
#     obsdf : pandas.DataFrame()
#         A MultiIndex dataframe with name -- datetime as index.

#     Returns
#     -------
#     obsdf : pandas.DataFrame()
#         The same dataframe with records inside gaps removed.

#     """
#     # Create index for gaps records in the obsdf
#     expanded_gabsidx = init_multiindex()
#     for gap in gaplist:
#         sta_records = xs_save(obsdf, gap.name, level="name").index  # filter by name

#         gaps_dt = sta_records[
#             (sta_records >= gap.startgap)
#             & (sta_records <= gap.endgap)  # filter if the observations are within a gap
#         ]

#         gaps_multiidx = pd.MultiIndex.from_arrays(
#             arrays=[[gap.name] * len(gaps_dt), gaps_dt], names=["name", "datetime"]
#         )

#         expanded_gabsidx = expanded_gabsidx.append(gaps_multiidx)

#     # remove gaps idx from the obsdf
#     obsdf = obsdf.drop(index=expanded_gabsidx)
#     return obsdf


# def remove_gaps_from_outliers(gaplist, outldf):
#     """Remove station - datetime records that are in the gaps from the outliersdf.

#     This will ignore the observation types! So all outliers of any observation
#     type, that are in a gap period, are removed.

#     Parameters
#     ----------
#     obsdf : pandas.DataFrame()
#         A MultiIndex dataframe with name -- datetime -- as index.

#     Returns
#     -------
#     obsdf : pandas.DataFrame()
#         The same dataframe with records inside gaps removed.

#     """
#     # to multiindex
#     outldf = outldf.reset_index().set_index(["name", "datetime"])

#     # remove records inside the gaps
#     suboutldf = remove_gaps_from_obs(gaplist=gaplist, obsdf=outldf)

#     # restet to triple index
#     outldf = suboutldf.reset_index().set_index(["name", "datetime", "obstype"])

#     return outldf


# # =============================================================================
# # Helpers
# # =============================================================================
# def apply_debias_era5_gapfill(
#     gapslist,
#     dataset,
#     eraModelData,
#     debias_settings,
#     obstype="temp",
#     overwrite_fill=False,
# ):
#     """Fill all gaps using ERA5 debiaset modeldata.

#     Parameters
#     ----------
#     gapslist : list
#         list of all gaps.
#     dataset : metobs_toolkit.Dataset
#         Dataset to fill the gaps of.
#     eraModelData : metobs_toolkit.Modeldata
#         Modeldata to use for gapfilling.
#     debias_settings : dict
#         Debias settings.
#     obstype : str, optional
#         MetObs observationtype to fill gaps for. The default is "temp".
#     overwrite_fill : bool, optional
#         If True, the filled values are overwritten. The default is False.

#     Returns
#     -------
#     None.

#     """
#     gapfill_settings = dataset.settings.gap["gaps_fill_info"]

#     # Convert modeldata to the same timzone as the data
#     targettz = str(dataset.df.index.get_level_values("datetime").tz)
#     eraModelData._conv_to_timezone(targettz)

#     for gap in gapslist:
#         if (not overwrite_fill) & (not gap.gapfill_df.empty):
#             logger.warning(
#                 f"Gap {gap.name} is already filled with {gap.gapfill_technique} and will not be overwirtten. Set overwrite_fill to True to overwrite."
#             )
#             continue

#         logger.info(f" Era5 gapfill for {gap}")
#         gap.gapfill_technique = gapfill_settings["label"]["model_debias"]

#         # avoid passing full dataset around
#         station = dataset.get_station(gap.name)

#         # Update gap attributes
#         gap.update_gaps_indx_in_obs_space(
#             obsdf=station.df,
#             outliersdf=station.outliersdf,
#             dataset_res=station.metadf["dataset_resolution"].squeeze(),
#         )

#         # get leading and trailing period
#         leading_obs, trailing_obs = create_leading_trailing_debias_periods(
#             station=station,
#             gap=gap,
#             debias_period_settings=debias_settings["debias_period"],
#             obstype=obstype,
#         )

#         # check if leading/trailing is valid
#         if leading_obs.empty | trailing_obs.empty:
#             logger.info(
#                 "No suitable leading or trailing period found. Gapfill not possible"
#             )
#             gap.gapfill_errormessage[obstype] = (
#                 "gapfill not possible: no leading/trailing period"
#             )

#             default_return = pd.Series(
#                 index=gap.exp_gap_idx, name=obstype, dtype="object"
#             )
#             gap.gapfill_errormessage[obstype] = (
#                 "gapfill not possible: no leading/trailing period"
#             )

#             default_return.name = obstype
#             gapfill_df = default_return.to_frame()
#             gapfill_df[obstype + "_" + gapfill_settings["label_columnname"]] = (
#                 gapfill_settings["label"]["model_debias"]
#             )

#             # update the gaps attributes
#             gap.gapfill_df = gapfill_df

#             continue

#         # extract model values at leading and trailing period
#         leading_model = eraModelData.interpolate_modeldata(leading_obs.index)
#         trailing_model = eraModelData.interpolate_modeldata(trailing_obs.index)

#         # TODO check if there is modeldata for the leading and trailing + obs period
#         if (leading_model[obstype].isnull().any()) | (
#             trailing_model[obstype].isnull().any()
#         ):
#             logger.info(
#                 "No modeldata for the full leading/trailing period found. Gapfill not possible"
#             )
#             gap.gapfill_errormessage[obstype] = (
#                 "gapfill not possible: not enough modeldata"
#             )

#             default_return = pd.Series(
#                 index=gap.exp_gap_idx, name=obstype, dtype="object"
#             )
#             default_return.name = obstype
#             gapfill_df = default_return.to_frame()
#             gapfill_df[obstype + "_" + gapfill_settings["label_columnname"]] = (
#                 gapfill_settings["label"]["model_debias"]
#             )

#             # update the gaps attributes
#             gap.gapfill_df = gapfill_df
#             continue

#         # Get model data for gap timestamps
#         gap_model = eraModelData.interpolate_modeldata(gap.exp_gap_idx)

#         # apply bias correction
#         filled_gap_series, fill_info, err_message = make_era_bias_correction(
#             leading_model=leading_model,
#             trailing_model=trailing_model,
#             gap_model=gap_model,
#             leading_obs=leading_obs,
#             trailing_obs=trailing_obs,
#             obstype=obstype,
#         )

#         filled_gap_series.name = obstype
#         gapfill_df = filled_gap_series.to_frame()
#         gapfill_df[obstype + "_" + gapfill_settings["label_columnname"]] = (
#             gapfill_settings["label"]["model_debias"]
#         )

#         # update the gaps attributes
#         gap.gapfill_df = gapfill_df
#         gap.gapfill_technique = gapfill_settings["label"]["model_debias"]
#         gap.gapfill_info = fill_info
#         if bool(err_message):
#             gap.gapfill_errormessage = err_message


# def apply_interpolate_gaps(
#     gapslist,
#     obsdf,
#     outliersdf,
#     dataset_res,
#     gapfill_settings,
#     obstype="temp",
#     method="time",
#     max_consec_fill=100,
#     overwrite_fill=False,
# ):
#     """Fill all gaps with interpolation and update attributes.

#     Parameters
#     ----------
#     gapslist : list
#         list of all gaps.
#     obsdf : pandas.DataFrame
#         Dataframe with the observations.
#     outliersdf : pandas.DataFrame
#         Dataframe with the outliers (to find leading/trailing records).
#     dataset_res : pandas.Series
#         Frequency for all stations in a series.
#     gapfill_settings : dict
#         Gapfill settings.
#     obstype : str, optional
#         MetObs observationtype to fill gaps for. The default is "temp".
#     method : str, optional
#         Numpy interpolation method. The default is "time".
#     max_consec_fill : int, optional
#         Maximum number of consecutive records to fill. The default is 100.
#     overwrite_fill : bool, optional
#         If True, the filled values are overwritten. The default is False.

#     Returns
#     -------
#     None.

#     """
#     for gap in gapslist:
#         if (not overwrite_fill) & (not gap.gapfill_df.empty):
#             logger.warning(
#                 f"Gap {gap.name} is already filled with {gap.gapfill_technique} and will not be overwirtten. Set overwrite_fill to True to overwrite."
#             )
#             continue
#         gapfill_series = interpolate_gap(
#             gap=gap,
#             obsdf=xs_save(obsdf, gap.name, level="name", drop_level=False),
#             outliersdf=xs_save(outliersdf, gap.name, level="name", drop_level=False),
#             dataset_res=dataset_res.loc[gap.name],
#             obstype=obstype,
#             method=method,
#             max_consec_fill=max_consec_fill,
#         )

#         gapfill_series.name = obstype
#         gapfill_df = gapfill_series.to_frame()
#         gapfill_df[obstype + "_" + gapfill_settings["label_columnname"]] = (
#             gapfill_settings["label"]["linear"]
#         )

#         # update the gaps attributes
#         gap.gapfill_df = gapfill_df
#         gap.gapfill_technique = gapfill_settings["label"]["linear"]


# def make_gapfill_df(gapslist):
#     """Create a dataframe with all filled values of all gaps."""
#     if not bool(gapslist):
#         # no gaps (will be in automatic gapfill if method is not triggerd)
#         return pd.DataFrame()
#     concatlist = []
#     for gap in gapslist:
#         subgapfill = gap.gapfill_df.reset_index()
#         subgapfill["name"] = gap.name
#         subgapfill = subgapfill.set_index(["name", "datetime"])

#         concatlist.append(subgapfill)

#     filldf = concat_save(concatlist).sort_index()

#     # When gapfill could (paritally) not been fulfilled,
#     # their values (=Nan) must be removed from gapfill,
#     # so they will be plotted as gaps
#     filldf = filldf.dropna()

#     return filldf


def find_gaps(df):
    """Find missing timestamps and gaps in the observations.

    Looking for missing timestaps by assuming an observation frequency. The assumed frequency is the highest occuring frequency PER STATION.
    If missing observations are detected, they can be catogirized as a missing timestamp or as gap.

    A gap is define as a sequence of missing values with more than N repetitive missing values. N is define in the QC settings.



    Parameters
    ----------
    df : pandas.DataFrame
        The observations dataframe of the dataset object (Dataset.df)
    gapsize_n : int
        The minimum number of consecutive missing observations to identify the
        period as a gap.

    Returns
    -------
    gap_list : list
        A list of tuples with (name, startdt, enddt) of the gaps.

    """
    gap_list = []

    # missing timestamp per station (because some stations can have other frequencies!)
    stationnames = df.index.get_level_values(level="name").unique()
    for station in stationnames:
        # find missing timestamps
        timestamps = xs_save(df, station, level="name").index
        likely_freq = get_likely_frequency(timestamps, method="highest", simplify=False)
        assert likely_freq.seconds > 0, "The frequency is not positive!"

        missing_datetimeseries = (
            pd.date_range(
                start=timestamps.min(), end=timestamps.max(), freq=likely_freq
            )
            .difference(timestamps)
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
            gapdef = tuple(
                [station, datetime_of_gap_records.min(), datetime_of_gap_records.max()]
            )

            gap_list.append(gapdef)

    return gap_list
