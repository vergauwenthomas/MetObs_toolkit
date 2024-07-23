#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 17:05:26 2023

@author: thoverga
"""

import numpy as np
import pandas as pd
from datetime import timedelta
import logging

from metobs_toolkit.df_helpers import (
    remove_outliers_from_obs,
    init_multiindexdf,
    format_outliersdf_to_doubleidx,
    xs_save,
)

logger = logging.getLogger(__name__)


# =============================================================================
# Helpers
# =============================================================================
def _add_diurnal_timestamps(df):
    if df.empty:
        # when the gap is at the beginning/end, there are no ankers --> df is empty
        df["hour"] = np.nan
        df["minutes"] = np.nan
        df["seconds"] = np.nan

    else:
        df["hour"] = df.index.hour
        df["minutes"] = df.index.minute
        df["seconds"] = df.index.second

    return df


def _create_anchor_df_for_leading_trailing_periods(
    Gap, Dataset, leading_period_duration, trailing_period_duration
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
    obsname = Gap.obstype.name
    sta_obs_series = xs_save(Dataset.df, Gap.name, "name", drop_level=True)
    sta_obs_series = xs_save(sta_obs_series, obsname, "obstype", drop_level=True)
    sta_obs_series = sta_obs_series["value"]

    # 1. Get leading and trailing info
    # get leading record, check validity and add to the gapfilldf
    (_, lead_period, lead_vals) = Gap.get_leading_period(
        observations_series=sta_obs_series,
        leading_period_duration=leading_period_duration,
    )

    (_, trail_period, trail_vals) = Gap.get_trailing_period(
        observations_series=sta_obs_series,
        trailing_period_duration=trailing_period_duration,
    )

    # 2. Create anchordf
    _leaddf_idx = pd.MultiIndex.from_arrays(
        arrays=[[Gap.name] * lead_period.shape[0], lead_period],
        names=["name", "datetime"],
    )

    _leaddf = pd.DataFrame(
        data={
            Gap.obstype.name: lead_vals,
            "fill_method": ["leading period"],
            # "msg": lead_msg,
        },
        index=_leaddf_idx,
    )

    _traildf_idx = pd.MultiIndex.from_arrays(
        arrays=[[Gap.name] * trail_period.shape[0], trail_period],
        names=["name", "datetime"],
    )

    _traildf = pd.DataFrame(
        data={
            Gap.obstype.name: trail_vals,
            "fill_method": ["trailing period"],
            # "msg": trail_msg,
        },
        index=_traildf_idx,
    )

    anchor_df = pd.concat([_leaddf, _traildf]).sort_index()

    # ensure a datetime index
    anchor_df = anchor_df.reset_index()
    anchor_df["datetime"] = pd.to_datetime(anchor_df["datetime"])
    anchor_df = anchor_df.set_index(["name", "datetime"])

    return anchor_df


def _combine_learning_and_gap_to_one_df(Modeldata, anchordf, gapdf, obsname):
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
    # anchordf = Gap.anchordf
    # obsname = Gap.obstype.name

    debiasdf = Modeldata.interpolate_modeldata(anchordf.index)
    assert (
        obsname in debiasdf.columns
    ), f"{obsname} not present in the modeldata: {Modeldata}"
    debiasdf = debiasdf[[obsname]]
    debiasdf = debiasdf.rename(columns={obsname: "modelvalues"})
    debiasdf["obsvalues"] = anchordf[obsname]
    debiasdf["fill_method"] = anchordf["fill_method"]

    # add the gap period
    gapdf = Modeldata.interpolate_modeldata(gapdf.index)
    gapdf = gapdf[[obsname]]
    gapdf = gapdf.rename(columns={obsname: "modelvalues"})
    gapdf["obsvalues"] = np.nan
    gapdf["fill_method"] = "gap"

    filldf = pd.concat([debiasdf, gapdf])
    filldf = filldf.sort_index()

    return filldf


def _label_anchors_for_diurnal_gapfilling(
    anchordf, gapdf, obsname, min_anchors_per_diurnal_timestamp
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

    anchordf["msg"] = "_init_"  # these will be overwritten
    # aggregate to diurnal groups and count number of samples For all anchors
    anchordf = anchordf.dropna()
    anchordf = anchordf.reset_index().set_index("datetime")
    anchordf = _add_diurnal_timestamps(anchordf)
    anchordf = anchordf.reset_index()  # else the datetimes will be lost when mergeing

    samplesizes = anchordf.groupby(["hour", "minutes", "seconds"])[obsname].count()
    samplesizes.name = "anchors"

    # aggregate the gaps records to diurnal groups, to find the diurnal timestamps that will be used
    gap_records = gapdf.reset_index().set_index("datetime")
    gap_records = _add_diurnal_timestamps(gap_records)
    gap_records = gap_records.groupby(["hour", "minutes", "seconds"])[obsname].count()

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

    # ensure a datetime index
    anchordf = anchordf.reset_index()
    anchordf["datetime"] = pd.to_datetime(anchordf["datetime"])
    anchordf = anchordf.set_index(["name", "datetime"])

    return anchordf


# =============================================================================
# Modeldata gapfillers
# =============================================================================


def _raw_modeldata_fill(filldf):
    """Return series of raw modeldata for gap records"""

    # Construct returns
    fillvalues = filldf["modelvalues"]
    msg = filldf.apply(
        lambda x: f"Modelvalue: {x['modelvalues']:.2f}, without correction.", axis=1
    )
    return fillvalues, msg


# dit nu overal toepassen (de msg return + opvangen in gap module)
def _debias_modeldata_fill(filldf):
    # calculate one bias for full anchor period
    anchors = filldf[filldf["fill_method"] != "gap"]
    bias = (anchors["modelvalues"] - anchors["obsvalues"]).mean()

    # correct the modelvalues with a single bias value
    gapfill = filldf[filldf["fill_method"] == "gap"]

    # Construct returns
    fillvalues = gapfill["modelvalues"] - bias
    msg = gapfill.apply(
        lambda x: f"Modelvalue: {x['modelvalues']:.2f} corrected with a {bias:.2f} bias.",
        axis=1,
    )

    return fillvalues, msg


def _diurnal_debias_modeldata_fill(filldf):
    # subset to anchors and gap
    anchordf = filldf[filldf["fill_method"] != "gap"]
    gapdf = filldf[filldf["fill_method"] == "gap"]
    anchordf = anchordf.reset_index().set_index("datetime")[
        ["modelvalues", "obsvalues"]
    ]
    gapdf = gapdf.reset_index().set_index("datetime")[["modelvalues", "obsvalues"]]

    # create diurnal record categories (hour, minutes, seconds)
    anchordf = _add_diurnal_timestamps(anchordf)
    gapdf = _add_diurnal_timestamps(gapdf)

    # calculate instantanious biases
    anchordf["diff"] = anchordf["modelvalues"] - anchordf["obsvalues"]
    # calculate biases for dirunal records
    diurnalbias = (
        anchordf.groupby(["hour", "minutes", "seconds"]).agg("mean").reset_index()
    )

    # merge biases to the gapdf
    gapdf = gapdf.reset_index().merge(
        diurnalbias[["diff", "hour", "minutes", "seconds"]],
        how="left",
        on=["hour", "minutes", "seconds"],
    )

    # revert the multiindex
    gapdf["name"] = filldf.index.get_level_values("name")[0]
    gapdf = gapdf.set_index(["name", "datetime"])

    # Construct returns
    fillvalues = gapdf["modelvalues"] - gapdf["diff"]  # correct with diurnal biasses

    def _msg_writer(row):
        if not np.isnan(row["diff"]):
            return f"Modelvalue: {row['modelvalues']:.2f} corrected with a {row['diff']:.2f} bias."
        else:
            return f"Modelvalue: {row['modelvalues']:.2f} cannot be corrected because anchor samples does not meet requirements"

    msg = gapdf.apply(
        lambda x: _msg_writer(x),
        axis=1,
    )

    return fillvalues, msg


def _weighted_diurnal_debias_modeldata(filldf):

    # subset to leading, trailing and gap
    leadingdf = filldf[filldf["fill_method"] == "leading period"]
    trailingdf = filldf[filldf["fill_method"] == "trailing period"]
    gapdf = filldf[filldf["fill_method"] == "gap"]

    leadingdf = leadingdf.reset_index().set_index("datetime")[
        ["modelvalues", "obsvalues"]
    ]
    trailingdf = trailingdf.reset_index().set_index("datetime")[
        ["modelvalues", "obsvalues"]
    ]
    gapdf = gapdf.reset_index().set_index("datetime")[["modelvalues", "obsvalues"]]

    # create diurnal record categories (hour, minutes, seconds)
    leadingdf = _add_diurnal_timestamps(leadingdf)
    trailingdf = _add_diurnal_timestamps(trailingdf)
    gapdf = _add_diurnal_timestamps(gapdf)

    # calculate instantanious biases
    leadingdf["lead_diff"] = leadingdf["modelvalues"] - leadingdf["obsvalues"]

    # calculate biases for dirunal records
    leadingbias = (
        leadingdf.groupby(["hour", "minutes", "seconds"]).agg("mean").reset_index()
    )

    trailingdf["trail_diff"] = trailingdf["modelvalues"] - trailingdf["obsvalues"]
    # calculate biases for dirunal records
    trailingbias = (
        trailingdf.groupby(["hour", "minutes", "seconds"]).agg("mean").reset_index()
    )

    # merge biases to the gapdf
    gapdf = gapdf.reset_index()
    gapdf = gapdf.merge(
        leadingbias[["lead_diff", "hour", "minutes", "seconds"]],
        how="left",
        on=["hour", "minutes", "seconds"],
    )

    gapdf = gapdf.merge(
        trailingbias[["trail_diff", "hour", "minutes", "seconds"]],
        how="left",
        on=["hour", "minutes", "seconds"],
    )

    # compute weights for each gap record
    gapstart = gapdf["datetime"].min()
    gapend = gapdf["datetime"].max()
    gapduration = gapend - gapstart
    gapdf["lead_weight"] = 1.0 - ((gapdf["datetime"] - gapstart) / gapduration)
    gapdf["trail_weight"] = 1.0 - gapdf["lead_weight"]

    # compute a finall diff (correction factor)
    gapdf["diff"] = (gapdf["lead_diff"] * gapdf["lead_weight"]) + (
        gapdf["trail_diff"] * gapdf["trail_weight"]
    )

    # revert the multiindex
    gapdf["name"] = filldf.index.get_level_values("name")[0]
    gapdf = gapdf.set_index(["name", "datetime"])

    # Construct returns
    fillvalues = gapdf["modelvalues"] - gapdf["diff"]  # correct with diurnal biasses

    def _msg_writer(row):
        if not np.isnan(row["diff"]):
            return f"Modelvalue: {row['modelvalues']:.2f} corrected with ({row['lead_diff']:.2f} x {row['lead_weight']:.2f} + {row['trail_diff']:.2f} x {row['trail_weight']:.2f}) bias."
        else:
            return f"Modelvalue: {row['modelvalues']:.2f} cannot be corrected because anchor samples does not meet requirements"

    msg = gapdf.apply(
        lambda x: _msg_writer(x),
        axis=1,
    )

    # msg = gapdf.apply(
    #     lambda x: f"Modelvalue: {x['modelvalues']:.2f} corrected with ({x['lead_diff']:.2f} x {x['lead_weight']:.2f} + {x['trail_diff']:.2f} x {x['trail_weight']:.2f}) bias.",
    #     axis=1,
    # )
    # # If the trailing/leading periods does not contain all diurnal categories
    # # as the the gap, the corresponding diff cannot be merge and Nan's are added.
    # # We specify this in the message
    # missing_lead_cat = gapdf[~gapdf["lead_diff"].notnull()].index  # find the indices
    # msg.loc[missing_lead_cat] = (
    #     "This diurnal record is not represented in the leading period."
    # )
    # missing_trail_cat = gapdf[~gapdf["trail_diff"].notnull()].index  # find the indices
    # msg.loc[missing_trail_cat] = (
    #     "This diurnal record is not represented in the trailing period."
    # )

    return fillvalues, msg


# # =============================================================================
# # Debiasing period
# # =============================================================================


# def get_sample_size(sample_duration_hours, sta):
#     """Get the number of records for a sample duration."""
#     stares = sta.metadf["dataset_resolution"].squeeze()
#     sample_size = timedelta(hours=sample_duration_hours) / stares
#     return int(sample_size)


# def create_leading_trailing_debias_periods(
#     station, gap, debias_period_settings, obstype
# ):
#     """Get the leading and trailing periods of a gap."""
#     # Get samplesizes
#     debias_pref_sample_size_leading = get_sample_size(
#         debias_period_settings["prefered_leading_sample_duration_hours"], station
#     )
#     debias_pref_sample_size_trailing = get_sample_size(
#         debias_period_settings["prefered_trailing_sample_duration_hours"], station
#     )
#     debias_min_sample_size_leading = get_sample_size(
#         debias_period_settings["minimum_leading_sample_duration_hours"], station
#     )
#     debias_min_sample_size_trailing = get_sample_size(
#         debias_period_settings["minimum_trailing_sample_duration_hours"], station
#     )

#     # get all observations that can be used for debias training
#     obs = station.df

#     # remove blacklist
#     # TODO
#     obs = remove_outliers_from_obs(
#         obs, format_outliersdf_to_doubleidx(station.outliersdf)
#     )

#     # add whitelist
#     # TODO

#     # only datetimes are relevant
#     obs = obs.reset_index()
#     obs = obs[["name", "datetime", obstype]]

#     # Select all leading and all trailing obs
#     leading_period = obs[obs["datetime"] < gap.startgap]
#     trailing_period = obs[obs["datetime"] > gap.endgap]
#     logger.debug(
#         f"   {leading_period.shape[0]} leading records, {trailing_period.shape[0]} trailing records."
#     )

#     # some derived integers
#     poss_shrinkage_leading = leading_period.shape[0] - debias_min_sample_size_leading
#     poss_shrinkage_trailing = trailing_period.shape[0] - debias_min_sample_size_trailing
#     poss_extention_leading = leading_period.shape[0] - debias_pref_sample_size_leading
#     poss_extention_trailing = (
#         trailing_period.shape[0] - debias_pref_sample_size_trailing
#     )

#     # check if desired sample sizes for leading and trailing are possible
#     if (leading_period.shape[0] >= debias_pref_sample_size_leading) & (
#         trailing_period.shape[0] >= debias_pref_sample_size_trailing
#     ):
#         logger.debug("leading and trailing periods are both available for debiassing.")
#         # both periods are oke
#         leading_df = leading_period[-debias_pref_sample_size_leading:]
#         trailing_df = trailing_period[:debias_pref_sample_size_trailing]

#     elif (leading_period.shape[0] <= debias_pref_sample_size_leading) & (
#         trailing_period.shape[0] >= debias_pref_sample_size_trailing
#     ):
#         logger.debug(
#             "leading periods for debiassing does not have a preferable size. Try translation/shrinkage ..."
#         )

#         # leading period to small, trailing period is OK

#         missing_records = debias_pref_sample_size_leading - leading_period.shape[0]

#         # 1 if the leading period is smaller thatn the minimum leading size --> return default
#         if poss_shrinkage_leading < 0:
#             leading_df = init_multiindexdf()
#             trailing_df = init_multiindexdf()  # this might be to strict
#             logger.debug(
#                 "The available leading debias samplesize is smaller than the minimum. A translation/shrinking is not possible."
#             )

#         # 2 Try translation without shrinkage

#         elif missing_records <= poss_extention_trailing:
#             # translation without shrinkage is possible
#             translation_trailing = missing_records

#             leading_df = leading_period
#             trailing_df = trailing_period[
#                 0 : (debias_pref_sample_size_trailing + translation_trailing)
#             ]

#             logger.debug(
#                 f"A translation of {translation_trailing} records is done towards the trailing period. (n_leading + n_trailing is conserved: {leading_df.shape[0] + trailing_df.shape[0]}"
#             )

#         # 3. Try if a translation is within the limits of shrinkage
#         elif (missing_records - poss_extention_trailing) <= poss_shrinkage_leading:
#             translation_trailing = poss_extention_trailing

#             leading_df = leading_period
#             trailing_df = trailing_period[
#                 0 : debias_pref_sample_size_trailing + translation_trailing
#             ]
#             logger.debug(
#                 f"A translation of {translation_trailing} records is done towards the trailing period. Since there was not engough translation space for the trailing obs, the condition n_leading + n_trailing is NOT conserved: {leading_df.shape[0] + trailing_df.shape[0]}. \
#                   Both leading and trailing sizes still achieves minimal size restrictions."
#             )
#         # 4. If all else fails, it is not possible to make a leading period
#         else:
#             logger.info(
#                 "The available leading samplesize can not reach minimal size restrictions."
#             )
#             # no translation is possible, even with shrinking
#             leading_df = init_multiindexdf()
#             trailing_df = init_multiindexdf()  # this might be to strict

#     elif (leading_period.shape[0] >= debias_pref_sample_size_leading) & (
#         trailing_period.shape[0] <= debias_pref_sample_size_trailing
#     ):
#         # leading period is ok, trailing period is to short
#         logger.debug(
#             "trailing periods for debiassing does not have a preferable size. Try translation/shrinkage ..."
#         )
#         missing_records = debias_pref_sample_size_trailing - trailing_period.shape[0]

#         # 1 if the trailing period is smaller thatn the minimum trailing size --> return default
#         if poss_shrinkage_trailing < 0:
#             leading_df = init_multiindexdf()  # might be to strict
#             trailing_df = init_multiindexdf()
#             logger.debug(
#                 "The available trailing debias samplesize is smaller than the minimum. A translation/shrinking is not possible."
#             )
#             # return

#         # 2 Try translation without shrinkage
#         elif missing_records <= poss_extention_leading:
#             # translation without shrinkage is possible
#             translation_leading = missing_records

#             leading_df = leading_period[
#                 -(debias_pref_sample_size_leading + translation_leading) :
#             ]
#             trailing_df = trailing_period
#             logger.debug(
#                 f"A translation of {translation_leading} records is done towards the leading period. (n_leading + n_trailing is conserved: {leading_df.shape[0] + trailing_df.shape[0]}"
#             )

#         # 3. Try if a translation is within the limits of shrinkage
#         elif (missing_records - poss_extention_leading) <= poss_shrinkage_trailing:
#             translation_leading = poss_extention_leading

#             leading_df = leading_period[
#                 -(debias_pref_sample_size_leading + translation_leading)
#             ]
#             trailing_df = trailing_period
#             logger.debug(
#                 f"A translation of {translation_leading} records is done towards the leading period. Since there was not engough translation space for the leading obs, the condition n_leading + n_trailing is NOT conserved: {leading_df.shape[0] + trailing_df.shape[0]}. \
#                   Both leading and trailing sizes still achieves minimal size restrictions."
#             )
#         # 4. If all else fails, it is not possible to make a trailing period
#         else:
#             # no translation is possible, even with shrinking
#             logger.info(
#                 "The available trailing samplesize can not reach minimal size restrictions."
#             )
#             leading_df = init_multiindexdf()  # this might be to strict
#             trailing_df = init_multiindexdf()

#     else:
#         # Both leading and trailing periods are not to small

#         # 1 does both (leading and trailing) still acchieves the minimal size condition for shrinking?
#         if (poss_shrinkage_leading >= 0) & (poss_shrinkage_trailing >= 0):
#             logger.debug(
#                 "Both leading and trailing periods do not have a prefered size, but still meet the minimal conditions."
#             )
#             leading_df = leading_period
#             trailing_df = trailing_period

#         else:
#             logger.info(
#                 "Both leading and trailing periods do not have a prefered size, and eighter of them does NOT meet minimal condition."
#             )
#             # either one of the periods does not reach minimal condition, so return default
#             leading_df = init_multiindexdf()
#             trailing_df = init_multiindexdf()

#     # convert to multiindex
#     if not leading_df.empty:
#         leading_df = leading_df.set_index(["name", "datetime"])
#     if not trailing_df.empty:
#         trailing_df = trailing_df.set_index(["name", "datetime"])

#     return leading_df, trailing_df


# def get_time_specific_biases(model, obs, obstype, period):
#     """Get hourly biases."""
#     diff = model - obs
#     diff = diff.reset_index().set_index("datetime")
#     diff["hours"] = diff.index.hour
#     diff["minutes"] = diff.index.minute
#     diff["seconds"] = diff.index.second

#     biases = diff.groupby(["name", "hours", "minutes", "seconds"])[obstype].mean()
#     biases.name = obstype + "_bias_" + period

#     biases = biases.reset_index()
#     return biases


# def make_era_bias_correction(
#     leading_model, trailing_model, gap_model, leading_obs, trailing_obs, obstype
# ):
#     """Make debias correction of the modeldata for a gap."""
#     error_message = ""
#     # 1. get lead timestamp biases
#     lead_biases = get_time_specific_biases(
#         model=leading_model, obs=leading_obs, obstype=obstype, period="lead"
#     )

#     # 2. get trailing timestamp biases
#     trail_biases = get_time_specific_biases(
#         model=trailing_model, obs=trailing_obs, obstype=obstype, period="trail"
#     )

#     # 3. apply bias correction on modeldata in gap

#     # linear interpolation of bias along the gap method:
#     gap_model["trail_weight"] = np.linspace(0.0, 1.0, gap_model.shape[0])
#     gap_model["lead_weight"] = 1.0 - gap_model["trail_weight"]

#     # aggregate to timestamps
#     gap_model["hours"] = gap_model.index.get_level_values("datetime").hour
#     gap_model["minutes"] = gap_model.index.get_level_values("datetime").minute
#     gap_model["seconds"] = gap_model.index.get_level_values("datetime").second

#     gap_model = gap_model.reset_index()

#     gap_model = gap_model.merge(
#         right=lead_biases[["hours", "minutes", "seconds", obstype + "_bias_lead"]],
#         how="left",
#         on=["hours", "minutes", "seconds"],
#     )

#     gap_model = gap_model.merge(
#         right=trail_biases[["hours", "minutes", "seconds", obstype + "_bias_trail"]],
#         how="left",
#         on=["hours", "minutes", "seconds"],
#     )

#     gap_model = gap_model.set_index(["name", "datetime"])

#     # Idea: if BOTH leadin and trailing (hourly) biases is available, than use
#     # use the debias corection (even if it is for a part of the gap!).
#     # If either one or both are missing, than no bias correction is applied
#     no_debias = gap_model[
#         (gap_model[obstype + "_bias_lead"].isnull())
#         | (gap_model[obstype + "_bias_trail"].isnull())
#     ].index
#     if not no_debias.empty:
#         error_message = f"No debias possible for these gap records: {no_debias},the gap will be filled by model data without bias correction. "
#         logger.warning(error_message)

#     # set weights to zero if not debias correction can be applied on that record
#     gap_model.loc[no_debias, obstype + "_bias_trail"] = 0.0
#     gap_model.loc[no_debias, obstype + "_bias_lead"] = 0.0

#     # 5. compute the debiased fill value
#     # leave this dataframe for debugging
#     gap_model[obstype + "_debiased_value"] = gap_model[obstype] - (
#         (gap_model["lead_weight"] * gap_model[obstype + "_bias_lead"])
#         + (gap_model["trail_weight"] * gap_model[obstype + "_bias_trail"])
#     )

#     # 7. format gapmodel
#     gap_model["time"] = (
#         gap_model["hours"].astype(str).str.zfill(2)
#         + ":"
#         + gap_model["minutes"].astype(str).str.zfill(2)
#         + ":"
#         + gap_model["seconds"].astype(str).str.zfill(2)
#     )
#     gap_model = gap_model.rename(columns={obstype: f"{obstype}_model_value"})

#     # 6. make returen
#     returnseries = gap_model[obstype + "_debiased_value"]
#     returnseries.name = obstype
#     return returnseries, gap_model, error_message
