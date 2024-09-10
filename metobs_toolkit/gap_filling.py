#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 17:05:26 2023

@author: thoverga
"""

import numpy as np
import pandas as pd

import logging
from metobs_toolkit.df_helpers import (
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


def _create_anchor_df_for_leading_trailing_periods_by_size(
    Gap, Dataset, n_lead_records, n_trail_records, max_lead_duration, max_trail_duration
):
    """Helper for constructing an anchorsdf (used by interpolation)

    The anchordf is constructed per gap, and is done by locating non-nan N records (record space),
    adjacent to the gap. Where N is specified for the leading and trailing period.

    An extra filter to maximum lead/trail duration is applied to limit the
    duration (time-space) of the learning periods.

    This method is similar to _create_anchor_df_for_leading_trailing_periods(),
    and the record-space arguments are converted to time-space arguments so that
    this function can be maximally reused.


    Parameters
    ----------
    Gap : metobs_toolkit.Gap
        The gap to construct the anchorsdf for.
    Dataset : metobs_toolkit.Dataset
        The Dataset where the Gap is located in.
    n_lead_records : int
        The number of required leading records.
    n_trail_records : int
        The number of required trailing records.
    max_lead_duration : Timedelta or None
        The maximum timedifference between the leading and gap.startdt timestamps.
    max_trail_duration : TYPE
        The maximum timedifference between the trailing and gap.enddt timestamps.

    Returns
    -------
    anchor_df : pandas.Dataframe
        The anchorddf dataframe
    lead_msg : str
        Status message of the leading period.
    trail_msg : str
        Status message of the trailing period.

    """

    obsname = Gap.obstype.name
    sta_obs_series = xs_save(Dataset.df, Gap.name, "name", drop_level=True)
    sta_obs_series = xs_save(sta_obs_series, obsname, "obstype", drop_level=True)
    sta_obs_series = sta_obs_series["value"]

    # Reuse the _create_anchor_df_leading_trailing_periods method as much as possible,
    # therefore, convert the n_.. to timedelta

    # covert nlead_records to leading_period duration
    lead_period = (
        sta_obs_series[sta_obs_series.index < Gap.startdt].dropna().sort_index()
    )
    if max_lead_duration is not None:
        lead_period = lead_period[
            lead_period.index >= (Gap.startdt - max_lead_duration)
        ]

    lead_period_startdt = lead_period[-n_lead_records:].index.min()
    lead_duration = Gap.startdt - lead_period_startdt

    # covert ntrail_records to trailing_period duration
    trail_period = (
        sta_obs_series[sta_obs_series.index > Gap.enddt].dropna().sort_index()
    )
    if max_trail_duration is not None:
        trail_period = trail_period[
            trail_period.index <= (Gap.enddt + max_trail_duration)
        ]

    trail_period_enddt = trail_period[:n_trail_records].index.max()
    trail_duration = trail_period_enddt - Gap.enddt

    anchor_df = _create_anchor_df_for_leading_trailing_periods(
        Gap=Gap,
        Dataset=Dataset,
        leading_period_duration=lead_duration,
        trailing_period_duration=trail_duration,
    )

    # subset to max lead and trail durations
    if max_lead_duration is not None:
        anchor_df = anchor_df

    # create msg's
    if (
        anchor_df[anchor_df["fill_method"] == "leading period"].shape[0]
        < n_lead_records
    ):
        lead_msg = f"to few leading records ({anchor_df[anchor_df['fill_method'] == 'leading period'].shape[0]} found but {n_lead_records} needed)."
    else:
        lead_msg = "ok"

    if (
        anchor_df[anchor_df["fill_method"] == "trailing period"].shape[0]
        < n_trail_records
    ):
        trail_msg = f"to few trailing records ({anchor_df[anchor_df['fill_method'] == 'trailing period'].shape[0]} found but {n_lead_records} needed)."
    else:
        trail_msg = "ok"

    # write msg column in anchordf
    anchor_df.loc[anchor_df["fill_method"] == "leading period", "msg"] = lead_msg
    anchor_df.loc[anchor_df["fill_method"] == "trailing period", "msg"] = trail_msg
    return anchor_df, lead_msg, trail_msg


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
    (_, lead_period, lead_vals) = Gap._get_leading_period(
        observations_series=sta_obs_series,
        leading_period_duration=leading_period_duration,
    )

    (_, trail_period, trail_vals) = Gap._get_trailing_period(
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
    Helper method to combine the anchorsdf and the gapdf into one dataframe.
    This combined dataframe is then used by the fill method.

    Modeldata is extracted for all these records (and interpolated to the
    to match time resolution).

    Each record is labeled by 'leading period', 'trailing period' or 'gap'

    Parameters
    ----------
    Modeldata : metobs_toolkit.Modeldata
        The modeldata to use for the gap filling

    Returns
    -------
    filldf : pandas.DataFrame
        The dataframe structured as the Gap.gapdf, that combines the
        anchors and the gap records and the corresponding modelvalues.

    """
    # anchordf = Gap.anchordf
    # obsname = Gap.obstype.name

    debiasdf = Modeldata._interpolate_modeldata(anchordf.index)
    assert (
        obsname in debiasdf.columns
    ), f"{obsname} not present in the modeldata: {Modeldata}"
    debiasdf = debiasdf[[obsname]]
    debiasdf = debiasdf.rename(columns={obsname: "modelvalues"})
    debiasdf["obsvalues"] = anchordf[obsname]
    debiasdf["fill_method"] = anchordf["fill_method"]

    # add the gap period
    gapdf = Modeldata._interpolate_modeldata(gapdf.index)
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
    """Helper method to label (in the msg column) all records of the anchorsdf.

    This method is applied when filling uses diurnal timestamps, thus the
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

    return fillvalues, msg


# =============================================================================
# Docstring test
# =============================================================================
if __name__ == "__main__":
    from metobs_toolkit.doctest_fmt import setup_and_run_doctest

    setup_and_run_doctest()
