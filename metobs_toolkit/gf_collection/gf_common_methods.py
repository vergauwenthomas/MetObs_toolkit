import logging
from typing import Union
import pandas as pd
import numpy as np

logger = logging.getLogger("<metobs_toolkit>")


def create_a_combined_df(leadseries, trailseries, gap):

    leaddf = pd.DataFrame(
        index=leadseries.index, data={"value": leadseries, "label": "lead"}
    )

    traildf = pd.DataFrame(
        index=trailseries.index, data={"value": trailseries, "label": "trail"}
    )
    gapdf = pd.DataFrame(
        index=gap.records.index, data={"value": gap.records, "label": "gap"}
    )
    return pd.concat([leaddf, traildf, gapdf]).sort_index()


def add_modeldata_to_combdf(combineddf, modeltimeseries):
    modelseries = modeltimeseries.series

    # 1. Ensure both series have the same timezone
    if modelseries.index.tz != combineddf.index.tz:
        modelseries = modelseries.tz_convert(combineddf.index.tz)

    modeldf = pd.DataFrame(index=modelseries.index, data={"modelvalue": modelseries})
    # 2. Reindex modelseries to match combineddf, interpolating if necessary
    combdf_reindexed = (
        pd.concat([combineddf, modeldf])
        .sort_index()
        .interpolate(method="time", limit_area="inside")
    )
    combdf_reindexed = combdf_reindexed.loc[
        combineddf.index
    ]  # subset to gap + leading + trailing period

    # A duplicated row at round-hours is introduced, but the label is nan -> remove these rows
    combdf_reindexed = combdf_reindexed.dropna(subset=["label"])

    # The iterpolate call did also interpolate the 'value' column, revert this
    combdf_reindexed.loc[combdf_reindexed["label"] == "gap", "value"] = np.nan

    # dulicates are introduced when timestamps are both in modelseries and gapseries
    combdf_reindexed = combdf_reindexed[
        ~combdf_reindexed.index.duplicated(keep="first")
    ]

    return combdf_reindexed


def check_if_modeltimeseries_is_compatible(
    gap, modeltimeseries, lp_duration: pd.Timedelta, tp_duration: pd.Timedelta
):

    # check if start of modeldata is before gapstart
    if modeltimeseries.start_datetime <= (gap.start_datetime - lp_duration):
        pass
    else:
        return (
            False,
            f"start of modeltimeseries is not compatible with the start of the gap (minus the leading period size): {modeltimeseries.start_datetime} <= ({gap.start_datetime} - {lp_duration}) == False.",
        )

    # check if end of modeldata is after gapend
    if modeltimeseries.end_datetime >= (gap.end_datetime + tp_duration):
        pass
    else:
        return (
            False,
            f"end of modeltimeseries is not compatible with the end of the gap (plus the trailing period size): {modeltimeseries.end_datetime} >= ({gap.end_datetime} + {tp_duration}) == False.",
        )

    # Check if the model represents the same obstype as the gap
    if modeltimeseries.obstype.is_compatible_with(gap.obstype):
        pass
    else:
        return (
            False,
            f"The obstypes of the modeltimeseries is not compatible to that of the gap: {modeltimeseries.obstype} == {gap.obstype} == False.",
        )

    return True, "_"


def get_trailing_period(
    gap,
    sensordata,
    n_records: int,
    duration: Union[pd.Timedelta, None] = None,
    fixed_by_records=True,
    fixed_by_duration=False,
):

    # Compute the leading period from non-nan records
    sta_obs_series = sensordata.series

    potential_trail_period = (
        sta_obs_series[sta_obs_series.index > gap.end_datetime].dropna().sort_index()
    )

    if (fixed_by_duration, fixed_by_records) == (True, False):

        # fixed_by_duration case
        tp = potential_trail_period[
            potential_trail_period.index <= (gap.end_datetime + duration)
        ]
    elif (fixed_by_duration, fixed_by_records) == (False, True):
        # fixed by N records case
        tp = potential_trail_period[:n_records]

        # Extra filter based on duration
        if duration is not None:
            tp = tp[tp.index <= (gap.end_datetime + duration)]

    else:
        # invalid combination
        raise NotImplementedError(
            f"The combination fixed_by_records: {fixed_by_records} and fixed_by_duration:{fixed_by_duration} is not implemented."
        )

    # Check validity
    if tp.shape[0] < n_records:
        logger.debug(f"to few trailing records: {tp}, but needed {n_records}")
        msg = f"to few trailing records ({tp.shape[0]} found but {n_records} needed"
        continueflag = False
    else:
        msg = "ok"
        continueflag = True
    return tp, continueflag, msg


def get_leading_period(
    gap,
    sensordata,
    n_records: int,
    duration: Union[pd.Timedelta, None] = None,
    fixed_by_records=True,
    fixed_by_duration=False,
):

    # Compute the leading period from non-nan records
    sta_obs_series = sensordata.series

    potential_lead_period = (
        sta_obs_series[sta_obs_series.index < gap.start_datetime].dropna().sort_index()
    )

    if (fixed_by_duration, fixed_by_records) == (True, False):
        # fixed_by_duration case
        lp = potential_lead_period[
            potential_lead_period.index >= (gap.start_datetime - duration)
        ]
    elif (fixed_by_duration, fixed_by_records) == (False, True):
        # fixed by N records case
        lp = potential_lead_period[-n_records:]

        # Extra filter based on duration
        if duration is not None:
            lp = lp[lp.index >= (gap.start_datetime - duration)]

    else:
        # invalid combination
        raise NotImplementedError(
            f"The combination fixed_by_records: {fixed_by_records} and fixed_by_duration:{fixed_by_duration} is not implemented."
        )

    # lp should have a minimum samplesize
    if lp.shape[0] < n_records:
        if duration is not None:
            msg = f"to few leading records ({lp.shape[0]} (with a maximum distance of {duration} wrt {gap.start_datetime}) found but {n_records} needed"
        else:
            msg = f"to few leading records ({lp.shape[0]} found but {n_records} needed"
        logger.debug(msg)
        continueflag = False
    else:
        msg = "ok"
        continueflag = True

    return lp, continueflag, msg
