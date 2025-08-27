import logging
from typing import Union
import pandas as pd

from metobs_toolkit.backend_collection.loggingmodule import log_entry

logger = logging.getLogger("<metobs_toolkit>")


@log_entry
def create_a_combined_df(
    leadseries: pd.Series, trailseries: pd.Series, gap
) -> pd.DataFrame:
    """
    Create a combined DataFrame from leading, trailing, and gap series.

    Parameters
    ----------
    leadseries : pd.Series
        The leading period time series.
    trailseries : pd.Series
        The trailing period time series.
    gap : object
        An object with a 'records' attribute containing the gap period as
        a pd.Series.

    Returns
    -------
    pd.DataFrame
        Combined DataFrame with columns 'value' and 'label', indexed by
        datetime.
    """
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


@log_entry
def add_modeldata_to_combdf(combineddf: pd.DataFrame, modeltimeseries) -> pd.DataFrame:
    """
    Add model data to a combined DataFrame, interpolating model values to
    match the index.

    Parameters
    ----------
    combineddf : pd.DataFrame
        DataFrame containing 'value' and 'label' columns for lead, trail,
        and gap periods.
    modeltimeseries : object
        An object with a 'series' attribute (pd.Series) representing
        model data.

    Returns
    -------
    pd.DataFrame
        DataFrame with an additional 'modelvalue' column, aligned and
        interpolated to the combineddf index.
    """
    modelseries = modeltimeseries.series

    # 1. Ensure both series have the same timezone
    if modelseries.index.tz != combineddf.index.tz:
        modelseries = modelseries.tz_convert(combineddf.index.tz)

    modeldf = pd.DataFrame(index=modelseries.index, data={"modelvalue": modelseries})
    # Interpolate modeldf to target timestamps

    # Create a new index that is the union of combineddf and modeldf
    # indices without duplicates
    all_timestamps_index = combineddf.index.union(modeldf.index).drop_duplicates()
    # Note: The reindexing and interpolation of modeldata is needed
    # because the model data is often only defined at a resolution of 1h,
    # while the combineddf is not always at 1h resolution.

    # Combine all timestamps of model and combineddf, because if
    # only using the combineddf index, and limit_area="inside", no modeldata
    # could be set for the extreme leading and trailing timestamps
    # if they appear not on the hour. Thus combine indexes, and after
    # the data is combined, subset to the combineddf index.

    modeldf = (
        modeldf.reindex(
            all_timestamps_index
        )  # set the index for observation frequencies
        .sort_index()  # sort the index
        .interpolate(method="time", limit_area="inside")  # interp sub hourly
    )
    combdf_reindexed = pd.concat([combineddf, modeldf], axis=1)

    combdf_reindexed = combdf_reindexed.loc[
        combineddf.index
    ]  # subset to gap + leading + trailing period

    return combdf_reindexed


@log_entry
def check_if_modeltimeseries_is_compatible(
    gap,
    modeltimeseries,
    lp_duration: pd.Timedelta,
    tp_duration: pd.Timedelta,
) -> tuple[bool, str]:
    """
    Check if a model time series is compatible with a gap and
    its leading/trailing periods.

    Parameters
    ----------
    gap : object
        An object with 'start_datetime', 'end_datetime',
        and 'obstype' attributes.
    modeltimeseries : object
        An object with 'start_datetime', 'end_datetime',
        and 'obstype' attributes.
    lp_duration : pd.Timedelta
        Duration of the leading period.
    tp_duration : pd.Timedelta
        Duration of the trailing period.

    Returns
    -------
    tuple of (bool, str)
        (True, "_") if compatible, otherwise (False, reason).
    """
    # Check if start of model data is before gap start
    if modeltimeseries.start_datetime <= (gap.start_datetime - lp_duration):
        pass
    else:
        return (
            False,
            f"Start of modeltimeseries is not compatible with the start of \
the gap (minus the leading period size): {modeltimeseries.start_datetime} <= \
({gap.start_datetime} - {lp_duration}) == False.",
        )

    # Check if end of model data is after gap end
    if modeltimeseries.end_datetime >= (gap.end_datetime + tp_duration):
        pass
    else:
        return (
            False,
            f"End of modeltimeseries is not compatible with the end of the \
gap (plus the trailing period size): {modeltimeseries.end_datetime} >= \
({gap.end_datetime} + {tp_duration}) == False.",
        )

    # Check if the model represents the same obstype as the gap
    if modeltimeseries.modelobstype.is_compatible_with(gap.obstype):
        pass
    else:
        return (
            False,
            f"The obstypes of the modeltimeseries is not compatible with that \
of the gap: {modeltimeseries.modelobstype} == {gap.obstype} == False.",
        )

    return True, "_"


@log_entry
def get_trailing_period(
    gap,
    sensordata,
    n_records: int,
    duration: Union[pd.Timedelta, None] = None,
    fixed_by_records: bool = True,
    fixed_by_duration: bool = False,
) -> tuple[pd.Series, bool, str]:
    """
    Get the trailing period after a gap from sensor data.

    Parameters
    ----------
    gap : object
        An object with 'end_datetime' attribute.
    sensordata : object
        An object with a 'series' attribute (pd.Series) containing sensor data.
    n_records : int
        Minimum number of records required in the trailing period.
    duration : pd.Timedelta or None, optional
        Maximum duration for the trailing period.
    fixed_by_records : bool, default True
        Whether to fix the trailing period by number of records.
    fixed_by_duration : bool, default False
        Whether to fix the trailing period by duration.

    Returns
    -------
    tuple
        (trailing_period: pd.Series, continueflag: bool, msg: str)
    """
    # Compute the trailing period from non-NaN records
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
            f"The combination fixed_by_records: {fixed_by_records} and \
fixed_by_duration: {fixed_by_duration} is not implemented."
        )

    # Check validity
    if tp.shape[0] < n_records:
        logger.debug(f"Too few trailing records: {tp}, but needed {n_records}")
        msg = f"Too few trailing records ({tp.shape[0]} found \
but {n_records} needed)"
        continueflag = False
    else:
        msg = "ok"
        continueflag = True
    return tp, continueflag, msg


@log_entry
def get_leading_period(
    gap,
    sensordata,
    n_records: int,
    duration: Union[pd.Timedelta, None] = None,
    fixed_by_records: bool = True,
    fixed_by_duration: bool = False,
) -> tuple[pd.Series, bool, str]:
    """
    Get the leading period before a gap from sensor data.

    Parameters
    ----------
    gap : object
        An object with 'start_datetime' attribute.
    sensordata : object
        An object with a 'series' attribute (pd.Series) containing sensor data.
    n_records : int
        Minimum number of records required in the leading period.
    duration : pd.Timedelta or None, optional
        Maximum duration for the leading period.
    fixed_by_records : bool, default True
        Whether to fix the leading period by number of records.
    fixed_by_duration : bool, default False
        Whether to fix the leading period by duration.

    Returns
    -------
    tuple
        (leading_period: pd.Series, continueflag: bool, msg: str)
    """
    # Compute the leading period from non-NaN records
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
            f"The combination fixed_by_records: {fixed_by_records} and \
fixed_by_duration: {fixed_by_duration} is not implemented."
        )

    # lp should have a minimum sample size
    if lp.shape[0] < n_records:
        if duration is not None:
            msg = f"Too few leading records ({lp.shape[0]} (with a \
maximum distance of {duration} wrt {gap.start_datetime}) found \
but {n_records} needed)"
        else:
            msg = f"Too few leading records ({lp.shape[0]} found \
but {n_records} needed)"
        logger.debug(msg)
        continueflag = False
    else:
        msg = "ok"
        continueflag = True

    return lp, continueflag, msg
