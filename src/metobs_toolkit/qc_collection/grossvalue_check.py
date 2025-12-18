import logging
from typing import Union
import pandas as pd


from .whitelist import SensorWhiteSet
from metobs_toolkit.backend_collection.decorators import log_entry

logger = logging.getLogger("<metobs_toolkit>")


@log_entry
def gross_value_check(
    records: pd.Series,
    lower_threshold: Union[int, float],
    upper_threshold: Union[int, float],
    sensorwhiteset: SensorWhiteSet,
) -> pd.DatetimeIndex:
    """
    Identify outliers in a time series based on lower and upper thresholds.

    Parameters
    ----------
    records : pd.Series
        Series with a datetime-like index to check.
    lower_threshold : int or float
        Threshold below which records are flagged as outliers.
    upper_threshold : int or float
        Threshold above which records are flagged as outliers.
    sensorwhiteset : SensorWhiteSet
        A SensorWhiteSet instance containing timestamps that should be excluded from outlier detection.
        Records matching the whiteset criteria will not be flagged as outliers even if they fall
        outside the threshold range.

    Returns
    -------
    pd.DatetimeIndex
        Timestamps of outlier records.


    """

    # Drop NaN values
    records = records.dropna()
    # Identify outliers
    outliers_idx = records[
        (records < lower_threshold) | (records > upper_threshold)
    ].index

    # Exclude white records if provided
    outliers_idx = sensorwhiteset.catch_white_records(outliers_idx=outliers_idx)

    logger.debug("Exiting function gross_value_check.")
    return outliers_idx
