import logging
from typing import Union
import pandas as pd

logger = logging.getLogger("<metobs_toolkit>")


def gross_value_check(
    records: pd.Series,
    lower_threshold: Union[int, float],
    upper_threshold: Union[int, float],
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

    Returns
    -------
    pd.DatetimeIndex
        Timestamps of outlier records.


    """
    logger.debug("Entering function gross_value_check.")

    # Drop NaN values
    records = records.dropna()
    # Identify outliers
    outliers = records[(records < lower_threshold) | (records > upper_threshold)].index
    logger.debug("Exiting function gross_value_check.")
    return outliers
