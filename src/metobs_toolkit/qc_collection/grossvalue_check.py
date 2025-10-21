import logging
from typing import Union
import pandas as pd

from .common_functions import catch_white_records
from metobs_toolkit.backend_collection.loggingmodule import log_entry

logger = logging.getLogger("<metobs_toolkit>")


@log_entry
def gross_value_check(
    records: pd.Series,
    lower_threshold: Union[int, float],
    upper_threshold: Union[int, float],
    white_records: Union[pd.DatetimeIndex, None] = None,
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
    white_records : pd.DatetimeIndex, optional
        A DatetimeIndex containing timestamps that should be excluded from outlier detection.
        These "white records" are known valid values and will not be flagged as outliers
        even if they fall outside the threshold range. The default is None.

    Returns
    -------
    pd.DatetimeIndex
        Timestamps of outlier records.


    """

    # Drop NaN values
    records = records.dropna()
    # Identify outliers
    outliers_idx = records[(records < lower_threshold) | (records > upper_threshold)].index
    
    # Catch the white records
    if white_records is not None:
        outliers_idx = catch_white_records(outliers_idx, white_records)
    
    logger.debug("Exiting function gross_value_check.")
    return outliers_idx
