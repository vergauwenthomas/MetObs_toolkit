import logging
import pandas as pd
from typing import Union

logger = logging.getLogger(__name__)

# Set up logging
logger = logging.getLogger(__file__)
logger.setLevel(logging.INFO)


def gross_value_check(
    records: pd.Series, lower_threshold: int | float, upper_threshold: int | float
) -> pd.DatetimeIndex:
    """Identify outliers in a timeseries based on thresholds.

    Parameters
    ----------
    records : pd.Series
        records with datetime-like index to check.
    lower_threshold: int or float
        Thresholds to flag records below as outliers.
    upper_threshold : int or float
        Thresholds to flag records above as outliers.

    Returns
    -------
    pd.DatetimeIndex
        Timestamps of outlier records.

    """

    logger.info("Entering function gross_value_check.")

    # Validate argument types
    if not isinstance(records, pd.Series):
        raise TypeError("Argument 'records' must be of type pandas.Series.")
    if not isinstance(lower_threshold, (float, int)):  # Allow int for thresholds
        raise TypeError("Argument 'lower_threshold' must be of type float or int.")
    if not isinstance(upper_threshold, (float, int)):  # Allow int for thresholds
        raise TypeError("Argument 'upper_threshold' must be of type float or int.")

    # Drop NaN values
    records = records.dropna()
    # Identify outliers
    outliers = records[(records < lower_threshold) | (records > upper_threshold)].index
    logger.info("Exiting function gross_value_check.")
    return outliers
