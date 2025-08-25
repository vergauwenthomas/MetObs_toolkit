import logging
from typing import Union
import pandas as pd

from metobs_toolkit.backend_collection.loggingmodule import log_entry

logger = logging.getLogger("<metobs_toolkit>")


@log_entry
def step_check(
    records: pd.Series,
    max_increase_per_second: Union[int, float],
    max_decrease_per_second: Union[int, float],
) -> pd.DatetimeIndex:
    """
    Check for 'spikes' and 'dips' in a time series.

    Tests if observations produce spikes in the time series. The maximum
    allowed increase and decrease per second is set in the arguments,
    and is tested for each record (with respect to the previous record).

    If the difference between two consecutive records (i.e., the spike or dip) is larger than the
    threshold, the record is flagged as an outlier.

    Parameters
    ----------
    records : pd.Series
        A pandas Series containing the time series data to be checked. The index should be datetime-like.
    max_increase_per_second : int or float,
        The maximum allowed increase (per second). This value is extrapolated to the time resolution of records.
        This value must be positive.
    max_decrease_per_second : int or float
        The maximum allowed decrease (per second). This value is extrapolated to the time resolution of records.
        This value must be negative.

    Returns
    -------
    pd.DatetimeIndex
        Timestamps of outlier records.

    Notes
    -----
    In general, for temperatures, the decrease threshold is set less stringent than the increase
    threshold. This is because a temperature drop is meteorologically more
    common than a sudden increase, which is often the result of a radiation error.
    """

    # Validate argument values
    if max_decrease_per_second > 0:
        raise ValueError("max_decrease_per_second must be negative!")
    if max_increase_per_second < 0:
        raise ValueError("max_increase_per_second must be positive!")

    # Drop outliers from the series (these are NaNs)
    input_series = records.dropna()

    # Calculate timedelta between rows
    time_diff = input_series.index.to_series().diff()

    # Define filter
    step_filter = (
        # Step increase
        (
            (input_series - input_series.shift(1))
            > (float(max_increase_per_second) * time_diff.dt.total_seconds())
        )  # or
        |
        # Step decrease
        (
            (input_series - input_series.shift(1))
            < (max_decrease_per_second * time_diff.dt.total_seconds())
        )
    )

    logger.debug("Exiting function step_check")
    return step_filter[step_filter].index
