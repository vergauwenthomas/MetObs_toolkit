import logging
import pandas as pd

logger = logging.getLogger(__file__)
logging.basicConfig(level=logging.INFO)


def step_check(
    records: pd.Series,
    max_increase_per_second: int | float,
    max_decrease_per_second: int | float,
) -> pd.DatetimeIndex:
    """Check for 'spikes' and 'dips' in a timeseries.

    Test if observations do not produce spikes in timeseries. The maximum
    allowed increase and decrease per second is set in the argument,
    and is tested to each record (with respect to the previous record).

    If the difference between two consecutive records (i.e., the spike/dip) is larger than the
    threshold, the record is flagged as an outlier.

    Parameters
    ----------
    records : pd.Series
        A pandas Series containing the time series data to be checked. The index should be datetime-like.
    max_increase_per_second : int | float, >0
        The maximum allowed increase (per second). This value is extrapolated to the time resolution of records.
        This value must be positive!
    max_decrease_per_second : int | float, <0
        The maximum allowed decrease (per second). This value is extrapolated to the time resolution of records.
        This value must be negative!

    Returns
    -------
    pd.DatetimeIndex
        Timestamps of outlier records.

    Note
    -----
    In general, for temperatures, the decrease threshold is set less stringent than the increase
    threshold. This is because a temperature drop is meteorologically more
    common than a sudden increase which is often the result of a radiation error.
    """
    logger.info("Entering function step_check")

    # Validate argument types
    if not isinstance(records, pd.Series):
        raise TypeError("Argument 'records' must be of type pandas.Series")
    if not isinstance(max_increase_per_second, (int, float)):
        raise TypeError(
            "Argument 'max_increase_per_second' must be of type int or float"
        )
    if not isinstance(max_decrease_per_second, (int, float)):
        raise TypeError(
            "Argument 'max_decrease_per_second' must be of type int or float"
        )

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
        )
        |  # or
        # Step decrease
        (
            (input_series - input_series.shift(1))
            < (max_decrease_per_second * time_diff.dt.total_seconds())
        )
    )

    logger.info("Exiting function step_check")
    return step_filter[step_filter].index
