import logging
import pandas as pd

from metobs_toolkit.backend_collection.loggingmodule import log_entry

logger = logging.getLogger("<metobs_toolkit>")


@log_entry
def repetitions_check(records: pd.Series, max_N_repetitions: int) -> pd.DatetimeIndex:
    """
    Test if an observation changes after a number of repetitions.

    This function checks if the observation changes after a number of repetitions.
    If a value is repeated more than the specified number of times, all the repeated
    records are flagged as outliers.

    Be aware that the performance of this check depends on the `max_N_repetitions`
    and the time resolution of the observations.

    Parameters
    ----------
    records : pd.Series
        A pandas Series containing the time series data to be checked. The index should be datetime-like.
    max_N_repetitions : int
        The maximum number of repetitions allowed before the records are flagged as outliers.
        If the number of repetitions exceeds this value, all repeated records are flagged as outliers.

    Returns
    -------
    pd.DatetimeIndex
        Timestamps of outlier records.

    Notes
    -----
    The repetitions check is similar to the persistence check, but not identical.
    The persistence check uses thresholds that are meteorologically based (e.g., the moving window is defined by a duration),
    in contrast to the repetitions check whose thresholds are instrumentally based (e.g., the "window" is defined by a number of records).
    """

    # Drop outliers from the series (these are NaNs)
    input_series = records.dropna()

    # Create group definitions for repeating values that do not change
    persistence_filter = ((input_series.shift() != input_series)).cumsum()
    persdf = pd.DataFrame(
        data={"value": input_series, "persistgroup": persistence_filter},
        index=input_series.index,
    )

    # Find outlier groups
    groups = persdf.groupby(["persistgroup"])
    # The above line groups the observations which have the same value and consecutive datetimes.
    group_sizes = groups.size()
    outlier_groups = group_sizes[group_sizes > max_N_repetitions]

    # Combine all outlier groups
    if outlier_groups.empty:
        logger.debug("No outliers detected. Exiting repetitions_check function.")
        return pd.DatetimeIndex([])

    outliers = pd.concat(
        [
            groups.get_group(
                outlgroup,
            )
            for outlgroup in outlier_groups.index
        ]
    )
    logger.debug("Outliers detected. Exiting repetitions_check function.")
    return outliers.index
