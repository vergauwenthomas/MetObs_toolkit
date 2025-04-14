import logging  # Python default package
import pandas as pd  # Dependency package

logger = logging.getLogger(__file__)
logger.setLevel(logging.INFO)


def repetitions_check(records: pd.Series, max_N_repetitions: int) -> pd.DatetimeIndex:
    """Test if an observation changes after a number of repetitions.

    Perform a check that tests if the observation changes after a number of repetitions.
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
    in contrast to the repetitions check whose thresholds are instrumentally based (e.g., the "window" is defined by a number of records.)
    """

    # Validate argument types
    if not isinstance(records, pd.Series):
        raise TypeError("Argument 'records' must be of type 'pd.Series'.")
    if not isinstance(max_N_repetitions, int):
        raise TypeError("Argument 'max_N_repetitions' must be of type 'int'.")

    logger.info("Entering repetitions_check function.")

    # Drop outliers from the series (these are NaNs)
    input_series = records.dropna()

    # Create group definitions for repeating values that do not change
    persistence_filter = ((input_series.shift() != input_series)).cumsum()  # TYPO
    persdf = pd.DataFrame(  # TYPO
        data={"value": input_series, "persistgroup": persistence_filter},  # TYPO
        index=input_series.index,
    )

    # Find outlier groups
    groups = persdf.groupby(["persistgroup"])  # TYPO
    # The above line groups the observations which have the same value and consecutive datetimes.
    group_sizes = groups.size()
    outlier_groups = group_sizes[group_sizes > max_N_repetitions]

    # Combine all outlier groups
    if outlier_groups.empty:
        logger.info("No outliers detected. Exiting repetitions_check function.")
        return pd.DatetimeIndex([])

    outliers = pd.concat(
        [groups.get_group(outlgroup) for outlgroup in outlier_groups.index]  # TYPO
    )
    logger.info("Outliers detected. Exiting repetitions_check function.")
    return outliers.index
