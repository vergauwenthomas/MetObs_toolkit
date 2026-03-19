import pandas as pd
import logging

from metobs_toolkit.backend_collection.decorators import log_entry
from metobs_toolkit.qcresult import (
    pass_cond,
    flagged_cond,
    saved_cond,
    unchecked_cond,
    unmet_cond,
)

logger = logging.getLogger("<metobs_toolkit>")


@log_entry
def test_moving_window_condition(
    records: pd.Series, windowsize: pd.Timedelta, min_records_per_window: int
) -> bool:
    """
    Test if the resolution of the records meets the window constraints.

    Parameters
    ----------
    records : pd.Series
        Series with a datetime-like index.
    windowsize : pd.Timedelta
        Size of the moving window.
    min_records_per_window : int
        Minimum number of records required per window.

    Returns
    -------
    bool
        True if the minimum window members condition is met, False otherwise.

    Raises
    ------
    TypeError
        If any argument is not of the expected type.
    Exception
        If the input records do not have a perfectly regular timestamp.
    """

    # Get frequency of records
    freqstr = pd.infer_freq(records.index)
    if freqstr is None:
        raise Exception("The input records do not have a perfectly regular timestamp.")
    # Convert to timedelta
    # Note: sometimes 'h' is returned, and this gives issues, so add a 1 in front
    if not freqstr[0].isdigit():
        freqstr = "1" + freqstr

    freq = pd.Timedelta(freqstr)

    # Test if minimum window members condition is met
    ismet = (windowsize / freq) >= min_records_per_window
    logger.debug("Exiting function test_moving_window_condition.")
    return ismet


def create_qcresult_flags(
    all_input_records: pd.Series,
    unmet_cond_idx: pd.DatetimeIndex,
    outliers_before_white_idx: pd.DatetimeIndex,
    outliers_after_white_idx: pd.DatetimeIndex,
) -> pd.Series:
    """Create quality control flags series for all input records.

    This function generates a pandas Series containing QC flags for all timestamps
    in the input records. Records are categorized as: unchecked (NaN values),
    passed (valid non-outliers), unmet condition, saved (whitelisted outliers),
    or flagged (detected outliers).

    Parameters
    ----------
    all_input_records : pd.Series
        Complete series of records with datetime index to flag.
    unmet_cond_idx : pd.DatetimeIndex
        Timestamps where QC check conditions were not met.
    outliers_before_white_idx : pd.DatetimeIndex
        Timestamps of all detected outliers before whitelist filtering.
    outliers_after_white_idx : pd.DatetimeIndex
        Timestamps of outliers remaining after whitelist filtering.

    Returns
    -------
    pd.Series
        Series with same index as all_input_records containing QC flag strings:
        'unchecked', 'passed', 'condition_unmet', 'saved', or 'flagged'.
    """
    flags = pd.Series(data=unchecked_cond, index=all_input_records.index)
    flags.loc[all_input_records.dropna().index] = pass_cond
    flags.loc[unmet_cond_idx] = unmet_cond

    saved_records_idx = outliers_before_white_idx.difference(outliers_after_white_idx)
    flags.loc[saved_records_idx] = saved_cond
    flags.loc[outliers_after_white_idx] = flagged_cond

    return flags
