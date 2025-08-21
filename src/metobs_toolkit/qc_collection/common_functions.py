import pandas as pd
import logging

from metobs_toolkit.backend_collection.loggingmodule import log_entry

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
