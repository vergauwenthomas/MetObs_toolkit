import logging
import pandas as pd

logger = logging.getLogger(__name__)


def test_moving_window_condition(
    records: pd.Series, windowsize: pd.Timedelta, min_records_per_window: int
) -> bool:
    """Test if the resolution of the records mets the window costrains"""
    # Get frequency of records
    freqstr = pd.infer_freq(records.index)
    if freqstr is None:
        raise Exception("The input records do not have a perfect frequent timestamp.")
    # Convert to timedelta
    freq = pd.Timedelta(freqstr)

    # Test if minimum window members condition is met
    ismet = (windowsize / freq) >= min_records_per_window
    return ismet
