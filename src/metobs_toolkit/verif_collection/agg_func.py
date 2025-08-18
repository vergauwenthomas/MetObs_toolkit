import numpy as np
import pandas as pd


def mean(records, **kwargs) -> float:
    return np.mean(records, **kwargs)

def median(records, **kwargs) -> float:
    return np.median(records, **kwargs)

def std(records, **kwargs) -> float:
    return np.std(records, **kwargs)

def count_notna(records) -> int:
    """
    Count the number of non-NaN elements in a pandas Series.

    Parameters
    ----------
    series : pd.Series
        The pandas Series to count non-NaN elements in.

    Returns
    -------
    int
        The count of non-NaN elements.
    """
    return pd.Series(records).notna().sum()
