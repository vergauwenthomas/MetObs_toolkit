import logging
import pandas as pd
import numpy as np

from .common_functions import test_moving_window_condition

logger = logging.getLogger(__name__)


def persistence_check(
    records: pd.Series,
    timewindow: pd.Timedelta,
    min_records_per_window: int,
) -> pd.DatetimeIndex:
    """Test observations to change over a specific period.

    Looking for values of an observation type that do not change during a timewindow. These are flagged as outliers.

    In order to perform this check, at least N observations should be in that time window.

    Schematically:
    1. Find the stations that have a maximum assumed observation frequency
       that does not exceed the minimum number of records for moving window
       size. The window size is defined by a duration.
    2. Subset to those stations.
    3. For each station, a moving window scan is applied that validates if
       there is variation in the observations (NaN's are excluded). The
       validation is only applied when a sufficient amount of records are
       found in the window specified by a threshold.
    4. After the scan, all records found in the windows without variation
       are labeled as outliers.

    Parameters
    ------------
    station_frequencies : pandas.Series
        The frequencies of all the stations. This is a column in the metadf
        attribute of the Dataset.
    obsdf : pandas.DataFrame
        The observations dataframe (Dataset.df) to check. Must have a triple
        index (name, obstype, datetime).
    obstype : str
        The observation type to check for outliers.
    checks_settings : dict
        The dictionary containing the settings for the quality control checks.


    Returns
    ----------
    obsdf : pandas.DataFrame()
        The observations dataframe updated for persistence outliers. These are
        represented by Nan values.
    outl_df : pandas.DataFrame
        The updated outliersdf.

    """

    # test if the conditions for the moving window are met by the records frequency
    ismet = test_moving_window_condition(
        records=records,
        windowsize=timewindow,
        min_records_per_window=min_records_per_window,
    )
    if not ismet:
        logger.warning(
            "The minimum number of window members for the persistance check is not met!"
        )
        return pd.DatetimeIndex()

    # apply persistance
    def is_unique(
        window,
    ):  # comp order of N (while using the 'unique' function is Nlog(N))
        a = window.values
        a = a[~np.isnan(a)]
        return (a[0] == a).all()

    # TODO: This is very expensive if no coarsening is applied !!!! Can we speed this up?
    window_is_constant = (
        records.dropna()  # exclude outliers and gaps
        .rolling(
            window=timewindow,
            closed="both",
            center=True,
            min_periods=min_records_per_window,
        )
        .apply(is_unique)
    )
    # the returns are numeric values (0--> false, nan --> not checked(members/window condition not met), 1 --> outlier)
    window_is_constant = window_is_constant.map({0.0: False, np.nan: False, 1.0: True})

    return window_is_constant[window_is_constant].index
