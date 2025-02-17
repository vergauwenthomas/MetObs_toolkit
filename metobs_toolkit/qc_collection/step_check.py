import logging
import pandas as pd

logger = logging.getLogger(__name__)


def step_check(
    records: pd.Series,
    max_increase_per_second: int | float,
    max_decrease_per_second: int | float,
) -> pd.DatetimeIndex:

    assert max_decrease_per_second < 0, f"max_decrease_per_second must be negative!"

    # drop outliers from the series (these are Nan's)
    input_series = records.dropna()

    # calculate timedelta between rows
    time_diff = input_series.index.to_series().diff()

    # define filter
    step_filter = (
        # step increase
        (
            (input_series - input_series.shift(1))
            > (float(max_increase_per_second) * time_diff.dt.total_seconds())
        )
        |  # or
        # step decrease
        (
            (input_series - input_series.shift(1))
            < (max_decrease_per_second * time_diff.dt.total_seconds())
        )
    )

    return step_filter[step_filter].index


# def step_check(obsdf, obstype, checks_settings):
#     """Test if observations do not produce spikes in timeseries.

#     Looking for jumps of the values of an observation type that are larger than
#     the limit specified in the qc_settings. These values are removed from the
#     input series and combined in the outlier df.

#     The purpose of this check is to flag observations with a value that is too
#     much different compared to the previous (not flagged) recorded value.

#     Schematically:

#     1. Iterate over all the stations.
#     2. Get the observations of the stations (i.g. drop the previously labeled outliers represented by NaN's).
#     3. Find the observations for which:

#        * The increase between two consecutive records is larger than the
#          threshold. This threshold is defined by a maximum increase per second
#          multiplied by the timedelta (in seconds) between the consecutive
#          records.
#        * Similar filter for a decrease.
#     4. The found observations are labeled as outliers.

#     Parameters
#     ------------
#     obsdf : pandas.DataFrame
#         The observations dataframe (Dataset.df) to check. Must have a triple
#         index (name, obstype, datetime).
#     obstype : str
#         The observation type to check for outliers.
#     outlierlabel: str
#         The label to use for flagged records.
#     checks_settings : dict
#         The dictionary containing the settings for the quality control checks.


#     Returns
#     ----------
#     obsdf : pandas.DataFrame()
#         The observations dataframe updated for step outliers. These are
#         represented by Nan values.
#     outl_df : pandas.DataFrame
#         The updated outliersdf.

#     Note
#     -----
#       In general, for temperatures, the decrease threshold is set less stringent than the increase
#       threshold. This is because a temperature drop is meteorologically more
#       common than a sudden increase which is often the result of a radiation error.

#     """
