import logging
import pandas as pd

logger = logging.getLogger(__name__)


def gross_value_check(
    records: pd.Series, lower_threshold: float, upper_threshold: float
) -> pd.DatetimeIndex:
    # drop nan's
    records = records.dropna()
    outliers = records[(records < lower_threshold) | (records > upper_threshold)].index
    return outliers


# def gross_value_check(obsdf, obstype, checks_settings):
#     """Filter out gross outliers from the observations.

#     Looking for values of an observation type that are not physical. These
#     values are labeled and the physical limits are specified in the qc_settings.

#     This check looks for outliers based on unrealistic values

#     1. Find observations that exceed a minimum and maximum value threshold.
#     2. These observations are labeled as outliers.

#     Parameters
#     ------------
#     obsdf : pandas.DataFrame
#         The observations dataframe (Dataset.df) to check. Must have a triple
#         index (name, obstype, datetime).
#     obstype : str
#         The observation type to check for outliers.
#     checks_settings : dict
#         The dictionary containing the settings for the quality control checks.


#     Returns
#     ----------
#     obsdf : pandas.DataFrame()
#         The observations dataframe updated for gross values. These are
#         represented by Nan values.
#     outl_df : pandas.DataFrame
#         The updated outliersdf.

#     """
