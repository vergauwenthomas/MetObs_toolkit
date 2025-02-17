import logging
import pandas as pd

logger = logging.getLogger(__name__)


def repetitions_check(records, max_N_repetitions):
    """Test if observation changes after a number of records.

    Looking for values of an observation type that are repeated at least with
    the frequency specified in the qc_settings. These values are labeled.

    Schematically:

    1. For each station, make a group of consecutive records for which
       the values do not change.
    2. Filter those groups that have more records than the maximum valid
       repetitions.
    3. All the records in these groups are labeled as outliers

    Parameters
    ------------
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
        The observations dataframe updated for repetitions outliers. These are
        represented by Nan values.
    outl_df : pandas.DataFrame
        The updated outliersdf.

    Note
    -----
      The repetitions check is similar to the persistence check, but not identical.
      The persistence check uses thresholds that are meteorologically based (i.g. the moving window is defined by a duration),
      in contrast to the repetitions check whose thresholds are instrumentally based (i.g. the "window" is defined by a number of records.)

    """

    # drop outliers from the series (these are Nan's)
    input_series = records.dropna()
    # Create group defenitions for repeting values that do not change
    persistence_filter = ((input_series.shift() != input_series)).cumsum()
    persdf = pd.DataFrame(
        data={"value": input_series, "persistgroup": persistence_filter},
        index=input_series.index,
    )

    # find outlier groups
    groups = persdf.groupby(["persistgroup"])
    # the above line groups the observations which have the same value and consecutive datetimes.
    group_sizes = groups.size()
    outlier_groups = group_sizes[group_sizes > max_N_repetitions]

    # combine all outlier groups
    outliers = pd.concat(
        [groups.get_group((outlgroup,)) for outlgroup in outlier_groups.index]
    )
    return outliers.index
