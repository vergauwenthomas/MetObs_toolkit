import pandas as pd

# from metobs_toolkit.backend_collection.loggingmodule import logentry
import logging

from metobs_toolkit.backend_collection.loggingmodule import log_entry

logger = logging.getLogger(__name__)


possible_time_aggregates = [
    "year",
    "month",
    "hour",
    "minute",
    "second",
    "day_of_year",
    "season",
]


@log_entry
def get_time_derivates(datetimes) -> pd.DataFrame:
    """
    Construct a dataframe where all columns are time derivatives,
    and the index is the datetimeindex of self.df.

    Parameters
    ----------
    datetimes : pandas.DatetimeIndex
        The datetime index to derive time features from.

    Returns
    -------
    pd.DataFrame
        DataFrame with time derivative columns.
    """
    timesdf = pd.DataFrame(index=datetimes)
    for deriv in possible_time_aggregates:
        if deriv == "season":
            # custom method
            timesdf[deriv] = get_season(timesdf.index)
        else:
            timesdf[deriv] = getattr(timesdf.index, deriv)
    return timesdf


@log_entry
def get_season(datetimeindex: pd.DatetimeIndex) -> pd.Series:
    """
    Assign a season label to each date in a DatetimeIndex.

    Parameters
    ----------
    datetimeindex : pandas.DatetimeIndex
        The datetime index to assign seasons to.

    Returns
    -------
    pandas.Series
        Series of season labels indexed by datetimeindex.

    Raises
    ------
    TypeError
        If datetimeindex is not a pandas.DatetimeIndex.
    """

    summer_start = pd.Timestamp("2020-06-01")
    autumn_start = pd.Timestamp("2020-09-01")
    winter_start = pd.Timestamp("2020-12-01")
    spring_start = pd.Timestamp("2020-03-01")

    summer_end = autumn_start
    autumn_end = winter_start
    winter_end = spring_start
    spring_end = summer_start

    def _bin_season(startday, endday, seasonname):
        if endday > startday:
            return {(startday, endday): seasonname}
        else:
            # split in two bins
            return {(startday, 366): seasonname, (0, endday): seasonname}

    bindict = {}
    bindict.update(
        _bin_season(summer_start.day_of_year, summer_end.day_of_year, "summer")
    )
    bindict.update(
        _bin_season(autumn_start.day_of_year, autumn_end.day_of_year, "autumn")
    )
    bindict.update(
        _bin_season(winter_start.day_of_year, winter_end.day_of_year, "winter")
    )
    bindict.update(
        _bin_season(spring_start.day_of_year, spring_end.day_of_year, "spring")
    )

    seasonbins = pd.IntervalIndex.from_tuples(list(bindict.keys()))
    binlabels = list(bindict.values())

    seasons = pd.cut(
        x=datetimeindex.day_of_year,
        bins=seasonbins,
        right=False,
        include_lowest=True,
        # labels=binlabels, #ignored when bins is intervalindex
        retbins=False,
    ).map(
        dict(zip(seasonbins, binlabels)),
        na_action=None,
    )  # convert categories to seasonstrings
    # convert to series
    seasons = pd.Series(index=datetimeindex, data=seasons, name="season")

    return seasons
