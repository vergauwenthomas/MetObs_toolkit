import logging
from typing import Union

import pandas as pd
import numpy as np


logger = logging.getLogger("<metobs_toolkit>")


from metobs_toolkit.backend_collection.errorclasses import *
from metobs_toolkit.backend_collection.argumentcheckers import (
    fmt_timedelta_arg,
    fmt_datetime_arg,
)
import metobs_toolkit.plot_collection as plotting

from metobs_toolkit.dataset import Dataset
from metobs_toolkit.station import Station

possible_time_derivates = [
    "year",
    "month",
    "hour",
    "minute",
    "second",
    "day_of_year",
    "season",
]


class Analysis:
    def __init__(self, Dataholder: Union[Dataset, Station]):
        if isinstance(Dataholder, Dataset):
            df = Dataholder.df
            metadf = Dataholder.metadf
            obstypes = Dataholder.obstypes
        elif isinstance(Dataholder, Station):
            df = Dataholder.df
            # add the name to the index --> same structure as for Dataset.df
            df = (
                df.assign(name=Dataholder.name)
                .reset_index()
                .set_index(["datetime", "obstype", "name"])
            )
            metadf = Dataholder.metadf
            obstypes = {
                sens.obstype.name: sens.obstype
                for sens in Dataholder.sensordata.values()
            }
        else:
            raise TypeError(
                f"Dataholder is not a Dataset or Station, but a {type(Dataholder)}"
            )

        # set data
        self.df = df
        self.metadf = metadf

        # extra data
        self._obstypes = obstypes  # for displayin units in plots

    def __eq__(self, other):
        if not isinstance(other, Analysis):
            return False
        return self.df.equals(other.df) and self.metadf.equals(other.metadf)

    def get_info(self):
        """Print basic details of the Analysis instance."""
        print("Analysis Instance Information:")
        print(f"Number of records: {len(self.df)}")
        print(
            f"Observation types: {self.df.index.get_level_values('obstype').unique().tolist()}"
        )
        print(f"Stations: {self.df.index.get_level_values('name').unique().tolist()}")
        print(f"Available metadata columns: {self.metadf.columns.tolist()}")
        print(f"All knonw time-derivatives: {possible_time_derivates}")

    def get_tz(self):
        # Set timezone of 'fake_datetime' to match 'datetime' column (tz is displayed as xlabel)
        datetime_tz = self.df.index.get_level_values("datetime").tz
        return datetime_tz

    def aggregate_df(
        self, trgobstype="temp", agg=["lcz", "hour"], method=np.nanmean
    ) -> pd.DataFrame:
        """Aggregate all 'values' to specific groups, and return the dataframe."""
        # test if trgobstype is known
        self._obstype_is_known(trgobstype)

        # test if agg categories are valid
        for aggcat in agg:
            if aggcat not in self._all_possible_agg_categories():
                raise ValueError(
                    f"{aggcat} is not a possible agg category for {self}. These are all the possible agg categories: {self._all_possible_agg_categories()}."
                )

        # create a fulldf
        fulldf = self._combine_data()
        # Filter to obstype
        fulldf = fulldf.loc[fulldf["obstype"] == trgobstype]
        # Filter fulldf to relevant columns
        fulldf = fulldf[agg + ["value"]]

        # Aggregate the df
        agg_df = fulldf.groupby(agg, observed=True).agg(method)  # descrepation warning
        # sort index
        agg_df = agg_df.reset_index()
        agg_df = agg_df.set_index(agg)
        return agg_df

    def plot_diurnal_cycle(
        self,
        trgobstype="temp",
        colorby="name",
        title=None,
        ax=None,
        colordict=None,
        legend=True,
        return_data=False,
        figkwargs={},
    ):

        # test if trgobstype is known
        self._obstype_is_known(trgobstype)

        # test if colorby is valid
        if colorby not in self._all_possible_agg_categories():
            raise ValueError(
                f"{colorby} is not a possible agg category for {self}. These are all the possible agg categories: {self._all_possible_agg_categories()}."
            )

        colorby_blacklist = ["hour", "minute", "second"]
        if colorby in colorby_blacklist:
            raise ValueError(f"colorby cannot be one of {colorby_blacklist}")

        # create plotdf
        aggdf = self.aggregate_df(
            trgobstype=trgobstype,
            agg=[colorby, "hour", "minute", "second"],
            method=np.nanmean,
        )
        # create fake_datetime axes (fake because it is one day, but needed for representations of data sub hourly)
        aggdf = aggdf.reset_index()
        aggdf["fake_datetime"] = pd.to_datetime(
            "1990-01-01 "
            + aggdf["hour"].astype(str)
            + ":"
            + aggdf["minute"].astype(str)
            + ":"
            + aggdf["second"].astype(str)
        )
        aggdf["fake_datetime"] = aggdf["fake_datetime"].dt.tz_localize(self.get_tz())

        # construct plotdf (index: fake_datetime, columns: colorby-groups, values: value)
        plotdf = aggdf[[colorby, "fake_datetime", "value"]]
        plotdf = plotdf.set_index(["fake_datetime", colorby]).unstack(level=colorby)
        plotdf = plotdf["value"]  # selection by level0 of columns

        # Styling
        default_style = plotting.default_plot_settings["cycle_plot"]

        if ax is None:
            default_kwargs = {"figsize": default_style["figsize"]}
            default_kwargs.update(figkwargs)  # overwrite default with user defined
            ax = plotting.create_axes(**default_kwargs)

        # plot the cycles
        ax = plotting.make_diurnal_plot(
            plotdf=plotdf,
            colordict=colordict,
            ax=ax,
            refstation=None,
            figkwargs=figkwargs,
        )

        # Styling

        # Set ylabel
        obstype = self._obstypes[trgobstype]
        plotting.set_ylabel(ax, obstype._get_plot_y_label())

        # Set xlabel
        plotting.set_xlabel(ax, f"Diurnal Timestamps (in {self.get_tz()})")

        # Format timestamp ticks
        plotting.format_datetime_axes(ax, set_diurnal_format=True)

        # Add legend
        if legend:
            plotting.set_legend(ax, ncols=default_style["legend_n_columns"])

        # set title
        if title is None:
            title = f"Diurnal cycle of {obstype.name} grouped per {colorby}."
        plotting.set_title(ax, titlestr=title)

        if return_data:
            return ax, aggdf
        return ax

    def plot_diurnal_cycle_with_refernce_station(
        self,
        ref_station: str,
        trgobstype="temp",
        colorby="name",
        title=None,
        ax=None,
        colordict=None,
        legend=True,
        return_data=False,
        figkwargs={},
    ):

        # test if trgobstype is known
        self._obstype_is_known(trgobstype)

        # test if colorby is valid
        if colorby not in self._all_possible_agg_categories():
            raise ValueError(
                f"{colorby} is not a possible agg category for {self}. These are all the possible agg categories: {self._all_possible_agg_categories()}."
            )

        colorby_blacklist = ["hour", "minute", "second"]
        if colorby in colorby_blacklist:
            raise ValueError(f"colorby cannot be one of {colorby_blacklist}")

        # check if ref_station has records
        if ref_station not in self.df.index.get_level_values("name"):
            raise MetObsStationNotFound(
                f"No records are found for ref_station: {ref_station}"
            )

        # create a fulldf
        fulldf = self._combine_data()
        # Filter to obstype
        fulldf = fulldf.loc[fulldf["obstype"] == trgobstype]

        # Create reference records
        refrecords = fulldf[fulldf["name"] == ref_station][["datetime", "value"]]
        refrecords = refrecords.rename(columns={"value": "ref_value"})
        # merge refrecords on fulldf (by timestamp)
        fulldf = pd.merge(fulldf, refrecords, how="left", on="datetime")

        fulldf["diff_value"] = fulldf["value"] - fulldf["ref_value"]

        # Filter fulldf to relevant columns
        fulldf = fulldf[[colorby, "diff_value", "hour", "minute", "second"]]

        # Aggregate the df
        aggdf = fulldf.groupby(
            [colorby, "hour", "minute", "second"], observed=True
        ).agg(np.nanmean)
        # sort index
        aggdf = aggdf.reset_index()
        aggdf = aggdf.set_index(colorby)

        # create fake_datetime axes (fake because it is one day, but needed for representations of data sub hourly)
        aggdf = aggdf.reset_index()
        aggdf["fake_datetime"] = pd.to_datetime(
            "1990-01-01 "
            + aggdf["hour"].astype(str)
            + ":"
            + aggdf["minute"].astype(str)
            + ":"
            + aggdf["second"].astype(str)
        )
        aggdf["fake_datetime"] = aggdf["fake_datetime"].dt.tz_localize(self.get_tz())

        # construct plotdf (index: fake_datetime, columns: colorby-groups, values: value)
        plotdf = aggdf[[colorby, "fake_datetime", "diff_value"]]
        plotdf = plotdf.set_index(["fake_datetime", colorby]).unstack(level=colorby)
        plotdf = plotdf["diff_value"]  # selection by level0 of columns

        # Styling
        default_style = plotting.default_plot_settings["cycle_plot"]

        if ax is None:
            default_kwargs = {"figsize": default_style["figsize"]}
            default_kwargs.update(figkwargs)  # overwrite default with user defined
            ax = plotting.create_axes(**default_kwargs)

        # plot the cycles
        ax = plotting.make_diurnal_plot(
            plotdf=plotdf,
            colordict=colordict,
            ax=ax,
            refstation=ref_station,
            figkwargs=figkwargs,
        )

        # Styling
        # Set ylabel
        obstype = self._obstypes[trgobstype]
        plotting.set_ylabel(ax, obstype._get_plot_y_label())

        # Set xlabel
        plotting.set_xlabel(ax, f"Diurnal Timestamps (in {self.get_tz()})")

        # Format timestamp ticks
        plotting.format_datetime_axes(ax, set_diurnal_format=True)

        # Add legend
        if legend:
            plotting.set_legend(ax, ncols=default_style["legend_n_columns"])

        # set title
        if title is None:
            title = f"Diurnal cycle of {obstype.name} differnces with {ref_station} as reference, grouped per {colorby}."
        plotting.set_title(ax, titlestr=title)

        if return_data:
            return ax, aggdf
        return ax

    def subset_period(self, startdt, enddt):

        startdt = fmt_datetime_arg(startdt)
        enddt = fmt_datetime_arg(enddt)

        self.df = self.df.loc[
            (self.df.index.get_level_values("datetime") >= startdt)
            & (self.df.index.get_level_values("datetime") <= enddt)
        ]
        if self.df.empty:
            logger.warning(
                "The resulting DataFrame is empty after subsetting the period."
            )

        return self.df

    # ------------------------------------------
    #    Helpers
    # ------------------------------------------

    def _combine_data(self) -> pd.DataFrame:
        """Combines the df, metadf and time aggregations into one big dataframe"""

        # get time derivates
        timederivdf = self._get_time_derivates()
        # combine both on 'datatime'
        fulldf = pd.merge(
            left=self.df.reset_index(),
            right=timederivdf,
            how="left",
            left_on="datetime",
            right_index=True,
        )

        # combine data and metadata by merge on stationnames
        fulldf = pd.merge(left=fulldf, right=self.metadf, how="left", on="name")
        return fulldf

    def _obstype_is_known(self, trgobstype) -> None:
        """Test if the trgobstype has records, else Exception is raised"""
        if trgobstype in self.df.index.get_level_values("obstype"):
            return
        raise MetObsObstypeNotFound(f"{trgobstype} is not present in {self}.")

    def _get_time_derivates(self) -> pd.DataFrame:
        """Construct a dataframe where all columns are time derivatives,
        and the index is the datetimeindex of self.df.
        """
        datetimes = self.df.index.get_level_values("datetime")

        timesdf = pd.DataFrame(index=datetimes)
        for deriv in possible_time_derivates:
            if deriv == "season":
                # custom method
                timesdf[deriv] = get_season(timesdf.index)
            else:
                timesdf[deriv] = getattr(timesdf.index, deriv)
        return timesdf

    def _all_possible_agg_categories(self) -> list:
        """return all possible categories for aggregating"""
        metacategories = list(self.metadf.reset_index().columns)
        return list(set(metacategories + possible_time_derivates))


def get_season(datetimeindex):

    if not isinstance(datetimeindex, pd.DatetimeIndex):
        raise TypeError("The datetimeindex is not a pandas.DatetimeIndex.")
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
        _bin_season(autumn_start.day_of_year, autumn_end.day_of_year, "autmn")
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
        dict(zip(seasonbins, binlabels))
    )  # convert categories to seasonstrings
    # convert to series
    seasons = pd.Series(index=datetimeindex, data=seasons, name="season")

    return seasons
