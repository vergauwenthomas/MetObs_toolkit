import logging
from typing import Union, Tuple, Callable

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime as datetypeclass

from metobs_toolkit.backend_collection.errorclasses import *
from metobs_toolkit.backend_collection.argumentcheckers import (
    fmt_timedelta_arg,
    fmt_datetime_arg,
)
import metobs_toolkit.backend_collection.printing_collection as printing
import metobs_toolkit.plot_collection as plotting

from metobs_toolkit.dataset import Dataset
from metobs_toolkit.station import Station

logger = logging.getLogger("<metobs_toolkit>")

possible_time_aggregates = [  # TYPO
    "year",
    "month",
    "hour",
    "minute",
    "second",
    "day_of_year",
    "season",
]


class Analysis:
    """
    Analysis class for handling and analyzing meteorological data from a Dataset or Station.

    Parameters
    ----------
    Dataholder : Union[Dataset, Station]
        The data source to initialize the analysis object. It can either be a `Dataset`
        or a `Station` object.

    Attributes
    ----------
    fulldf : pandas.DataFrame
        A wide-format DataFrame containing the main data and time derivatives.
    metadf : pandas.DataFrame
        Metadata associated with the data.
    """

    def __init__(self, Dataholder: Union[Dataset, Station]):
        
        if not isinstance(Dataholder, (Dataset, Station)):
            raise TypeError(
                f"Dataholder is not a Dataset or Station, but a {type(Dataholder)}"
            )

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

        # Transform data to wide structure
        widedf = df.unstack("obstype")["value"]

        # set attributes for constructing the df
        self._df_cols = widedf.columns  # to reconstruct the df from fulldf

        # add time derivatives
        timederivdf = _get_time_derivates(
            datetimes=widedf.index.get_level_values("datetime").unique()
        )
        # combine both on 'datetime'
        fulldf = pd.merge(
            left=widedf.reset_index(),
            right=timederivdf,
            how="left",
            left_on="datetime",
            right_index=True,
        )
        self.fulldf = fulldf
        self.metadf = metadf

        # extra data
        self._obstypes = obstypes  # for displaying units in plots

    @property
    def df(self) -> pd.DataFrame:
        """
        Returns the full DataFrame without the time derivatives.

        Returns
        -------
        pd.DataFrame
            DataFrame indexed by ['datetime', 'name'] and containing observation columns.
        """
        logger.debug(f"Entering {self.__class__.__name__}.df property")
        return self.fulldf.set_index(["datetime", "name"])[self._df_cols]

    def __eq__(self, other: object) -> bool:
        """
        Check equality with another Analysis object.

        Parameters
        ----------
        other : object
            The object to compare with.

        Returns
        -------
        bool
            True if both the data and metadata are equal, False otherwise.
        """
        logger.debug(f"Entering {self.__class__.__name__}.__eq__")
        if not isinstance(other, Analysis):
            return False
        return self.df.equals(other.df) and self.metadf.equals(other.metadf)

    def get_info(self, printout: bool = True) -> Union[None, str]:
        """
        Provides information about the Analysis instance, including the number of records,
        observation types, metadata columns, station names, and known time derivatives.

        Parameters
        ----------
        printout : bool, optional
            If True, prints the information to the console. If False, returns the information
            as a string. Default is True.

        Returns
        -------
        None or str
            Returns None if `printout` is True. Returns the information string if `printout`
            is False.
        """
        logger.debug(f"Entering {self.__class__.__name__}.get_info")

        infostr = ""
        infostr += printing.print_fmt_title('General info of Analysis')
        infostr += printing.print_fmt_line(f"Number of records: {len(self.df)}")
        infostr += printing.print_fmt_line(f"Observation types: {list(self._df_cols)}")
        infostr += printing.print_fmt_line(f"Available metadata columns: {self.metadf.columns.tolist()}")
        infostr += printing.print_fmt_line(f"Stations: {self.fulldf['name'].unique().tolist()}")
        infostr += printing.print_fmt_line(f"All known time-derivatives: {possible_time_aggregates}")

        if printout:
            print(infostr)
        else:
            return infostr

    def get_tz(self) -> str:
        """
        Retrieve the timezone information of the 'datetime' index.

        Returns
        -------
        str
            The timezone of the 'datetime' index.
        """
        logger.debug(f"Entering {self.__class__.__name__}.get_tz")
        return self.df.index.get_level_values("datetime").tz

    def apply_filter_on_metadata(self, filter_string: str) -> "Analysis":
        """
        Apply a filter expression to the metadata and update the records accordingly.

        Parameters
        ----------
        filter_string : str
            A string representation of a filter expression to be applied on the metadata. For more
            details, see the documentation of pandas.DataFrame.query. The query is applied on
            the `.metadf` attribute.

        Warning
        -------
        This function modifies the data in place, so filtered-out data will be lost.

        Note
        -----
        A common error (UndefinedVariableError) is raised when a string value is not put
        in quotation marks.

        """
        logger.debug(f"Entering {self.__class__.__name__}.apply_filter_on_metadata")

        # Apply the filter string
        try:
            self.metadf.query(filter_string, inplace=True)
        except Exception as e:
            raise ValueError(f"Invalid filter string: {filter_string}. Error: {e}")

        # Subset the fulldf by the leftover stations
        self.fulldf = self.fulldf.loc[self.fulldf["name"].isin(self.metadf.index)]

        # Check if the filtered DataFrame is empty
        if self.fulldf.empty:
            logger.warning(
                "The resulting DataFrame is empty after applying the filter."
            )

    def apply_filter_on_records(self, filter_string: str) -> "Analysis":
        """
        Apply a filter expression to the records.

        Parameters
        ----------
        filter_string : str
            A string representation of a filter expression to be applied on the records. For more
            details, see the documentation of pandas.DataFrame.query. The query is applied on
            the `.fulldf` attribute.

        Warning
        -------
        This function modifies the data in place, so filtered-out data will be lost.

        Note
        -----
        A common error (UndefinedVariableError) is raised when a string value is not put
        in quotation marks.

        """
        logger.debug(f"Entering {self.__class__.__name__}.apply_filter_on_records")
        # Apply the filter string
        try:
            self.fulldf.query(filter_string, inplace=True)
        except Exception as e:
            raise ValueError(f"Invalid filter string: {filter_string}. Error: {e}")

        # Check if the filtered DataFrame is empty
        if self.fulldf.empty:
            logger.warning(
                "The resulting DataFrame is empty after applying the filter."
            )

    def subset_period(
        self,
        startdt: Union[pd.Timestamp, datetypeclass, str],
        enddt: Union[pd.Timestamp, datetypeclass, str],
    ) -> pd.DataFrame:
        """
        Subset the DataFrame to a specific time period.

        Parameters
        ----------
        startdt : Union[pd.Timestamp, datetypeclass, str]
            The start date and time of the desired period.
        enddt : Union[pd.Timestamp, datetypeclass, str]
            The end date and time of the desired period.

        Returns
        -------
        pd.DataFrame
            A DataFrame containing only the rows within the specified time period.

        Warning
        -------
        This function modifies the data in place, so filtered-out data will be lost.
        """
        logger.debug(f"Entering {self.__class__.__name__}.subset_period")
        startdt = fmt_datetime_arg(startdt)
        enddt = fmt_datetime_arg(enddt)

        # Subset the fulldf
        self.fulldf = self.fulldf[
            (self.fulldf["datetime"] >= startdt) & (self.fulldf["datetime"] <= enddt)
        ]
        if self.fulldf.empty:
            logger.warning(
                "The resulting DataFrame is empty after subsetting the period."
            )

    def aggregate_df(
        self,
        trgobstype: str = "temp",
        agg: list[str] = ["LCZ", "hour"],
        method: Callable = np.nanmean,
    ) -> pd.DataFrame:
        """
        Aggregate all 'values' to specific groups and return the resulting DataFrame.

        Parameters
        ----------
        trgobstype : str, optional
            The target observation type to aggregate, by default "temp".
        agg : list of str, optional
            List of column names to group by for aggregation, by default ["LCZ", "hour"].
        method : callable, optional
            Aggregation method to apply, by default `np.nanmean`.

        Returns
        -------
        pd.DataFrame
            A DataFrame with aggregated values, indexed by the specified grouping columns.
        """
        logger.debug(f"Entering {self.__class__.__name__}.aggregate_df")
    
        if not callable(method):
            raise TypeError("method must be callable.")

        # test if trgobstype is known
        self._obstype_is_known(trgobstype)

        # test if agg categories are valid
        for aggcat in agg:
            if aggcat not in self._all_possible_agg_categories():
                raise ValueError(
                    f"{aggcat} is not a possible agg category for {self}. These are all the possible agg categories: {self._all_possible_agg_categories()}."
                )

        # create a fulldf
        fulldf = pd.merge(left=self.fulldf, right=self.metadf, how="left", on="name")
        # Filter fulldf to relevant columns
        fulldf = fulldf[agg + [trgobstype]]

        # Aggregate the df
        agg_df = fulldf.groupby(agg, observed=True).agg(method)
        # sort index
        agg_df = agg_df.reset_index()
        agg_df = agg_df.set_index(agg)
        return agg_df

    def plot_diurnal_cycle(
        self,
        trgobstype: str = "temp",
        colorby: str = "name",
        title: Union[str, None] = None,
        ax: Union[plt.Axes, None] = None,
        colordict: Union[dict, None] = None,
        legend: bool = True,
        return_data: bool = False,
        figkwargs: dict = {},
    ) -> Union[plt.Axes, Tuple[plt.Axes, pd.DataFrame]]:
        """
        Plot the diurnal cycle of a specified observation type, grouped by a given category.

        Parameters
        ----------
        trgobstype : str, optional
            The target observation type to plot (e.g., "temp"). Default is "temp".
        colorby : str, optional
            The category by which to group and color the data (e.g., "name"). Default is "name".
        title : str or None, optional
            The title of the plot. If None, a default title is generated. Default is None.
        ax : matplotlib.axes.Axes or None, optional
            The axes on which to plot. If None, a new figure and axes are created. Default is None.
        colordict : dict or None, optional
            A dictionary mapping group names to colors. If None, default colors are used. Default is None.
        legend : bool, optional
            Whether to include a legend in the plot. Default is True.
        return_data : bool, optional
            If True, returns the plot axes and the aggregated data used for plotting. Default is False.
        figkwargs : dict, optional
            Additional keyword arguments for figure creation. Default is an empty dictionary.

        Returns
        -------
        matplotlib.axes.Axes or tuple
            The plot axes. If `return_data` is True, returns a tuple of the plot axes and the aggregated DataFrame.
        """
        logger.debug(f"Entering {self.__class__.__name__}.plot_diurnal_cycle")
      

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
        plotdf = aggdf[[colorby, "fake_datetime", trgobstype]]
        plotdf = plotdf.set_index(["fake_datetime", colorby]).unstack(level=colorby)
        plotdf = plotdf[trgobstype]  # selection by level0 of columns

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

    def plot_diurnal_cycle_with_reference_station(
        self,
        ref_station: str,
        trgobstype: str = "temp",
        colorby: str = "name",
        title: Union[str, None] = None,
        ax: Union[plt.Axes, None] = None,
        colordict: Union[dict, None] = None,
        legend: bool = True,
        return_data: bool = False,
        figkwargs: dict = {},
    ) -> Union[plt.Axes, Tuple[plt.Axes, pd.DataFrame]]:
        """
        Plot the diurnal cycle of differences between observations and a reference station.

        Parameters
        ----------
        ref_station : str
            The name of the reference station to compare against.
        trgobstype : str, optional
            The observation type to analyze (e.g., "temp"). Default is "temp".
        colorby : str, optional
            The column name to group data by for coloring the plot. Default is "name".
        title : str or None, optional
            The title of the plot. If None, a default title is generated. Default is None.
        ax : matplotlib.axes.Axes or None, optional
            The axes object to plot on. If None, a new figure and axes are created. Default is None.
        colordict : dict or None, optional
            A dictionary mapping group names (from `colorby`) to colors. Default is None.
        legend : bool, optional
            Whether to include a legend in the plot. Default is True.
        return_data : bool, optional
            If True, returns the aggregated data used for plotting along with the axes. Default is False.
        figkwargs : dict, optional
            Additional keyword arguments for figure creation. Default is an empty dictionary.

        Returns
        -------
        Union[matplotlib.axes.Axes, Tuple[matplotlib.axes.Axes, pandas.DataFrame]]
            The axes object with the plot. If `return_data` is True, also returns the aggregated DataFrame.

        Notes
        -----
        The method computes the differences between the target observation type (`trgobstype`)
        and the reference station's values, then aggregates the data by the specified grouping
        (`colorby`) and time components (hour, minute, second). The resulting diurnal cycle
        is plotted, with options for customization.
        """
        logger.debug(f"Entering {self.__class__.__name__}.plot_diurnal_cycle_with_reference_station")
      
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

        # create a differences in time
        fulldf = pd.merge(
            left=self.fulldf, right=self.metadf, how="left", on="name"
        ).set_index(["datetime", "name"])[
            [trgobstype, colorby, "hour", "minute", "second"]
        ]

        refrecords = fulldf.xs(ref_station, level="name", drop_level=True)[trgobstype]
        fulldf["diff_value"] = fulldf[trgobstype] - refrecords

        # Aggregate the df
        aggdf = fulldf.groupby(
            [colorby, "hour", "minute", "second"], observed=True
        ).agg(np.nanmean)

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
        plotting.set_ylabel(ax, f"{obstype._get_plot_y_label()} difference")

        # Set xlabel
        plotting.set_xlabel(ax, f"Diurnal Timestamps (in {self.get_tz()})")

        # Format timestamp ticks
        plotting.format_datetime_axes(ax, set_diurnal_format=True)

        # Add legend
        if legend:
            plotting.set_legend(ax, ncols=default_style["legend_n_columns"])

        # set title
        if title is None:
            title = f"Diurnal cycle of {obstype.name} differences with {ref_station} as reference, grouped per {colorby}."
        plotting.set_title(ax, titlestr=title)

        if return_data:
            return ax, aggdf
        return ax

    # ------------------------------------------
    #    Helpers
    # ------------------------------------------

    def _obstype_is_known(self, trgobstype: str) -> None:
        """
        Test if the trgobstype has records, else Exception is raised.

        Parameters
        ----------
        trgobstype : str
            The observation type to check.

        Raises
        ------
        MetObsObstypeNotFound
            If the observation type is not present.
        """
        logger.debug(f"Entering {self.__class__.__name__}._obstype_is_known")
        if trgobstype in self.df.columns:
            return
        raise MetObsObstypeNotFound(f"{trgobstype} is not present in {self}.")

    def _all_possible_agg_categories(self) -> list:
        """
        Return all possible categories for aggregating.

        Returns
        -------
        list
            List of possible aggregation categories.
        """
        logger.debug(f"Entering {self.__class__.__name__}._all_possible_agg_categories")
        metacategories = list(self.metadf.reset_index().columns)
        return list(set(metacategories + possible_time_aggregates))


# ------------------------------------------
#    Helping methods
# ------------------------------------------
def _get_time_derivates(datetimes) -> pd.DataFrame:
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
    logger.debug("Entering _get_time_derivates function.")
    timesdf = pd.DataFrame(index=datetimes)
    for deriv in possible_time_aggregates:
        if deriv == "season":
            # custom method
            timesdf[deriv] = get_season(timesdf.index)
        else:
            timesdf[deriv] = getattr(timesdf.index, deriv)
    return timesdf


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
    logger.debug("Entering get_season function.")

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
        dict(zip(seasonbins, binlabels))
    )  # convert categories to seasonstrings
    # convert to series
    seasons = pd.Series(index=datetimeindex, data=seasons, name="season")

    return seasons
