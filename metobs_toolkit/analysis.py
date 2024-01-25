#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module contains the Analysis class and all its methods.

A Analysis holds a set of 'good' observations and the methods will analyse it.
"""
from datetime import datetime
import pandas as pd
import numpy as np
import logging
import copy
from scipy.stats import pearsonr

from metobs_toolkit.plotting_functions import (
    cycle_plot,
    heatmap_plot,
    correlation_scatter,
)

from metobs_toolkit.df_helpers import (
    datetime_subsetting,
    subset_stations,
    fmt_datetime_argument,
)

logger = logging.getLogger(__name__)


class Analysis:
    """The Analysis class contains methods for analysing observations."""

    def __init__(self, obsdf, metadf, settings, data_template):
        """Initialize an Analysis."""
        self.df = obsdf
        self.metadf = metadf
        self.settings = settings
        self.data_template = data_template

        # analysis objects
        self.lc_cor_dict = {}
        self._lc_cor_obstype = None
        self._lc_groupby_labels = None

        # add empty lcz column to metadf if it is not present
        if "lcz" not in self.metadf.columns:
            self.metadf["lcz"] = np.nan

    def __str__(self):
        """Print a overview of the analysis."""
        if self.df.empty:
            return "Empty Analysis instance."
        add_info = ""
        n_stations = self.df.index.get_level_values("name").unique().shape[0]
        n_obs_tot = self.df.shape[0]

        startdt = self.df.index.get_level_values("datetime").min()
        enddt = self.df.index.get_level_values("datetime").max()

        if (not self.metadf["lat"].isnull().all()) & (
            not self.metadf["lon"].isnull().all()
        ):
            add_info += "     *Coordinates are available for all stations. \n"

        if not self.metadf["lcz"].isnull().all():
            add_info += "     *LCZ's are available for all stations. \n"

        if bool(self.lc_cor_dict):
            add_info += f"     *landcover correlations are computed on group: {self._lc_groupby_labels}  \n"

        return (
            f"Analysis instance containing: \n \
    *{n_stations} stations \n \
    *{self.df.columns.to_list()} observation types \n \
    *{n_obs_tot} observation records \n{add_info} \n \
    *records range: {startdt} --> {enddt} (total duration:  {enddt - startdt})"
            + add_info
        )

    def __repr__(self):
        """Print a overview of the analysis."""
        return self.__str__()

    # =============================================================================
    #     Setters
    # =============================================================================

    def subset_period(self, startdt, enddt):
        """Subset the observations of the Analysis to a specific period.

        The same timezone is assumed as the data.

        Parameters
        ----------
        startdt : datetime.datetime
            The start datetime to filter the observations to.
        enddt : datetime.datetime
            The end datetime to filter the observations to.

        Returns
        -------
        None.

        Note
        --------
        If a timezone unaware datetime is given as an argument, it is interpreted
        as if it has the same timezone as the observations.
        """
        if not isinstance(startdt, type(datetime(2020, 1, 1))):
            logger.info(f" {startdt} not a datetime type. Ignore subsetting!")
            return
        if not isinstance(enddt, type(datetime(2020, 1, 1))):
            logger.info(f" {enddt} not a datetime type. Ignore subsetting!")
            return

        startdt = fmt_datetime_argument(
            startdt, self.settings.time_settings["timezone"]
        )
        enddt = fmt_datetime_argument(enddt, self.settings.time_settings["timezone"])

        self.df = datetime_subsetting(self.df, startdt, enddt)

    # =============================================================================
    #   Helpers
    # =============================================================================

    def apply_filter(self, expression):
        """Filter an Analysis by a user definde string expression.

        This can be used to filter the observation to specific meteorological
        conditions (i.e. low windspeeds, high humidity, cold temperatures, ...)

        The filter expression contains only columns present in the Analysis.df
        and/or the Analysis.metadf.

        A New Analysis object is returned.

        Parameters
        ----------
        expression : str
            A filter expression using columnnames present in either df or metadf.
            The following timestamp derivatives can be used as well: [minute, hour,
            month, year, day_of_year, week_of_year, season]. The quarry_str may
            contain number and expressions like <, >, ==, >=, \*, +, .... Multiple filters
            can be combine to one expression by using & (AND) and | (OR).

        Returns
        -------
        filtered_analysis : metobs_toolkit.Analysis
            The filtered Analysis.


        Note
        -------
        All timestamp derivative values are numeric except for 'season',
        possible values are ['winter', 'spring', 'summer', 'autumn'].

        Note
        ------
        Make shure to use \" of \' to indicate string values in the expression if
        needed.

        """
        child_df, child_metadf = filter_data(
            df=self.df, metadf=self.metadf, quarry_str=expression
        )

        return Analysis(
            obsdf=child_df,
            metadf=child_metadf,
            settings=self.settings,
            data_template=self.data_template,
        )

    def aggregate_df(self, df=None, agg=["lcz", "hour"], method="mean"):
        """Aggregate observations to a (list of) categories.

        The output will be a dataframe that is aggregated to one, or more
        categories. A commen example is aggregating to LCZ's.


        Parameters
        ----------
        df : pandas.DataFrame or None
            The observations to aggregate. If None, the df attribute of the
            Analysis instance is used. The default is None.
        agg : list, optional
            The list of columnnames to aggregate to. If 'lcz' is included, the
            lcz information is extracted from the Analysis.metadf. The default
            is ['lcz', 'datetime'].
        method : str, optional
            list of functions and/or function names, e.g. [np.sum, 'mean']. The
            default is 'mean'.

        Returns
        -------
        pandas.DataFrame
            A dataframe with the agg columns as an index. The values are the
            aggregated values.

        Note
        -------
        Present columns that ar non-numeric and are not in the agg list, are
        not present in the return, since these values cannot be aggregated.

        """
        if df is None:
            df = copy.deepcopy(self.df)
        df = df.reset_index()

        time_agg_keys = [
            "minute",
            "hour",
            "month",
            "year",
            "day_of_year",
            "week_of_year",
            "season",
        ]

        # scan trough the metadf for aggregation keys
        for agg_key in agg:
            if agg_key not in df.columns:
                # look in metadf
                if agg_key in self.metadf.columns:
                    df = pd.merge(
                        df,
                        self.metadf[[agg_key]],
                        how="left",
                        left_on="name",
                        right_index=True,
                    )

        # Check if all agg keys are present or defined:
        possible_agg_keys = time_agg_keys
        possible_agg_keys.extend(list(df.columns))
        unmapped = [agg_key for agg_key in agg if agg_key not in possible_agg_keys]
        assert len(unmapped) == 0, f"cannot aggregate to unknown labels: {unmapped}."

        # make time-derivate columns if required
        df = _make_time_derivatives(df, agg)

        # check if not all values are Nan
        for agg_name in agg:
            assert (
                not df[agg_name].isnull().all()
            ), f"Aggregation to {agg_name} not possible because no valid values found for {agg_name}."

        # remove datetime column if present, because no aggregation can be done on
        # datetime and it gives a descrepation warning
        if "datetime" in df.columns:
            df = df.drop(columns=["datetime"])

        # Remove name column if present and not in the aggregation scheme,
        # this happens because name was in the index
        if "name" not in agg:
            df = df.drop(columns=["name"], errors="ignore")

        # Aggregate the df
        agg_df = df.groupby(agg).agg(method, numeric_only=True)  # descrepation warning
        # sort index
        agg_df = agg_df.reset_index()
        agg_df = agg_df.set_index(agg)
        return agg_df

    # =============================================================================
    #   Analyse method
    # =============================================================================
    def get_anual_statistics(
        self,
        groupby=["name"],
        obstype="temp",
        agg_method="mean",
        stations=None,
        startdt=None,
        enddt=None,
        plot=True,
        errorbands=False,
        title=None,
        y_label=None,
        legend=True,
        _return_all_stats=False,
    ):
        """
        Create an anual cycle for aggregated groups.

        (In the plot, unique combination of groupby categories is presented
         as a line.)

        Parameters
        ----------
        groupby : list string, optional
            Variables to aggregate to. These can be columns in the metadf, or
            time aggregations ('hour', 'year', 'week_of_year', ...]. 'name' will
            aggregate to the stationnames. The default is ['name'].
        obstype : str, optional
            Element of the metobs_toolkit.observation_types The default is 'temp'.
        agg_method : str, optional
            Function names to use for aggregation, e.g. [np.sum, 'mean']. The
            default is 'mean'.
        stations : list, optional
             List of station names to use. If None, all present stations will be used. The default is None.
        startdt : datetime.datetime, optional
            The start datetime of the observations to use. If None, all timestamps will be used. The default is None.
        enddt : datetime.datetime, optional
            The end datetime of the observations to use. If None, all timestamps will be used. The default is None.
        plot : bool, optional
            If True, an anual plot is made. The default is True.
        errorbands : bool, optional
             If True, the std is representd in the plot by colored bands. The default is False.
        title : string, optional
             Title of the figure, if None a default title is generated. The default is None.
        y_label : string, optional
             y-axes label of the figure, if None a default label is generated. The default is None.
        legend : bool, optional
             I True, a legend is added to the plot. The default is True.

        Returns
        -------
        df : pandas.DataFrame()
            The dataframe containing the aggregated values.

        Note
        --------
        If a timezone unaware datetime is given as an argument, it is interpreted
        as if it has the same timezone as the observations.

        """
        # title
        desc_dict = self.data_template[obstype].to_dict()

        if "description" not in desc_dict:
            desc_dict["description"] = obstype
        if not isinstance(desc_dict["description"], str):
            desc_dict["description"] = obstype

        if title is None:
            title = f'Anual {desc_dict["description"]} cycle plot per {groupby}.'
        else:
            title = str(title)

        # ylabel
        if y_label is None:
            if "units" not in desc_dict:
                y_label = f'{desc_dict["description"]} (units unknown)'
            else:
                y_label = f'{desc_dict["description"]} ({desc_dict["units"]})'
        else:
            y_label = str(y_label)

        stats = self.get_aggregated_cycle_statistics(
            obstype=obstype,
            stations=stations,
            aggregation=groupby,
            aggregation_method=agg_method,
            horizontal_axis="month",
            startdt=startdt,
            enddt=enddt,
            plot=plot,
            title=title,
            y_label=y_label,
            legend=legend,
            errorbands=errorbands,
            verbose=_return_all_stats,
        )
        return stats

    def get_diurnal_statistics(
        self,
        colorby="name",
        obstype="temp",
        stations=None,
        startdt=None,
        enddt=None,
        plot=True,
        title=None,
        y_label=None,
        legend=True,
        errorbands=False,
        _return_all_stats=False,
    ):
        """
        Create an average diurnal cycle for the observations.

        (In the plot, each station is represed by a line.)


        Parameters
        ----------
        colorby : 'name' or 'lcz', optional
            If 'name' the plotted lines will be colored per station, if 'lcz' the colors represent the stations lcz. The default is 'name'.
        obstype : str, optional
            Element of the metobs_toolkit.observation_types The default is 'temp'.
        stations : list, optional
            List of station names to use. If None, all present stations will be used. The default is None.
        startdt : datetime.datetime, optional
            The start datetime of the observations to use. If None, all timestamps will be used. The default is None.
        enddt : datetime.datetime, optional
            The end datetime of the observations to use. If None, all timestamps will be used. The default is None.
        plot : bool, optional
            If True, an diurnal plot is made. The default is True.
        title : string, optional
             Title of the figure, if None a default title is generated. The default is None.
        y_label : string, optional
             y-axes label of the figure, if None a default label is generated. The default is None.
        legend : bool, optional
             I True, a legend is added to the plot. The default is True.
        errorbands : bool, optional
            If True, the std is representd in the plot by colored bands. The default is False.

        Returns
        -------
        df : pandas.DataFrame()
            The dataframe containing the aggregated values.

        Note
        --------
        If a timezone unaware datetime is given as an argument, it is interpreted
        as if it has the same timezone as the observations.

        """
        # title
        desc_dict = self.data_template[obstype].to_dict()

        if "description" not in desc_dict:
            desc_dict["description"] = obstype
        if not isinstance(desc_dict["description"], str):
            desc_dict["description"] = obstype

        if title is None:
            if startdt is None:
                if enddt is None:
                    title = f"Hourly average {obstype} diurnal cycle"
                else:
                    title = f"Hourly average {obstype} diurnal cycle until {enddt}"
            else:
                if enddt is None:
                    title = f"Hourly average {obstype} diurnal cycle from {startdt}"
                else:
                    title = f"Hourly average {obstype} diurnal cycle for period {startdt} - {enddt}"

        else:
            title = str(title)

        # ylabel
        if y_label is None:
            if "units" not in desc_dict:
                y_label = f'{desc_dict["description"]} (units unknown)'
            else:
                y_label = f'{desc_dict["description"]} ({desc_dict["units"]})'
        else:
            y_label = str(y_label)

        stats = self.get_aggregated_cycle_statistics(
            obstype=obstype,
            stations=stations,
            aggregation=[colorby],
            aggregation_method="mean",
            horizontal_axis="hour",
            startdt=startdt,
            enddt=enddt,
            plot=plot,
            title=title,
            y_label=y_label,
            legend=legend,
            errorbands=errorbands,
            verbose=_return_all_stats,
        )
        return stats

    def get_diurnal_statistics_with_reference(
        self,
        refstation,
        colorby="name",
        obstype="temp",
        tollerance="30T",
        stations=None,
        startdt=None,
        enddt=None,
        plot=True,
        title=None,
        y_label=None,
        legend=True,
        errorbands=False,
        show_zero_horizontal=True,
        _return_all_stats=False,
    ):
        """
        Create an average diurnal cycle for the observation differences of a reference station.

        All observational values are converted to differences with the closest
        (in time) reference observation. No reference observation is found when
        the time difference is larger than the tollerance.

        (In the plot, each station is represed by a line.)

        Parameters
        ----------
        refstation : str,
            Name of the station to use as a reference.
        colorby : 'name' or 'lcz', optional
            If 'name' the plotted lines will be colored per station, if 'lcz' the colors represent the stations lcz. The default is 'name'.
        obstype : str, optional
            Element of the metobs_toolkit.observation_types The default is 'temp'.
        tollerance : Timedelta or str, optional
            The tollerance string or object representing the maximum translation in time to find a reference
            observation for each observation. Ex: '5T' is 5 minutes, '1H', is one hour. The default is '30T'.
        stations : list, optional
            List of station names to use. If None, all present stations will be used. The default is None.
        startdt : datetime.datetime, optional
            The start datetime of the observations to use. If None, all timestamps will be used. The default is None.
        enddt : datetime.datetime, optional
            The end datetime of the observations to use. If None, all timestamps will be used. The default is None.
        plot : bool, optional
            If True, a diurnal plot is made. The default is True.
        title : string, optional
             Title of the figure, if None a default title is generated. The default is None.
        y_label : string, optional
             y-axes label of the figure, if None a default label is generated. The default is None.
        legend : bool, optional
             I True, a legend is added to the plot. The default is True.
        errorbands : bool, optional
            If True, the std is representd in the plot by colored bands. The upper bound represents +1 x std, the lower bound -1 x std. The default is False.
        show_zero_horizontal : bool, optional
            If True a horizontal line is drawn in the plot at zero. The default is True.

        Returns
        -------
        df : pandas.DataFrame()
            The dataframe containing the aggregated values.

        Note
        --------
        If a timezone unaware datetime is given as an argument, it is interpreted
        as if it has the same timezone as the observations.

        """
        obsdf = self.df
        obsdf = obsdf[obstype].reset_index()

        # extract refernce from observations
        refdf = obsdf[obsdf["name"] == refstation]
        obsdf = obsdf[obsdf["name"] != refstation]

        assert (
            not refdf.empty
        ), f"Error: No reference observation found (after filtering) for {refstation}"
        assert not obsdf.empty, "Error: No observation found (after filtering)"

        # Syncronize observations with the reference observations
        refdf = refdf.rename(
            columns={obstype: "ref_" + obstype, "datetime": "ref_datetime"}
        )
        mergedf = pd.merge_asof(
            left=obsdf.sort_values("datetime"),
            right=refdf[["ref_datetime", "ref_" + obstype]].sort_values("ref_datetime"),
            right_on="ref_datetime",
            left_on="datetime",
            direction="nearest",
            tolerance=pd.Timedelta(tollerance),
        )

        # Get differnces
        mergedf["temp"] = mergedf["temp"] - mergedf["ref_temp"]

        # Subset to relavent columns
        mergedf = mergedf.reset_index()
        mergedf = mergedf[["name", "datetime", obstype]]
        mergedf = mergedf.set_index(["name", "datetime"])

        # title
        desc_dict = self.data_template[obstype].to_dict()
        if "description" not in desc_dict:
            desc_dict["description"] = obstype
        if not isinstance(desc_dict["description"], str):
            desc_dict["description"] = obstype

        if title is None:
            if startdt is None:
                if enddt is None:
                    title = f"Hourly average {obstype} diurnal cycle, with {refstation} as reference,"
                else:
                    title = f"Hourly average {obstype} diurnal cycle, with {refstation} as reference, until {enddt}"
            else:
                if enddt is None:
                    title = f"Hourly average {obstype} diurnal cycle, with {refstation} as reference, from {startdt}"
                else:
                    title = f"Hourly average {obstype} diurnal cycle, with {refstation} as reference, for period {startdt} - {enddt}"

        else:
            title = str(title)

        # ylabel
        if y_label is None:
            if "units" not in desc_dict:
                y_label = f'{desc_dict["description"]} (units unknown)'
            else:
                y_label = f'{desc_dict["description"]} ({desc_dict["units"]})'
        else:
            y_label = str(y_label)

        stats = self.get_aggregated_cycle_statistics(
            obstype=obstype,
            stations=stations,
            aggregation=[colorby],
            aggregation_method="mean",
            horizontal_axis="hour",
            startdt=startdt,
            enddt=enddt,
            plot=plot,
            title=title,
            y_label=y_label,
            legend=legend,
            errorbands=errorbands,
            verbose=_return_all_stats,
            _obsdf=mergedf,
            _show_zero_line=show_zero_horizontal,
        )
        return stats

    def get_aggregated_cycle_statistics(
        self,
        obstype="temp",
        aggregation=["lcz", "datetime"],
        aggregation_method="mean",
        horizontal_axis="hour",
        stations=None,
        startdt=None,
        enddt=None,
        plot=True,
        title=None,
        y_label=None,
        legend=True,
        errorbands=False,
        verbose=False,
        _obsdf=None,
        _show_zero_line=False,
    ):
        """Create an average cycle for an aggregated categorie.

        A commen example is to aggregate to the LCZ's, so to get the diurnal
        cycle per LCZ rather than per station.

        (In the plot, each aggregated category different from datetime, is represed by a line.)

        Parameters
        ----------
        obstype : str, optional
            Element of the metobs_toolkit.observation_types The default is 'temp'.
        aggregation : list, optional
            List of variables to aggregate to. These variables should either a
            categorical observation type, a categorical column in the metadf or
            a time aggregation. All possible time aggreagetions are: ['minute',
            'hour', 'month', 'year', 'day_of_year',
            'week_of_year', 'season']. The default is ['lcz', 'datetime'].
        aggregation_method : str, optional
            Which (numpy) function is used to aggregate the observations. The default is 'mean'.
        horizontal_axis : str, optional
            Which aggregated value will be represented on the horizontal axis
            of the plot. The default is 'hour'.
        stations : list, optional
            List of station names to use. If None, all present stations will be used. The default is None.
        startdt : datetime.datetime, optional
            The start datetime of the observations to use. If None, all timestamps will be used. The default is None.
        enddt : datetime.datetime, optional
            The end datetime of the observations to use. If None, all timestamps will be used. The default is None.
        plot : bool, optional
            If True, a diurnal plot is made. The default is True.
        title : string, optional
             Title of the figure, if None a default title is generated. The default is None.
        y_label : string, optional
             y-axes label of the figure, if None a default label is generated. The default is None.
        legend : bool, optional
             I True, a legend is added to the plot. The default is True.
        errorbands : bool, optional
            If True, the std is representd in the plot by colored bands. The default is False.
        verbose : True, optional
            If True, an additional dataframe with aggregation information is returned . The default is False.

        Returns
        -------
        df : pandas.DataFrame()
            The dataframe containing the aggregated values.

        Note
        -------
        If a timezone unaware datetime is given as an argument, it is interpreted
        as if it has the same timezone as the observations.

        """
        if _obsdf is None:
            obsdf = self.df[[obstype]]
        else:
            obsdf = _obsdf

        assert not obsdf.empty, f"Error: No observations in the analysis.df: {self.df}"
        # Filter stations
        if stations is not None:
            if isinstance(stations, str):
                stations = [stations]

            obsdf = subset_stations(obsdf, stations)
        assert (
            not obsdf.empty
        ), f"Error: No more observations after subsetting to {stations}"

        # Filter datetimes
        obsdf = datetime_subsetting(df=obsdf, starttime=startdt, endtime=enddt)
        assert (
            not obsdf.empty
        ), f"Error: No more observations after subsetting to {startdt} and {enddt}"

        startdt = obsdf.index.get_level_values("datetime").min()
        enddt = obsdf.index.get_level_values("datetime").max()

        # add hour to aggregation (will be the x-axis)
        if horizontal_axis not in aggregation:
            aggregation.insert(0, horizontal_axis)

        # add other methods for errorbands and stats
        methods = ["mean", "std", "median"]
        methods.append(aggregation_method)
        methods = list(set(methods))

        # compute the aggregation statistics
        aggdf = self.aggregate_df(df=obsdf, agg=aggregation, method=methods)

        # since only one observation type is in the stats, drop the column
        # level with the obstye, this is not relevant
        aggdf = aggdf.droplevel(0, axis="columns")

        # format dataframe for plotting
        # Categories to string format
        aggdf = aggdf.reset_index()
        for idx_col in aggdf:
            if idx_col == horizontal_axis:
                continue  # if numeric, let it be numeric!
            aggdf[idx_col] = aggdf[idx_col].astype(str)
        aggdf = aggdf.set_index(aggregation)

        # sorting cateigories (months and seisons)

        seasons = ["winter", "spring", "summer", "autumn"]
        months = [
            "January",
            "February",
            "March",
            "April",
            "May",
            "June",
            "July",
            "August",
            "September",
            "October",
            "November",
            "December",
        ]

        season_order_dict = {}
        months_order_dict = {}
        for i, item in enumerate(seasons):
            season_order_dict[item] = i
        for i, item in enumerate(months):
            months_order_dict[item] = i

        # Sort columns
        aggdf = aggdf.reset_index()
        sort_list = aggregation.copy()
        if "season" in aggdf.columns:
            aggdf["season_num"] = aggdf["season"].map(season_order_dict)
            sort_list = ["season_num" if x == "season" else x for x in sort_list]
        if "month" in aggdf.columns:
            aggdf["month_num"] = aggdf["month"].map(months_order_dict)
            sort_list = ["month_num" if x == "month" else x for x in sort_list]
        # sort dataframe
        aggdf = aggdf.sort_values(sort_list, axis=0)
        # drop dummy num coluns (if they are present)
        aggdf = aggdf.drop(columns=["season_num", "month_num"], errors="ignore")
        # reset the index
        aggdf = aggdf.set_index(aggregation)

        # unstack aggregation
        aggregation.remove(horizontal_axis)  # let horizontal axes be the index
        all_stats = aggdf.unstack(aggregation)  # return on verbose

        # Sort index if categorical
        if all_stats.index.name == "season":
            all_stats = all_stats.reindex(seasons)
        if all_stats.index.name == "month":
            all_stats = all_stats.reindex(months)

        # split in values and std
        values_df = all_stats[aggregation_method]
        std_df = all_stats["std"]

        # make shure all data is numeric
        values_df = values_df.astype(float)
        std_df = std_df.astype(float)

        # squize all column levels to one category for plotting
        if len(aggregation) > 1:  # more than one level for the columns
            values_df.columns = [
                " ,".join(col).strip() for col in values_df.columns.values
            ]
            std_df.columns = [" ,".join(col).strip() for col in std_df.columns.values]

        if plot:
            # description of the obstype
            desc_dict = self.data_template[obstype].to_dict()
            if "description" not in desc_dict:
                desc_dict["description"] = obstype

            if not isinstance(desc_dict["description"], str):
                desc_dict["description"] = obstype

            description = desc_dict["description"]

            # generate title
            if title is None:
                startdtstr = datetime.strftime(
                    startdt, format=self.settings.app["print_fmt_datetime"]
                )
                enddtstr = datetime.strftime(
                    enddt, format=self.settings.app["print_fmt_datetime"]
                )
                title = f"{aggregation_method} - {horizontal_axis } {obstype} cycle for period {startdtstr} - {enddtstr} grouped by {aggregation}"

            # ylabel
            if y_label is None:
                if "units" not in desc_dict:
                    y_label = f'{desc_dict["description"]} (units unknown)'
                else:
                    y_label = f'{desc_dict["description"]} ({desc_dict["units"]})'
            else:
                y_label = str(y_label)

            # generate errorbands df
            if errorbands:
                stddf = std_df
            else:
                stddf = None

            # Make plot
            ax = cycle_plot(
                cycledf=values_df,
                errorbandsdf=stddf,
                title=title,
                plot_settings=self.settings.app["plot_settings"]["diurnal"],
                aggregation=aggregation,
                data_template=self.data_template,
                obstype=obstype,
                y_label=y_label,
                legend=legend,
                show_zero_horizontal=_show_zero_line,
            )

            ax.set_ylabel(y_label)
            if horizontal_axis == "hour":
                # extract timezone
                tzstring = str(self.df.index.get_level_values("datetime").tz)

                ax.xaxis.set_major_formatter("{x:.0f} h")
                ax.set_xlabel(f"Hours (timezone: {tzstring})")

        if verbose:
            if plot:
                return values_df, all_stats, ax
            return values_df, all_stats

        return values_df

    # =============================================================================
    # Correlations analysis
    # =============================================================================

    def get_lc_correlation_matrices(self, obstype=["temp"], groupby_labels=["hour"]):
        """Compute pearson correlation coeficients.

        A method to compute the Pearson correlation between an obervation type
        and present landcover fractions in the metadf.

        The correlations are computed per group as defined by unique combinations
        of the groupby_labels.

        A dictionary is returnd where each key represents a unique combination of
        the groupby_labels. The value is a dictionary with the following keys
        and values:

        * cor matrix: the Pearson correlation matrix
        * significance matrix: the significance (p-)values of the correlations.
        * combined matrix: A human readable combination of the correlations and their p values. Indicate by \*, \*\* or \*\*\* representing p-values < 0.05, 0.01 and 0.001 respectively.

        This dictionary is also stored as a lc_cor_dict attribute.

        Parameters
        ----------
        obstype : str, or list optional
            The observation type(s) to compute the correlations on. The default is ['temp'].
        groupby_labels : list, optional
            List of variables to form one group, resulting in one correlation.
            These variables should either a categorical observation type, a categorical column in the metadf or
            a time aggregation. All possible time aggreagetions are: ['minute',
            'hour', 'month', 'year', 'day_of_year',
            'week_of_year', 'season']. The default is ['hour'].

        Returns
        -------
        cor_dict : dict
            A nested dictionary with unique combinations of groupby values.

        """
        if not isinstance(obstype, list):
            obstype = [obstype]

        # get data
        df = self.df[obstype].reset_index()
        df = _make_time_derivatives(df, groupby_labels)

        for group_lab in groupby_labels:
            if group_lab in self.metadf.columns:
                df = df.merge(
                    self.metadf[[group_lab]],
                    how="left",
                    left_on="name",
                    right_index=True,
                )

        for group_lab in groupby_labels:
            assert (
                group_lab in df.columns
            ), f'"{group_lab}" is found in the observations of possible groupby_labels.'

        # subset columns
        relev_columns = [label for label in groupby_labels]  # to avoid deep copy import
        relev_columns.append("name")
        relev_columns.extend(obstype)
        df = df[relev_columns]

        # find landcover columnnames in the metadf
        lc_columns = [
            col for col in self.metadf.columns if (("_" in col) & (col.endswith("m")))
        ]

        # get landcover data
        lc_df = self.metadf[lc_columns]

        if lc_df.empty:
            logger.warning(
                "No landcover columns found in the metadf. Landcover correlations cannot be computed."
            )
            return None

        # merge together
        df = df.merge(lc_df, how="left", left_on="name", right_index=True)

        # remove name column if it is not explicit in the groupby labels
        if "name" not in groupby_labels:
            df = df.drop(columns=["name"])

        # create return
        cor_dict = {}

        # Iterate over all groups

        # avoid futur pandas warning for groupby labels of len==1
        if len(groupby_labels) == 1:
            groups = df.groupby(groupby_labels[0])
        else:
            groups = df.groupby(groupby_labels)

        for group_lab, groupdf in groups:
            # No correlations can be computed when no variance is found
            if groupdf.shape[0] <= 1:
                logger.warning(
                    f"No variance found in correlationd group {group_lab}. Correlation thus not be computed for this group: {groupdf}."
                )
                continue
            # drop groupby labels
            groupdf = groupdf.drop(columns=groupby_labels, errors="ignore")

            rho = groupdf.corr(method="pearson")
            pval = groupdf.corr(method=lambda x, y: pearsonr(x, y)[1]) - np.eye(
                *rho.shape
            )
            # represent p values by stars
            p_stars = pval.applymap(
                lambda x: "".join(["*" for t in [0.05, 0.01, 0.001] if x <= t])
            )

            # combined human readable df
            comb_df = pd.DataFrame(index=rho.index)
            for col in rho.columns:
                comb_df[col] = (
                    rho[col].apply(lambda x: f"{x:.02f}") + " " + p_stars[col]
                )

            cor_dict[group_lab] = {
                "cor matrix": rho,
                "significance matrix": pval,
                "combined matrix": comb_df,
            }

        # Update attribute
        self.lc_cor_dict = cor_dict
        self._lc_cor_obstype = obstype
        self._lc_groupby_labels = groupby_labels

        return cor_dict

    def plot_correlation_heatmap(
        self, groupby_value=None, title=None, _return_ax=False
    ):
        """Make a heatmap plot af a correaltion matrix.

        To specify which correlation matrix to plot, specify the group value
        using the groupby_value argument.

        All possible groupby_values are the keys of the lc_cor_dict attribute.

        Parameters
        ----------
        groupby_value : str, num, None, optional
            A groupby value to indicate which correlation matrix to visualise.
            If None is given, the first groupby value that is present is
            chosen.The default is None.
        title : str, optional
            Title of the figure. If None, a default title is constructed.The
            default is None.

        Returns
        -------
        None.

        Note
        ------
        To list all possible groupby_values, one can use
        ` print(Analysis_instance.lc_cor_dict.keys())`

        """
        # check if there are correlation matrices
        assert bool(
            self.lc_cor_dict
        ), "No correlation matrices found, use the metod get_lc_correlation_matrices first."

        if groupby_value is None:
            groupby_value = list(self.lc_cor_dict.keys())[0]
            logger.warning(
                "No groupby_value is given, so the first groupby value (={groupby_value}) will be used!"
            )
            logger.info(
                f"The correlations are computed over {self._lc_groupby_labels} with the following unique values: {list(self.lc_cor_dict.keys())}"
            )

        # check if groupby value exists
        assert (
            groupby_value in self.lc_cor_dict.keys()
        ), f"{groupby_value} not found as a groupby value. These are all the possible values: {self.lc_cor_dict.keys()}"

        if title is None:
            title = f"Correlation heatmap for group: {self._lc_groupby_labels} = {groupby_value}"

        ax = heatmap_plot(
            cor_dict=self.lc_cor_dict[groupby_value],
            title=title,
            heatmap_settings=self.settings.app["plot_settings"]["correlation_heatmap"],
        )

        if _return_ax:
            return ax

    def plot_correlation_variation(self, title=None):
        """Create correlation scatter plot.

        Make a scatter plot of the correlations to visualise differences between
        multiple group values.

        Group values are represented by the horizontal axes, and correlations
        on the vertical axe.

        All correlations, that are not constant, are plotted as scatters with
        unique colors.

        The scatter marker indicates the p-value of the correlations.

        Parameters
        ----------
        title : str, optional
            Title of the figure. If None, a default title is constructed. The
            default is None.

        Returns
        -------
        None.

        Note
        ------
        If to many possible group values exist, one can use the apply_filter()
        method to reduce the group values.
        """
        # check if there are correlation matrices
        assert bool(
            self.lc_cor_dict
        ), "No correlation matrices found, use the metod get_lc_correlation_matrices first."

        # check if correlation evolution information is available
        if len(self.lc_cor_dict.keys()) <= 1:
            logger.warning(
                f"Only one correlation group is found: {self.lc_cor_dict.keys()}"
            )
            logger.warning("The variance plot can not be made.")
            return

        if title is None:
            title = f"Correlation scatter for group: {self._lc_groupby_labels}"

        ax = correlation_scatter(
            full_cor_dict=self.lc_cor_dict,
            groupby_labels=self._lc_groupby_labels,
            obstypes=self._lc_cor_obstype,
            title=title,
            cor_scatter_settings=self.settings.app["plot_settings"][
                "correlation_scatter"
            ],
        )
        return ax


def _make_time_derivatives(df, required, get_all=False):
    """Construct time derivated columns if required.

    datetime must be a column.
    """
    if ("minute" in required) | (get_all):
        df["minute"] = df["datetime"].dt.minute
    if ("hour" in required) | (get_all):
        df["hour"] = df["datetime"].dt.hour
    if ("month" in required) | (get_all):
        df["month"] = df["datetime"].dt.month_name()
    if ("year" in required) | (get_all):
        df["year"] = df["datetime"].dt.year
    if ("day_of_year" in required) | (get_all):
        df["day_of_year"] = df["datetime"].dt.day_of_year
    if ("week_of_year" in required) | (get_all):
        df["week_of_year"] = df["datetime"].dt.isocalendar()["week"]
    if ("season" in required) | (get_all):
        df["season"] = get_seasons(df["datetime"])

    return df


def get_seasons(
    datetimeseries,
    start_day_spring="01/03",
    start_day_summer="01/06",
    start_day_autumn="01/09",
    start_day_winter="01/12",
):
    """Convert a datetimeseries to a season label (i.g. categorical).

    Parameters
    ----------
    datetimeseries : datetime.datetime
        The timeseries that you want to split up in seasons.
    start_day_spring : str , optional
        Start date for spring, default is '01/03' and if changed the input
        should have the same format as the default value.
    start_day_summer : str , optional
        Start date for summer, default is '01/06' and if changed the input
        should have the same format as the default value.
    start_day_autumn : str , optional
        Start date for autumn, default is '01/09' and if changed the input
        should have the same format as the default value.
    start_day_winter : str , optional
        Start date for winter, default is '01/12' and if changed the input
        should have the same format as the default value.

    Returns
    -------
    output : dataframe
        A obtained dataframe that has where a label for the seasons has been added.
    """
    spring_startday = datetime.strptime(start_day_spring, "%d/%m")
    summer_startday = datetime.strptime(start_day_summer, "%d/%m")
    autumn_startday = datetime.strptime(start_day_autumn, "%d/%m")
    winter_startday = datetime.strptime(start_day_winter, "%d/%m")

    seasons = pd.Series(
        index=["spring", "summer", "autumn", "winter"],
        data=[spring_startday, summer_startday, autumn_startday, winter_startday],
        name="startdt",
    ).to_frame()
    seasons["day_of_year"] = seasons["startdt"].dt.day_of_year - 1

    bins = [0]
    bins.extend(seasons["day_of_year"].to_list())
    bins.append(366)

    labels = ["winter", "spring", "summer", "autumn", "winter"]

    return pd.cut(
        x=datetimeseries.dt.day_of_year,
        bins=bins,
        labels=labels,
        ordered=False,
    )


def filter_data(df, metadf, quarry_str):
    """Filter a dataframe by a user definde string expression.

    This can be used to filter the observation to specific meteorological
    conditions (i.e. low windspeeds, high humidity, cold temperatures, ...)

    The filter expression contains only columns present in the df and/or the
    metadf.

    The filtered df and metadf are returned

    Parameters
    ----------
    df : pandas.DataFrame
        The dataframe containing all the observations to be filterd.
    metadf : pandas.DataFrame
        The dataframe containig all the metadata per station.
    quarry_str : str
        A filter expression using columnnames present in either df or metadf.
        The following timestamp derivatives can be used as well: [minute, hour,
        month, year, day_of_year, week_of_year, season]. The quarry_str may
        contain number and expressions like <, >, ==, >=, \*, +, .... Multiple filters
        can be combine to one expression by using & (AND) and | (OR).

    Returns
    -------
    filter_df : pandas.DataFrame
        The filtered df.
    filter_metadf : pandas.DataFrame
        The filtered metadf.

    """
    # save index order and names for reconstruction
    df_init_idx = list(df.index.names)
    metadf_init_idx = list(metadf.index.names)

    # easyer for sperationg them
    df = df.reset_index()
    metadf = metadf.reset_index()

    # save columns orders
    df_init_cols = df.columns
    metadf_init_cols = metadf.columns

    # create time derivative columns
    df = _make_time_derivatives(df, required=" ", get_all=True)

    # merge together on name
    mergedf = df.merge(metadf, how="left", on="name")

    # apply filter
    filtered = mergedf.query(expr=quarry_str)

    # split to df and metadf
    filter_df = filtered[df_init_cols]
    filter_metadf = filtered[metadf_init_cols]

    # set indexes
    filter_df = filter_df.set_index(df_init_idx)
    filter_metadf = filter_metadf.set_index(metadf_init_idx)

    return filter_df, filter_metadf
