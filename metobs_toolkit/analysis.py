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
from scipy.stats import pearsonr
from metobs_toolkit.settings_files.default_formats_settings import (
    gapfill_label_group,
    label_def,
)
from metobs_toolkit.plotting_functions import (
    cycle_plot,
    heatmap_plot,
    correlation_scatter,
)

from metobs_toolkit.df_helpers import (
    xs_save,
    datetime_subsetting,
    subset_stations,
    empty_outliers_df,
)

logger = logging.getLogger(__name__)

from metobs_toolkit import Dataset


class Analysis(Dataset):
    """The Analysis class contains methods for analyzing observations."""

    def __init__(self, orig_dataset, use_gapfilled_values=False):
        """Create an instance of Analysis.

        The Analysis is similar as a Dataset and contains observations,
        that are considered to be good. There are no outliers.

        The methods of the Analysis class target on scientific analysis, and
        is a feature-futured class.




        Parameters
        ----------
        orig_dataset : metobs_toolkit.Dataset
            The Dataset to create an Analysis of.
        use_gapfilled_values : bool, optional
            If True, the gapfilled records are used in the Analysis,
            else the gapfilled records are (like outliers) present
            as NaN values. The default is False.

        Returns
        -------
        metobs_toolkit.Analysis

        Examples
        --------
        An Analysis is always created from a Dataset, so we start by creating
        a Dataset and importing data into it.

        >>> import metobs_toolkit
        >>>
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> dataset.import_data_from_file(
        ...                         input_data_file=metobs_toolkit.demo_datafile,
        ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                         template_file=metobs_toolkit.demo_template,
        ...                         )

        Now we have a dataset containing records. We now create an Analysis from it.

        >>> ana = metobs_toolkit.Analysis(orig_dataset=dataset)
        >>> ana # doctest: +ELLIPSIS
        Instance of Analysis at ...

        """

        # analysis objects
        self.lc_cor_dict = {}
        self._lc_cor_obstype = None
        self._lc_groupby_labels = None
        self.use_gapfilled_values = use_gapfilled_values

        # check if there is data
        orig_dataset._data_is_required_check()

        # overload attributes
        self._set_df(orig_dataset.df)
        self._set_outliersdf(empty_outliers_df())  # NO OUTLIERS

        # becomes a list of gaps required for including gapfilled records
        self._set_gaps(orig_dataset.gaps)

        # Dataset with metadata (static)
        self._set_metadf(orig_dataset.metadf)
        # dictionary storing present observationtypes
        self._set_obstypes(orig_dataset.obstypes)  # init with all tlk obstypes
        self._istype = "Analysis"

        self._set_settings(orig_dataset.settings)

        # add empty lcz column to metadf if it is not present
        if "lcz" not in self.metadf.columns:
            self.metadf["lcz"] = np.nan

    def get_analysis_records(self):
        """
        Get the Analysis records in a DataFrame form.

        The Analysis records are the observations that are used by the
        Analysis methods.

        Returns
        -------
        widedf : pandas.DataFrame
            The observation records are presented in a wide structure where each
            column represents an observationtype.

        Examples
        --------
        An Analysis is always created from a Dataset, so we start by creating
        a Dataset and importing data into it.

        >>> import metobs_toolkit
        >>>
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> dataset.import_data_from_file(
        ...                         input_data_file=metobs_toolkit.demo_datafile,
        ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                         template_file=metobs_toolkit.demo_template,
        ...                         )

        Now we have a dataset containing records.We now create an Analysis from it.

        >>> ana = metobs_toolkit.Analysis(orig_dataset=dataset)
        >>> ana # doctest: +ELLIPSIS
        Instance of Analysis at ...

        To get all the records in the Analysis in the form of a pandas.DataFrame,
        we use `get_analysis_records()`.

        >>> df = ana.get_analysis_records()
        >>> df
        obstype                              humidity  temp  wind_direction  wind_speed
        name      datetime
        vlinder01 2022-09-01 00:00:00+00:00      65.0  18.8            65.0        1.56
                  2022-09-01 00:05:00+00:00      65.0  18.8            75.0        1.53
                  2022-09-01 00:10:00+00:00      65.0  18.8            75.0        1.42
                  2022-09-01 00:15:00+00:00      65.0  18.7            85.0        1.67
                  2022-09-01 00:20:00+00:00      65.0  18.7            65.0        1.39
        ...                                       ...   ...             ...         ...
        vlinder28 2022-09-15 23:35:00+00:00      77.0  13.4           275.0        0.00
                  2022-09-15 23:40:00+00:00      77.0  13.3           275.0        0.00
                  2022-09-15 23:45:00+00:00      77.0  13.2           275.0        0.00
                  2022-09-15 23:50:00+00:00      77.0  13.2           275.0        0.00
                  2022-09-15 23:55:00+00:00      77.0  13.0           285.0        0.00
        <BLANKLINE>
        [120960 rows x 4 columns]
        """

        if self.use_gapfilled_values:
            mergedf = self.get_full_status_df(return_as_wide=False)

            # get all ok and all fill-method labels
            fill_labels = [label_def[method]["label"] for method in gapfill_label_group]
            wanted_labels = [label_def["goodrecord"]["label"]]
            wanted_labels.extend(fill_labels)

            # subset to good and gapfilled labels
            df = mergedf[mergedf["label"].isin(wanted_labels)]

        else:
            df = self.df

        # TO WIDE
        widedf = df.unstack("obstype")["value"]
        return widedf

    # =============================================================================
    #     Setters
    # =============================================================================

    def subset_period(self, startdt, enddt):
        """Subset the observations of the Analysis to a specific period.


        Parameters
        ----------
        startdt : datetime.datetime
            The start datetime to filter the observations.
        enddt : datetime.datetime
            The end datetime to filter the observations.

        Returns
        -------
        Analysis()
            A new Analysis is returned that is the subset version of the one
            the filter is applied on.

        See Also
        --------
        Analysis: The Analysis class
        Analysis.apply_filter: Apply a filter on the analysis.
        Analysis.get_possible_filter_keywords: Get all possible filter keywords.
        Analysis.get_analysis_records : Get the records in a pandas DataFrame.

        Note
        --------
        If a timezone unaware datetime is given as an argument, it is interpreted
        as if it has the same timezone as the observations.

        Examples
        --------
        An Analysis is always created from a Dataset, so we start by creating
        a Dataset and importing data into it.

        >>> import metobs_toolkit
        >>>
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> dataset.import_data_from_file(
        ...                         input_data_file=metobs_toolkit.demo_datafile,
        ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                         template_file=metobs_toolkit.demo_template,
        ...                         )

        Now we have a dataset containing records. We now create an Analysis from it.

        >>> ana = metobs_toolkit.Analysis(orig_dataset=dataset)
        >>> ana # doctest: +ELLIPSIS
        Instance of Analysis at ...

        To subset the Analysis in time we use the `subset_period()` method.

        >>> import datetime
        >>>
        >>> small_ana = ana.subset_period(startdt = datetime.datetime(2022,9,5),
        ...                               enddt = datetime.datetime(2022,9,7))
        >>> small_ana.df
                                                        value
        name      obstype    datetime
        vlinder01 humidity   2022-09-05 00:00:00+00:00  80.00
                             2022-09-05 00:05:00+00:00  78.00
                             2022-09-05 00:10:00+00:00  77.00
                             2022-09-05 00:15:00+00:00  76.00
                             2022-09-05 00:20:00+00:00  75.00
        ...                                               ...
        vlinder28 wind_speed 2022-09-06 23:40:00+00:00   0.00
                             2022-09-06 23:45:00+00:00   0.00
                             2022-09-06 23:50:00+00:00   0.08
                             2022-09-06 23:55:00+00:00   1.00
                             2022-09-07 00:00:00+00:00   1.00
        <BLANKLINE>
        [64624 rows x 1 columns]

        """
        logger.info(f"Subsetting {self.__repr__()} to ({startdt}) -- ({enddt}) period.")

        # arg format checks
        startdt = self._datetime_arg_check(startdt)
        enddt = self._datetime_arg_check(enddt)

        df = datetime_subsetting(self.df, startdt, enddt)
        if df.empty:
            logger.warning("The period subset results in an empty Analysis!")
        # create new analysis
        Filtered_Dataset = Dataset()
        Filtered_Dataset._set_df(df)
        Filtered_Dataset._set_metadf(self.metadf)
        Filtered_Dataset._set_gaps(self.gaps)
        Filtered_Dataset._set_settings(self.settings)
        Filtered_Dataset._set_obstypes(self.obstypes)

        subsetted_an = Analysis(
            orig_dataset=Filtered_Dataset,
            use_gapfilled_values=self.use_gapfilled_values,
        )
        return subsetted_an

    # =============================================================================
    #   Helpers
    # =============================================================================

    def get_possible_filter_keywords(self):
        """Get all keywords that can be used for filtering.


        Returns
        -------
        list
            A list of string keywords that can be used in the
            `metobs_toolkit.Analysis.apply_filter()` method.

        See Also
        --------
        Analysis: The Analysis class
        Analysis.apply_filter: Apply a filter on the analysis.


        Examples
        --------
        An Analysis is always created from a Dataset, so we start by creating
        a Dataset and importing data into it.

        >>> import metobs_toolkit
        >>>
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> dataset.import_data_from_file(
        ...                         input_data_file=metobs_toolkit.demo_datafile,
        ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                         template_file=metobs_toolkit.demo_template,
        ...                         )

        Now we have a dataset containing records.We now create an Analysis from it.

        >>> ana = metobs_toolkit.Analysis(orig_dataset=dataset)

        To get all possible keywords for filtering the analysis we use
        `get_possible_filter_keywords()`

        >>> ana.get_possible_filter_keywords()
        ['dataset_resolution', 'datetime', 'day_of_year', 'dt_end', 'dt_start', 'geometry', 'hour', 'humidity', 'lat', 'lon', 'minute', 'month', 'name', 'school', 'season', 'temp', 'week_of_year', 'wind_direction', 'wind_speed', 'year']

        If you want to filter to, or aggregate to, LCZ for example then you need
        the LCZ information in the metadata. If you do not have the LCZ's per
        station, you can use the `Dataset.get_lcz()` method to extract them. So
        make sure that you get the LCZ's info in the Dataset, and then create
        an Analysis from it.

        >>> #Create your Dataset and add LCZ to it
        >>> dataset_with_lcz = metobs_toolkit.Dataset() #empty Dataset
        >>> dataset_with_lcz .import_data_from_file(
        ...                         input_data_file=metobs_toolkit.demo_datafile,
        ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                         template_file=metobs_toolkit.demo_template,
        ...                         )
        >>> # extract LCZ
        >>> lcz_series = dataset_with_lcz .get_lcz()
        >>>
        >>> # create an analysis from it
        >>> ana_with_lcz = metobs_toolkit.Analysis(dataset_with_lcz)

        Now the LCZ information is transferred to the ana_with_lcz.

        >>> ana_with_lcz.get_possible_filter_keywords()
        ['dataset_resolution', 'datetime', 'day_of_year', 'dt_end', 'dt_start', 'geometry', 'hour', 'humidity', 'lat', 'lcz', 'lon', 'minute', 'month', 'name', 'school', 'season', 'temp', 'week_of_year', 'wind_direction', 'wind_speed', 'year']

        >>> 'lcz' in ana_with_lcz.get_possible_filter_keywords()
        True

        """
        filter_keys = []
        # all present opbstypes
        filter_keys.extend(self.df.index.get_level_values("obstype").unique().to_list())

        # add 'name'
        filter_keys.append("name")

        # all metadata coluns that is not all Nan's
        for metacol in self.metadf.columns:
            if not self.metadf[metacol].isna().all():
                filter_keys.append(metacol)

        # all time derivatives
        dummydf = self.df[:5].reset_index()[["datetime"]]
        time_agg_cols = _make_time_derivatives(df=dummydf, required=[], get_all=True)
        filter_keys.extend(time_agg_cols)

        filter_keys = list(set(filter_keys))
        filter_keys.sort()
        return filter_keys

    def apply_filter(self, expression):
        """Filter an Analysis by a user defined string expression.

        This can be used to filter the observation to specific meteorological
        conditions (i.e. low windspeeds, high humidity, cold temperatures, ...)

        The filtered expression contains references to present observationtypes,
        metadata and timestamp aggregates. Use the
        `metobs_toolkit.Analysis.get_possible_filter_keywords()` method for a
        list of all possible keywords.

        Records that are not compliant with the filter expression are set to Nan,
        and are thus ignored in all other methods.


        Parameters
        ----------
        expression : str
            A filter expression using columnnames present in either df or metadf.
            The following timestamp derivatives can be used as well: [minute, hour,
            month, year, day_of_year, week_of_year, season]. The quarry_str may
            contain number and expressions like <, >, ==, >=, \*, +, .... Multiple filters
            can be combined to one expression by using & (AND) and | (OR).

        Returns
        -------
        Analysis()
            A new Analysis is returned that is the filtered version of the one
            the filter is applied on.

        See Also
        --------
        Analysis: The Analysis class.
        Analysis.get_possible_filter_keywords: Get all possible keywords for filtering.
        Analysis.subset_period: Subset the data to a time period.

        Note
        ------
        Make sure to use \" of \' to indicate string values in the expression if
        needed.

        Examples
        --------

        .. plot::
            :context: close-figs

            An Analysis is always created from a Dataset, so we start by creating
            a Dataset and importing data into it.

            >>> import metobs_toolkit
            >>>
            >>> #Create your Dataset
            >>> dataset = metobs_toolkit.Dataset() #empty Dataset
            >>> dataset.import_data_from_file(
            ...                         input_data_file=metobs_toolkit.demo_datafile,
            ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
            ...                         template_file=metobs_toolkit.demo_template,
            ...                         )

            Now we have a dataset containing records.We now create an Analysis from it.

            >>> ana = metobs_toolkit.Analysis(orig_dataset=dataset)
            >>> ana # doctest: +ELLIPSIS
            Instance of Analysis at ...

            To get all possible filter keywords we use the `get_possible_filter_keywords()`
            method. We can also inspect the records in a pandas.DataFrame format
            to look for thresholds for the filter.

            >>> ana.get_possible_filter_keywords()
            ['dataset_resolution', 'datetime', 'day_of_year', 'dt_end', 'dt_start', 'geometry', 'hour', 'humidity', 'lat', 'lon', 'minute', 'month', 'name', 'school', 'season', 'temp', 'week_of_year', 'wind_direction', 'wind_speed', 'year']

            >>> df = ana.get_analysis_records()
            >>> df
            obstype                              humidity  temp  wind_direction  wind_speed
            name      datetime
            vlinder01 2022-09-01 00:00:00+00:00      65.0  18.8            65.0        1.56
                      2022-09-01 00:05:00+00:00      65.0  18.8            75.0        1.53
                      2022-09-01 00:10:00+00:00      65.0  18.8            75.0        1.42
                      2022-09-01 00:15:00+00:00      65.0  18.7            85.0        1.67
                      2022-09-01 00:20:00+00:00      65.0  18.7            65.0        1.39
            ...                                       ...   ...             ...         ...
            vlinder28 2022-09-15 23:35:00+00:00      77.0  13.4           275.0        0.00
                      2022-09-15 23:40:00+00:00      77.0  13.3           275.0        0.00
                      2022-09-15 23:45:00+00:00      77.0  13.2           275.0        0.00
                      2022-09-15 23:50:00+00:00      77.0  13.2           275.0        0.00
                      2022-09-15 23:55:00+00:00      77.0  13.0           285.0        0.00
            <BLANKLINE>
            [120960 rows x 4 columns]

            We create a filter and apply it

            >>> filtered_ana = ana.apply_filter(
            ... expression="temp <= 16.3 & wind_direction > 180 & season=='autumn'")
            >>>
            >>> filtered_ana.get_analysis_records()
            obstype                              humidity  temp  wind_direction  wind_speed
            name      datetime
            vlinder01 2022-09-01 00:00:00+00:00       NaN   NaN             NaN         NaN
                      2022-09-01 00:05:00+00:00       NaN   NaN             NaN         NaN
                      2022-09-01 00:10:00+00:00       NaN   NaN             NaN         NaN
                      2022-09-01 00:15:00+00:00       NaN   NaN             NaN         NaN
                      2022-09-01 00:20:00+00:00       NaN   NaN             NaN         NaN
            ...                                       ...   ...             ...         ...
            vlinder28 2022-09-15 23:35:00+00:00      77.0  13.4           275.0         0.0
                      2022-09-15 23:40:00+00:00      77.0  13.3           275.0         0.0
                      2022-09-15 23:45:00+00:00      77.0  13.2           275.0         0.0
                      2022-09-15 23:50:00+00:00      77.0  13.2           275.0         0.0
                      2022-09-15 23:55:00+00:00      77.0  13.0           285.0         0.0
            <BLANKLINE>
            [120960 rows x 4 columns]

            All records that do not meet the filter are converted to NaN. These values
            will be ignored in all methods of the `Analysis` class.

            >>> filtered_ana.make_plot(obstype='temp')
            <Axes: ylabel='temp (Celsius)'>

        """
        logger.info(f"applying filter ({expression}) on {self.__repr__()}")
        # apply filter on wide df
        _initial_n_records = self.df["value"].count()  # for logging
        dfwide = self.get_analysis_records()
        present_obstypes = list(dfwide.columns)
        dfwide = dfwide.reset_index()

        # merge metadata
        dfwide = dfwide.merge(self.metadf, how="left", left_on="name", right_index=True)

        # merge daterelated
        dfwide = _make_time_derivatives(df=dfwide, required=[], get_all=True)

        # apply querry
        filtered = dfwide.query(expr=expression)

        # To triple index
        relev_columns = present_obstypes
        relev_columns.extend(["name", "datetime"])
        longfiltered = (
            filtered[relev_columns]
            .set_index(["name", "datetime"])
            .stack(future_stack=True)
        )
        longfiltered.index.rename("obstype", level=-1, inplace=True)
        longfiltered.name = "value"
        longfiltered = longfiltered.reset_index().set_index(
            ["name", "obstype", "datetime"]
        )
        # longfiltered = longfiltered.to_frame()

        _filtered_n_records = longfiltered["value"].count()  # for logging

        logger.info(
            f"Filtering resulted in a loss of {(1-(_filtered_n_records/_initial_n_records))*100:.2f}% of data."
        )
        if longfiltered.empty:
            logger.warning(f"No data is left after applying the {expression} filter!")

        # create new analysis
        Filtered_Dataset = Dataset()
        Filtered_Dataset._set_df(self.df.copy(deep=True))

        # (do not drop idexes, because then are no ideal timestamps, and the link
        # to dataset is broken)
        Filtered_Dataset.df.loc[
            Filtered_Dataset.df.index.difference(longfiltered.index), :
        ] = np.nan

        Filtered_Dataset._set_metadf(self.metadf)
        Filtered_Dataset._set_gaps(self.gaps)
        Filtered_Dataset._set_settings(self.settings)
        Filtered_Dataset._set_obstypes(self.obstypes)

        filtered_an = Analysis(
            orig_dataset=Filtered_Dataset,
            use_gapfilled_values=self.use_gapfilled_values,
        )

        return filtered_an

    def aggregate_df(self, df=None, agg=["lcz", "hour"], method="mean"):
        """Aggregate observations to a (list of) categories.

        The output will be a dataframe that is aggregated to one, or more
        categories. A common example is aggregating to LCZ's.


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

        See Also
        --------
        Analysis: The Analysis class.
        Analysis.get_diurnal_statistics: Aggregate to diurnal cycle.
        Analysis.get_annual_statistics: Aggregate to annual cycle.

        Note
        -------
        Present columns that are non-numeric and are not in the agg list, are
        not present in the return, since these values cannot be aggregated.

        Examples
        --------
        An Analysis is always created from a Dataset, so we start by creating
        a Dataset and importing data into it.

        >>> import metobs_toolkit
        >>>
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> dataset.import_data_from_file(
        ...                         input_data_file=metobs_toolkit.demo_datafile,
        ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                         template_file=metobs_toolkit.demo_template,
        ...                         )
        >>> # Get LCZ (so we can aggregate per LCZ)
        >>> lcz_series = dataset.get_lcz()

        Now we have a dataset containing records.We now create an Analysis from it.

        >>> ana = metobs_toolkit.Analysis(orig_dataset=dataset)
        >>> ana # doctest: +ELLIPSIS
        Instance of Analysis at ...

        As an example, we aggregate to LCZ and the day of the year.

        >>> aggdf = ana.aggregate_df(agg=['lcz', 'day_of_year'],
        ...                          method='median')
        >>> aggdf
                                     humidity   temp  wind_direction  wind_speed
        lcz             day_of_year
        Compact midrise 244              51.0  21.05            85.0        0.31
                        245              43.0  26.50            85.0        0.31
                        246              66.0  21.70            55.0        0.08
                        247              58.0  22.25            45.0        0.11
                        248              63.0  20.60            75.0        0.28
        ...                               ...    ...             ...         ...
        Water (LCZ G)   254              92.0  16.80           175.0        0.17
                        255              87.0  16.80           175.0        0.24
                        256              84.0  19.50            65.0        0.35
                        257              92.0  14.70            65.0        0.75
                        258              89.0  16.70           265.0        1.50
        <BLANKLINE>
        [150 rows x 4 columns]
        """
        if df is None:
            # df = copy.deepcopy(self.df)
            df = self.get_analysis_records()

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
        agg_df = df.groupby(agg, observed=True).agg(
            method, numeric_only=True
        )  # descrepation warning
        # sort index
        agg_df = agg_df.reset_index()
        agg_df = agg_df.set_index(agg)
        return agg_df

    # =============================================================================
    #   Analyse method
    # =============================================================================
    def get_annual_statistics(
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
        legend=True,
        _return_all_stats=False,
    ):
        """
        Create an annual cycle for aggregated groups.

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
            If True, an annual plot is made. The default is True.
        errorbands : bool, optional
             If True, the std is representd in the plot by colored bands. The default is False.
        title : string, optional
             Title of the figure, if None a default title is generated. The default is None.
        legend : bool, optional
             I True, a legend is added to the plot. The default is True.

        Returns
        -------
        df : pandas.DataFrame()
            The dataframe contains the aggregated values.

        See Also
        --------

        Analysis: The Analysis class.
        Analysis.aggregate_df: Data method for aggregating Analysis data.
        Analysis.get_diurnal_statistics: Aggregate to diurnal cycle.
        Analysis.get_annual_statistics: Aggregate to annual cycle.

        Note
        --------

        If a timezone unaware datetime is given as an argument, it is interpreted
        as if it has the same timezone as the observations.

        Examples
        --------

        .. plot::
            :context: close-figs

            An Analysis is always created from a Dataset, so we start by creating
            a Dataset and importing data into it.

            >>> import metobs_toolkit
            >>>
            >>> #Create your Dataset
            >>> dataset = metobs_toolkit.Dataset() #empty Dataset
            >>> dataset.import_data_from_file(
            ...                         input_data_file=metobs_toolkit.demo_datafile,
            ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
            ...                         template_file=metobs_toolkit.demo_template,
            ...                         )
            >>> # Get LCZ (so we can aggregate per LCZ)
            >>> lcz_series = dataset.get_lcz()

            Now we have a dataset containing records.We now create an Analysis from it.

            >>> ana = metobs_toolkit.Analysis(orig_dataset=dataset)
            >>> ana # doctest: +ELLIPSIS
            Instance of Analysis at ...

            We can create a annual cycle, and aggregate to the LCZ.

            >>> annual_aggdf = ana.get_annual_statistics(
            ...                     groupby=["lcz"],
            ...                     obstype="temp",
            ...                     agg_method="mean",
            ...                     plot=True,
            ...                     errorbands=False)
            >>> annual_aggdf[7:]
            lcz        Compact midrise  Dense Trees (LCZ A)  Large lowrise  Low plants (LCZ D)  ...  Open midrise  Scattered Trees (LCZ B)  Sparsely built  Water (LCZ G)
            month                                                                               ...
            August                 NaN                  NaN            NaN                 NaN  ...           NaN                      NaN             NaN            NaN
            September            19.32                17.52          18.76               18.21  ...         18.63                    17.99           18.21           19.1
            October                NaN                  NaN            NaN                 NaN  ...           NaN                      NaN             NaN            NaN
            November               NaN                  NaN            NaN                 NaN  ...           NaN                      NaN             NaN            NaN
            December               NaN                  NaN            NaN                 NaN  ...           NaN                      NaN             NaN            NaN
            <BLANKLINE>
            [5 rows x 10 columns]

        """
        # arg format checks
        startdt = self._datetime_arg_check(startdt)
        enddt = self._datetime_arg_check(enddt)

        assert (
            obstype in self.obstypes.keys()
        ), f"{obstype} is not a known obstype of the Analysis: {self.obstypes}"
        # convert obstype (string) to obstype(Obstype)
        obstype = self.obstypes[obstype]

        if title is None:
            title = f"Anual {obstype.name} ({obstype.get_orig_name()}) cycle plot per {groupby}."
        else:
            title = str(title)

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
            y_label=obstype._get_plot_y_label(),
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
            If True, a diurnal plot is made. The default is True.
        title : string, optional
             Title of the figure, if None a default title is generated. The default is None.
        legend : bool, optional
             I True, a legend is added to the plot. The default is True.
        errorbands : bool, optional
            If True, the std is represented in the plot by colored bands. The default is False.

        Returns
        -------
        df : pandas.DataFrame()
            The dataframe containing the aggregated values.

        See Also
        --------

        Analysis: The Analysis class
        Analysis.aggregate_df: Data method for aggregating Analysis data.
        Analysis.get_diurnal_statistics_with_reference: Diurnal difference to reference cycle.
        Analysis.get_annual_statistics: Aggregate to annual cycle.

        Note
        --------
        If a timezone unaware datetime is given as an argument, it is interpreted
        as if it has the same timezone as the observations.

        Examples
        --------

        .. plot::
            :context: close-figs

            An Analysis is always created from a Dataset, so we start by creating
            a Dataset and importing data into it.

            >>> import metobs_toolkit
            >>>
            >>> #Create your Dataset
            >>> dataset = metobs_toolkit.Dataset() #empty Dataset
            >>> dataset.import_data_from_file(
            ...                         input_data_file=metobs_toolkit.demo_datafile,
            ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
            ...                         template_file=metobs_toolkit.demo_template,
            ...                         )
            >>> # Get LCZ (so we can aggregate per LCZ)
            >>> lcz_series = dataset.get_lcz()

            Now we have a dataset containing records.We now create an Analysis from it.

            >>> ana = metobs_toolkit.Analysis(orig_dataset=dataset)
            >>> ana
            Instance of Analysis at ...

            We can create a diurnal cycle for each station

            >>> diurnal_aggdf = ana.get_diurnal_statistics(
            ...                      colorby="name",
            ...                      obstype="temp",
            ...                      plot=True)

        .. plot::
            :context: close-figs

            >>> diurnal_aggdf
            name  vlinder01  vlinder02  vlinder03  vlinder04  ...  vlinder25  vlinder26  vlinder27  vlinder28
            hour                                              ...
            0         16.19      16.95      16.35      14.81  ...      17.61      16.84      17.84      15.22
            1         15.59      16.45      15.76      14.34  ...      17.14      16.26      17.28      14.74
            2         15.09      15.90      15.13      13.47  ...      16.62      15.61      16.76      14.08
            3         15.21      15.80      14.87      13.32  ...      16.38      15.42      16.57      13.87
            4         15.24      15.64      14.72      13.28  ...      16.09      15.22      16.34      13.81
            ...         ...        ...        ...        ...  ...        ...        ...        ...        ...
            19        17.90      18.77      19.02      17.44  ...      19.05      19.56      19.69      17.55
            20        17.08      18.07      18.23      16.55  ...      18.53      18.65      19.01      16.61
            21        16.65      17.55      17.55      15.95  ...      18.18      18.04      18.52      16.00
            22        16.39      17.11      16.96      15.30  ...      17.80      17.51      18.12      15.56
            23        15.99      16.83      16.42      14.77  ...      17.51      16.95      17.81      15.22
            <BLANKLINE>
            [24 rows x 28 columns]

            We can also aggregate them in groups. As an example, we group them by LCZ.
            >>> diurnal_lcz_aggdf = ana.get_diurnal_statistics(
            ...                      colorby="lcz",
            ...                      obstype="temp",
            ...                      plot=True)

        .. plot::
            :context: close-figs

            >>> diurnal_lcz_aggdf
            lcz   Compact midrise  Dense Trees (LCZ A)  Large lowrise  Low plants (LCZ D)  ...  Open midrise  Scattered Trees (LCZ B)  Sparsely built  Water (LCZ G)
            hour                                                                           ...
            0               17.74                16.14          16.95               15.97  ...         16.60                    15.58           16.23          18.22
            1               17.23                15.57          16.45               15.54  ...         16.01                    15.11           15.81          17.88
            2               16.66                15.09          15.90               14.93  ...         15.37                    14.30           15.20          17.57
            3               16.42                14.91          15.80               14.82  ...         15.14                    14.16           15.00          17.39
            4               16.21                14.77          15.64               14.66  ...         14.97                    13.97           14.87          17.21
            ...               ...                  ...            ...                 ...  ...           ...                      ...             ...            ...
            19              19.69                17.29          18.77               17.80  ...         19.29                    17.63           18.02          19.38
            20              19.10                16.74          18.07               17.06  ...         18.44                    17.04           17.37          18.93
            21              18.54                16.46          17.55               16.60  ...         17.79                    16.52           16.97          18.55
            22              18.03                16.15          17.11               16.27  ...         17.24                    15.96           16.52          18.32
            23              17.70                15.89          16.83               15.88  ...         16.68                    15.54           16.16          18.13
            <BLANKLINE>
            [24 rows x 10 columns]
        """
        # arg format checks
        startdt = self._datetime_arg_check(startdt)
        enddt = self._datetime_arg_check(enddt)

        assert (
            obstype in self.obstypes.keys()
        ), f"{obstype} is not a known obstype of the Analysis: {self.obstypes}"
        # convert obstype (string) to obstype(Obstype)
        obstype = self.obstypes[obstype]

        if title is None:
            if startdt is None:
                if enddt is None:
                    title = f"Hourly average {obstype.name} ({obstype.get_orig_name()} diurnal cycle"
                else:
                    title = f"Hourly average {obstype.name} ({obstype.get_orig_name()}) diurnal cycle until {enddt}"
            else:
                if enddt is None:
                    title = f"Hourly average {obstype.name} ({obstype.get_orig_name()}) diurnal cycle from {startdt}"
                else:
                    title = f"Hourly average {obstype.name} ({obstype.get_orig_name()}) diurnal cycle for period {startdt} - {enddt}"

        else:
            title = str(title)

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
            y_label=obstype._get_plot_y_label(),
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
        tolerance="30min",
        stations=None,
        startdt=None,
        enddt=None,
        plot=True,
        title=None,
        legend=True,
        errorbands=False,
        show_zero_horizontal=True,
        _return_all_stats=False,
    ):
        """
        Create an average diurnal cycle for the observation differences of a reference station.

        All observational values are converted to differences with the closest
        (in time) reference observation. No reference observation is found when
        the time difference is larger than the tolerance.

        (In the plot, each station is represented by a line.)

        Parameters
        ----------
        refstation : str,
            Name of the station to use as a reference.
        colorby : 'name' or 'lcz', optional
            If 'name' the plotted lines will be colored per station, if 'lcz' the colors represent the stations lcz. The default is 'name'.
        obstype : str, optional
            Element of the metobs_toolkit.observation_types The default is 'temp'.
        tolerance : Timedelta or str, optional
            The tolerance string or object represents the maximum translation in time to find a reference
            observation for each observation. Ex: '5min' is 5 minutes, '1h', is one hour. The default is '30min'.
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
        legend : bool, optional
              I True, a legend is added to the plot. The default is True.
        errorbands : bool, optional
            If True, the std is represented in the plot by colored bands. The upper bound represents +1 x std, the lower bound -1 x std. The default is False.
        show_zero_horizontal : bool, optional
            If True a horizontal line is drawn in the plot at zero. The default is True.

        Returns
        -------
        df : pandas.DataFrame()
            The dataframe containing the aggregated values.

        See Also
        --------

        Analysis: The Analysis class
        Analysis.aggregate_df: Data method for aggregating Analysis data.
        Analysis.get_diurnal_statistics: Aggregate to diurnal cycle.
        Analysis.get_annual_statistics: Aggregate to annual cycle.

        Note
        --------
        If a timezone unaware datetime is given as an argument, it is interpreted
        as if it has the same timezone as the observations.

        Examples
        --------

        .. plot::
            :context: close-figs

            An Analysis is always created from a Dataset, so we start by creating
            a Dataset and importing data into it.

            >>> import metobs_toolkit
            >>>
            >>> #Create your Dataset
            >>> dataset = metobs_toolkit.Dataset() #empty Dataset
            >>> dataset.import_data_from_file(
            ...                         input_data_file=metobs_toolkit.demo_datafile,
            ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
            ...                         template_file=metobs_toolkit.demo_template,
            ...                         )

            Now we have a dataset containing records. We now create an Analysis from it.

            >>> ana = metobs_toolkit.Analysis(orig_dataset=dataset)
            >>> ana
            Instance of Analysis at ...

            We can create a diurnal cycle for each station but rather than displaying
            the average observation value, we present the instantanious difference
            with a reference station. By doing so, we can subtract 'the weather' from
            'the climate'.

            >>> diurnal_diff_aggdf = ana.get_diurnal_statistics_with_reference(
            ...                      refstation ="vlinder02",
            ...                      colorby="name",
            ...                      obstype="temp",
            ...                      tolerance="20min",
            ...                      plot=True)

            >>> diurnal_diff_aggdf
            name  vlinder01  vlinder03  vlinder04  vlinder05  ...  vlinder25  vlinder26  vlinder27  vlinder28
            hour                                              ...
            0         -0.76  -5.97e-01      -2.14       1.70  ...       0.66      -0.11       0.89      -1.73
            1         -0.87  -6.89e-01      -2.12       2.20  ...       0.69      -0.19       0.83      -1.71
            2         -0.81  -7.71e-01      -2.43       2.75  ...       0.72      -0.29       0.86      -1.82
            3         -0.59  -9.28e-01      -2.47       2.85  ...       0.58      -0.38       0.77      -1.92
            4         -0.40  -9.16e-01      -2.36       3.01  ...       0.45      -0.42       0.70      -1.83
            ...         ...        ...        ...        ...  ...        ...        ...        ...        ...
            19        -0.87   2.56e-01      -1.33       0.12  ...       0.28       0.79       0.93      -1.22
            20        -0.99   1.67e-01      -1.52       0.55  ...       0.46       0.59       0.94      -1.46
            21        -0.90   3.89e-03      -1.60       0.86  ...       0.63       0.49       0.97      -1.54
            22        -0.72  -1.52e-01      -1.81       1.29  ...       0.68       0.40       1.01      -1.56
            23        -0.84  -4.12e-01      -2.06       1.58  ...       0.68       0.12       0.98      -1.61
            <BLANKLINE>
            [24 rows x 27 columns]

        """
        # arg format checks
        startdt = self._datetime_arg_check(startdt)
        enddt = self._datetime_arg_check(enddt)

        assert (
            obstype in self.obstypes.keys()
        ), f"{obstype} is not a known obstype of the Analysis: {self.obstypes}"
        # convert obstype (string) to obstype(Obstype)
        obstype = self.obstypes[obstype]

        # filter df
        obsdf = xs_save(self.df, obstype.name, "obstype").reset_index()

        # extract refernce from observations
        refdf = obsdf[obsdf["name"] == refstation]
        obsdf = obsdf[obsdf["name"] != refstation]

        assert (
            not refdf.empty
        ), f"Error: No reference observation found (after filtering) for {refstation}"
        assert not obsdf.empty, "Error: No observation found (after filtering)"

        # Syncronize observations with the reference observations
        refdf = refdf.rename(columns={"value": "ref_value", "datetime": "ref_datetime"})

        mergedf = pd.merge_asof(
            left=obsdf.sort_values("datetime"),
            right=refdf[["ref_datetime", "ref_value"]].sort_values("ref_datetime"),
            right_on="ref_datetime",
            left_on="datetime",
            direction="nearest",
            tolerance=pd.Timedelta(tolerance),
        )

        # Get differnces (overwrite the values)
        mergedf["value"] = mergedf["value"] - mergedf["ref_value"]

        # Subset to relavent columns and reconstruct triple index
        mergedf = mergedf.reset_index()
        mergedf["obstype"] = obstype.name
        mergedf = mergedf[["name", "obstype", "datetime", "value"]]
        mergedf = mergedf.set_index(["name", "obstype", "datetime"]).sort_index()

        # title
        if title is None:
            if startdt is None:
                if enddt is None:
                    title = f"Hourly average {obstype.name} diurnal cycle, with {refstation} as reference,"
                else:
                    title = f"Hourly average {obstype.name} diurnal cycle, with {refstation} as reference, until {enddt}"
            else:
                if enddt is None:
                    title = f"Hourly average {obstype.name} diurnal cycle, with {refstation} as reference, from {startdt}"
                else:
                    title = f"Hourly average {obstype.name} diurnal cycle, with {refstation} as reference, for period {startdt} - {enddt}"

        else:
            title = str(title)

        y_label = (
            f"{obstype.name} difference to {refstation} ({obstype.get_standard_unit()})"
        )

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
        """Create an average cycle for an aggregated category.

        A commen example is to aggregate to the LCZ's, so to get the diurnal
        cycle per LCZ rather than per station.

        (In the plot, each aggregated category different from datetime, is represented by a line.)

        Parameters
        ----------
        obstype : str or metobs_toolkit.Obstype, optional
            The observation type to plot. The default is 'temp'.
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
            If True, the std is represented in the plot by colored bands. The default is False.
        verbose : True, optional
            If True, an additional dataframe with aggregation information is returned . The default is False.

        Returns
        -------
        df : pandas.DataFrame()
            The dataframe containing the aggregated values.

        See Also
        --------

        Analysis: The Analysis class
        Analysis.aggregate_df: Data method for aggregating Analysis data.
        Analysis.get_diurnal_statistics: Aggregate to diurnal cycle.
        Analysis.get_diurnal_statistics_with_reference: Diurnal difference to reference cycle.
        Analysis.get_annual_statistics: Aggregate to annual cycle.

        Note
        -------
        If a timezone unaware datetime is given as an argument, it is interpreted
        as if it has the same timezone as the observations.

        Examples
        --------

        .. plot::
            :context: close-figs

            An Analysis is always created from a Dataset, so we start by creating
            a Dataset and importing data into it.

            >>> import metobs_toolkit
            >>>
            >>> #Create your Dataset
            >>> dataset = metobs_toolkit.Dataset() #empty Dataset
            >>> dataset.import_data_from_file(
            ...                         input_data_file=metobs_toolkit.demo_datafile,
            ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
            ...                         template_file=metobs_toolkit.demo_template,
            ...                         )
            >>> # Get LCZ (so we can aggregate per LCZ)
            >>> lcz_series = dataset.get_lcz()

            Now we have a dataset containing records.We now create an Analysis from it.

            >>> ana = metobs_toolkit.Analysis(orig_dataset=dataset)
            >>> ana
            Instance of Analysis at ...

            If you are interested in the diurnal cycles of your stations grouped
            per LCZ and compare them per season, we can do so by:

            >>> custom_aggdf = ana.get_aggregated_cycle_statistics(
            ...     obstype="temp",
            ...     aggregation=["lcz", "season"],
            ...     aggregation_method="mean",
            ...     horizontal_axis="hour",
            ...     stations=None,
            ...     startdt=None,
            ...     enddt=None,
            ...     plot=True,
            ...     title=None,
            ...     y_label=None,
            ...     legend=True)

            >>> custom_aggdf
                  Compact midrise ,autumn  Dense Trees (LCZ A) ,autumn  Large lowrise ,autumn  Low plants (LCZ D) ,autumn  ...  Open midrise ,autumn  Scattered Trees (LCZ B) ,autumn  Sparsely built ,autumn  Water (LCZ G) ,autumn
            hour                                                                                                           ...
            0                       17.74                     16.14                     16.95                     15.97    ...                 16.60                     15.58                          16.23                  18.22
            1                       17.23                     15.57                     16.45                     15.54    ...                 16.01                     15.11                          15.81                  17.88
            2                       16.66                     15.09                     15.90                     14.93    ...                 15.37                     14.30                          15.20                  17.57
            3                       16.42                     14.91                     15.80                     14.82    ...                 15.14                     14.16                          15.00                  17.39
            4                       16.21                     14.77                     15.64                     14.66    ...                 14.97                     13.97                          14.87                  17.21
            ...                       ...                       ...                       ...                       ...    ...                   ...                       ...                            ...                    ...
            19                      19.69                     17.29                     18.77                     17.80    ...                 19.29                     17.63                          18.02                  19.38
            20                      19.10                     16.74                     18.07                     17.06    ...                 18.44                     17.04                          17.37                  18.93
            21                      18.54                     16.46                     17.55                     16.60    ...                 17.79                     16.52                          16.97                  18.55
            22                      18.03                     16.15                     17.11                     16.27    ...                 17.24                     15.96                          16.52                  18.32
            23                      17.70                     15.89                     16.83                     15.88    ...                 16.68                     15.54                          16.16                  18.13
            <BLANKLINE>
            [24 rows x 10 columns]
        """
        # arg format checks
        startdt = self._datetime_arg_check(startdt)
        enddt = self._datetime_arg_check(enddt)

        # Check if obstype is known
        if isinstance(obstype, str):
            if obstype not in self.obstypes.keys():
                logger.error(
                    f"{obstype} is not found in the known observation types: {list(self.obstypes.keys())}"
                )
                return None
            else:
                obstype = self.obstypes[obstype]

        if _obsdf is None:
            obsdf = xs_save(self.df, obstype.name, "obstype")
            # obsdf = self.df[[obstype]]
        else:
            obsdf = xs_save(_obsdf, obstype.name, "obstype")
        assert (
            not obsdf.empty
        ), f"Error: No {obstype} observations in the analysis.df: {self.df}"

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

        # make sure all data is numeric
        values_df = values_df.astype(float)
        std_df = std_df.astype(float)

        # squize all column levels to one category for plotting
        if len(aggregation) > 1:  # more than one level for the columns
            values_df.columns = [
                " ,".join(col).strip() for col in values_df.columns.values
            ]
            std_df.columns = [" ,".join(col).strip() for col in std_df.columns.values]

        if plot:
            # generate title
            if title is None:
                startdtstr = datetime.strftime(
                    startdt, format=self.settings.app["print_fmt_datetime"]
                )
                enddtstr = datetime.strftime(
                    enddt, format=self.settings.app["print_fmt_datetime"]
                )
                title = f"{aggregation_method} - {horizontal_axis } {obstype.name} cycle for period {startdtstr} - {enddtstr} grouped by {aggregation}"

            # ylabel
            if y_label is None:
                y_label = obstype._get_plot_y_label()
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
                y_label=y_label,
                legend=legend,
                show_zero_horizontal=_show_zero_line,
            )

            if horizontal_axis == "hour":
                # extract timezone
                tzstring = str(self._get_tz())

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
        """Compute Pearson correlation coefficients.

        A method to compute the Pearson correlation between an obervation type
        and present landcover fractions in the metadf.

        The correlations are computed per group as defined by unique combinations
        of the groupby_labels.

        A dictionary is returnd where each key represents a unique combination of
        the groupby_labels. The value is a dictionary with the following keys
        and values:

        * cor matrix: the Pearson correlation matrix
        * significance matrix: the significance (p-)values of the correlations.
        * combined matrix: A human-readable combination of the correlations and their p values. Indicate by \*, \*\* or \*\*\* representing p-values < 0.05, 0.01 and 0.001 respectively.

        This dictionary is also stored as a lc_cor_dict attribute.

        Parameters
        ----------
        obstype : str, or list optional
            The observation type(s) to compute the correlations on. The default is ['temp'].
        groupby_labels : list, optional
            List of variables to form one group, resulting in one correlation.
            These variables should be either a categorical observation type, a categorical column in the metadf or
            a time aggregation. All possible time aggreagetions are: ['minute',
            'hour', 'month', 'year', 'day_of_year',
            'week_of_year', 'season']. The default is ['hour'].

        Returns
        -------
        cor_dict : dict
            A nested dictionary with unique combinations of groupby values.

        See Also
        ---------
        Dataset.get_landcover: Get landcover fractions at your stations.
        Analysis.plot_correlation_heatmap: Make a heatmap of a correlation matrix
        Analysis.plot_correlation_variation: Scatter plot of correlations in time.

        Examples
        --------

        An Analysis is always created from a Dataset, so we start by creating
        a Dataset and importing data into it.

        >>> import metobs_toolkit
        >>>
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> dataset.import_data_from_file(
        ...                         input_data_file=metobs_toolkit.demo_datafile,
        ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                         template_file=metobs_toolkit.demo_template,
        ...                         )
        >>> # Get landcover fraction (to calculate correlations for)
        >>> landcover_df = dataset.get_landcover(buffers=[250, 500], aggregate=True)

        Now we have a dataset containing records.We now create an Analysis from it.

        >>> ana = metobs_toolkit.Analysis(orig_dataset=dataset)
        >>> ana
        Instance of Analysis at ...

        As a demo, we calculate landcover correlations with temperature, and we
        do this for each hour (thus aggregated to the hour --> diurnal characteristics).

        >>> cordict = ana.get_lc_correlation_matrices(obstype=["temp"], groupby_labels=["hour"])

        Inspect the cordict to find the correlation matrices (and the corresponding
        significance matrices) for each group (the hours in this example).

        >>> # The tempeature correlation matrix at 4 UTC
        >>> cordict[4]['cor matrix']
                         temp  water_250m  water_500m  pervious_250m  pervious_500m  impervious_250m  impervious_500m
        temp             1.00        0.12        0.13          -0.33          -0.32             0.25             0.23
        water_250m       0.12        1.00        0.97          -0.32          -0.18            -0.33            -0.35
        water_500m       0.13        0.97        1.00          -0.28          -0.20            -0.34            -0.35
        pervious_250m   -0.33       -0.32       -0.28           1.00           0.95            -0.79            -0.76
        pervious_500m   -0.32       -0.18       -0.20           0.95           1.00            -0.83            -0.85
        impervious_250m  0.25       -0.33       -0.34          -0.79          -0.83             1.00             0.98
        impervious_500m  0.23       -0.35       -0.35          -0.76          -0.85             0.98             1.00

        There is also a matrix combining the correlations and the significance,
        using a star-presentation.

        >>> cordict[4]['combined matrix']
                              temp water_250m water_500m pervious_250m pervious_500m impervious_250m impervious_500m
        temp              1.00 ***   0.12 ***   0.13 ***     -0.33 ***     -0.32 ***        0.25 ***        0.23 ***
        water_250m        0.12 ***   1.00 ***   0.97 ***     -0.32 ***     -0.18 ***       -0.33 ***       -0.35 ***
        water_500m        0.13 ***   0.97 ***   1.00 ***     -0.28 ***     -0.20 ***       -0.34 ***       -0.35 ***
        pervious_250m    -0.33 ***  -0.32 ***  -0.28 ***      1.00 ***      0.95 ***       -0.79 ***       -0.76 ***
        pervious_500m    -0.32 ***  -0.18 ***  -0.20 ***      0.95 ***      1.00 ***       -0.83 ***       -0.85 ***
        impervious_250m   0.25 ***  -0.33 ***  -0.34 ***     -0.79 ***     -0.83 ***        1.00 ***        0.98 ***
        impervious_500m   0.23 ***  -0.35 ***  -0.35 ***     -0.76 ***     -0.85 ***        0.98 ***        1.00 ***
        """
        if not isinstance(obstype, list):
            obstype = [obstype]

        # get data
        df = self.get_analysis_records()
        df = df[obstype].reset_index()
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
            p_stars = pval.map(
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
        """Make a heatmap plot of a correlation matrix.

        To specify which correlation matrix to plot, specify the group value
        using the groupby_value argument.

        All possible groupby_values are the keys of the lc_cor_dict attribute.

        Parameters
        ----------
        groupby_value : str, num, None, optional
            A groupby value to indicate which correlation matrix to visualize.
            If None is given, the first groupby value that is present is
            chosen. The default is None.
        title : str, optional
            Title of the figure. If None, a default title is constructed. The
            default is None.

        Returns
        -------
        None.

        See Also
        ---------
        Dataset.get_landcover: Get landcover fractions at your stations.
        Analysis.get_lc_correlation_matrices: Make a heatmap of a correlation matrix
        Analysis.plot_correlation_variation: Scatter plot of correlations in time.

        Note
        ------
        To list all possible groupby_values, one can use
        ` print(Analysis_instance.lc_cor_dict.keys())`

        Examples
        --------

        .. plot::
            :context: close-figs

            An Analysis is always created from a Dataset, so we start by creating
            a Dataset and importing data into it.

            >>> import metobs_toolkit
            >>>
            >>> #Create your Dataset
            >>> dataset = metobs_toolkit.Dataset() #empty Dataset
            >>> dataset.import_data_from_file(
            ...                         input_data_file=metobs_toolkit.demo_datafile,
            ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
            ...                         template_file=metobs_toolkit.demo_template,
            ...                         )
            >>> # Get landcover fraction (to calculate correlations for)
            >>> landcover_df = dataset.get_landcover(buffers=[250, 500], aggregate=True)

            Now we have a dataset containing records.We now create an Analysis from it.

            >>> ana = metobs_toolkit.Analysis(orig_dataset=dataset)
            >>> ana
            Instance of Analysis at ...

            As a demo, we calculate landcover correlations with temperature, and we
            do this for each hour (thus aggregated to the hour --> diurnal characteristics).

            >>> cordict = ana.get_lc_correlation_matrices(obstype=["temp"], groupby_labels=["hour"])

            Now we can make a heatmap plot of the correlations at an hour of choice.

            >>> ana.plot_correlation_heatmap(groupby_value=5) #cor plot at 5 UTC


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

        Make a scatter plot of the correlations to visualize differences between
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

        See Also
        ---------
        Dataset.get_landcover: Get landcover fractions at your stations.
        Analysis.get_lc_correlation_matrices: Make a heatmap of a correlation matrix
        Analysis.plot_correlation_heatmap: Make a heatmap of a correlation matrix

        Note
        ------
        If to many possible group values exist, one can use the
        get_full_dataframe(), filter the dataframe and set_data() method to
        reduce the group values.

        Examples
        --------

        .. plot::
            :context: close-figs

            An Analysis is always created from a Dataset, so we start by creating
            a Dataset and importing data into it.

            >>> import metobs_toolkit
            >>>
            >>> #Create your Dataset
            >>> dataset = metobs_toolkit.Dataset() #empty Dataset
            >>> dataset.import_data_from_file(
            ...                         input_data_file=metobs_toolkit.demo_datafile,
            ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
            ...                         template_file=metobs_toolkit.demo_template,
            ...                         )
            >>> # Get landcover fraction (to calculate correlations for)
            >>> landcover_df = dataset.get_landcover(buffers=[250, 500], aggregate=True)

            Now we have a dataset containing records.We now create an Analysis from it.

            >>> ana = metobs_toolkit.Analysis(orig_dataset=dataset)
            >>> ana
            Instance of Analysis at ...

            As a demo, we calculate landcover correlations with temperature, and we
            do this for each hour (thus aggregated to the hour --> diurnal caracteristics).

            >>> cordict = ana.get_lc_correlation_matrices(obstype=["temp"], groupby_labels=["hour"])

            Now we can make an avolution plot of the landcover correlations.

            >>> ana.plot_correlation_variation()
            <Axes: title={'center': "Correlation scatter for group: ['hour']"}, xlabel="Groups of ['hour']", ylabel='Pearson correlation'>

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
    """Construct time-derivated columns if required.

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
        An obtained dataframe where a label for the seasons has been added.
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


# =============================================================================
# Docstring test
# =============================================================================
if __name__ == "__main__":
    from metobs_toolkit.doctest_fmt import setup_and_run_doctest

    setup_and_run_doctest()
