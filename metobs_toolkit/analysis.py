#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module contains the Analysis class and all its methods.

A Analysis holds a set of 'good' observations and the methods will analyse them.
"""
from datetime import datetime
import pandas as pd
import numpy as np
import copy
from scipy.stats import pearsonr

from metobs_toolkit.plotting_functions import (diurnal_plot,
                                               heatmap_plot,
                                               correlation_scatter)

from metobs_toolkit.df_helpers import (init_multiindexdf,
                                        datetime_subsetting,
                                        subset_stations,
                                        fmt_datetime_argument)

class Analysis():
    """ The Analysis class contains methods for analysing diurnal cycles and landcover effects"""
    def __init__(self, obsdf, metadf, settings, data_template):
        self.df = obsdf
        self.metadf = metadf
        self.settings = settings
        self.data_template = data_template

        # analysis objects
        self.lc_cor_dict = {}
        self._lc_cor_obstype = None
        self._lc_groupby_labels = None

    def __str__(self):
        if self.df.empty:
            return f"Empty Analysis instance."
        add_info = ''
        n_stations = self.df.index.get_level_values('name').unique().shape[0]
        n_obs_tot = self.df.shape[0]

        startdt = self.df.index.get_level_values('datetime').min()
        enddt = self.df.index.get_level_values('datetime').max()

        if ((not self.metadf['lat'].isnull().all()) &
            (not self.metadf['lon'].isnull().all())):
            add_info += '     *Coordinates are available for all stations. \n'

        if (not self.metadf['lcz'].isnull().all()):
            add_info += "     *LCZ's are available for all stations. \n"

        if bool(self.lc_cor_dict):
            add_info += f"     *landcover correlations are computed on group: {self._lc_groupby_labels}  \n"




        return (f"Analysis instance containing: \n \
    *{n_stations} stations \n \
    *{self.df.columns.to_list()} observation types \n \
    *{n_obs_tot} observation records \n{add_info} \n \
    *records range: {startdt} --> {enddt} (total duration:  {enddt - startdt})" + add_info)

    def __repr__(self):
        return self.__str__()

    # =============================================================================
    #     Setters
    # =============================================================================


    def subset_period(self, startdt, enddt):
        """
        Subset the observations of the Analysis to a specific period. The same
        timezone is assumed as the data.

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
        if not isinstance(startdt, type(datetime(2020,1,1))):
            print(f' {startdt} not a datetime type. Ignore subsetting!')
            return
        if not isinstance(enddt, type(datetime(2020,1,1))):
            print(f' {enddt} not a datetime type. Ignore subsetting!')
            return

        startdt = fmt_datetime_argument(startdt, self.settings.time_settings['timezone'])
        enddt = fmt_datetime_argument(enddt, self.settings.time_settings['timezone'])

        self.df = datetime_subsetting(self.df, startdt, enddt)

    # =============================================================================
    #   Helpers
    # =============================================================================


    def apply_filter(self, expression):
        """
        Method to filter an Analysis by a user definde string expression. This
        can be used to filter the observation to specific meteorological conditions
        (i.e. low windspeeds, high humidity, cold temperatures, ...)

        The filter expression contains only columns present in the Analysis.df
        and/or the Analysis.metadf.

        A New Analysis object is returned.

        Parameters
        ----------

        expression : str
            A filter expression using columnnames present in either df or metadf.
            The following timestamp derivatives can be used as well: [minute, hour,
            month, year, day_of_year, week_of_year, season]. The quarry_str may
            contain number and expressions like <, >, ==, >=, *, +, .... Multiple filters
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
        Make shure to use " of ' to indicate string values in the expression if
        needed.


        """

        child_df, child_metadf = filter_data(df = self.df,
                                             metadf = self.metadf,
                                             quarry_str = expression)

        return Analysis(obsdf=child_df,
                        metadf = child_metadf,
                        settings = self.settings,
                        data_template = self.data_template)




    def aggregate_df(self, df=None, agg=['lcz', 'hour'], method='mean'):
        """
        Aggregate observations to a (list of) categories.

        The output will be a dataframe that is aggregated to one, or more categories.
        A commen example is aggregating to LCZ's.


        Parameters
        ----------
        df : pandas.DataFrame or None
            The observations to aggregate. If None, the df attribute of the
            Analysis instance is used. The default is None.
        agg : list, optional
            The list of columnnames to aggregate to. If 'lcz' is included, the
            lcz information is extracted from the Analysis.metadf. The default is ['lcz', 'datetime'].
        method : str, optional
            list of functions and/or function names, e.g. [np.sum, 'mean']. The default is 'mean'.

        Returns
        -------
        pandas.DataFrame
            A dataframe with the agg columns as an index. The values are the aggregated values.

        Note
        -------
        Present columns that ar non-numeric and are not in the agg list, are
        not present in the return, since these values cannot be aggregated.

        """

        if df is None:
            df = copy.deepcopy(self.df)
        df = df.reset_index()

        time_agg_keys = ['minute', 'hour', 'month', 'year', 'day_of_year',
                         'week_of_year', 'season']

        # scan trough the metadf for aggregation keys
        for agg_key in agg:
            if agg_key not in df.columns:
                # look in metadf
                if agg_key in self.metadf.columns:
                    df = pd.merge(df, self.metadf[[agg_key]],
                                      how='left', left_on='name',
                                      right_index=True)




        # Check if all agg keys are present or defined:
        possible_agg_keys = time_agg_keys
        possible_agg_keys.extend(list(df.columns))
        unmapped = [agg_key for agg_key in agg if agg_key not in possible_agg_keys]
        assert len(unmapped) == 0, f'cannot aggregate to unknown labels: {unmapped}.'


        # make time-derivate columns if required
        df = _make_time_derivatives(df, agg)


        # check if not all values are Nan
        for agg_name in agg:
            assert df[agg_name].isnull().all() == False, f'Aggregation to {agg_name} not possible because no valid values found for {agg_name}.'



        # Aggregate the df
        agg_df = df.groupby(agg).agg(method, numeric_only=True)
        # sort index
        agg_df = agg_df.reset_index()
        agg_df = agg_df.set_index(agg)
        return agg_df


    # =============================================================================
    #   Analyse method
    # =============================================================================
    def get_anual_statistics(self, groupby=['name'], obstype='temp',
                             agg_method='mean', stations=None,
                             startdt=None, enddt=None, plot=True,
                             errorbands=False, title=None, y_label=None,
                             legend=True):

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


        obsdf = self.df[[obstype]]

        # Filter stations
        if not stations is None:
            if isinstance(stations, str):
                stations = [stations]

            obsdf = subset_stations(obsdf, stations)

        # Filter datetimes
        obsdf = datetime_subsetting(df=obsdf,
                                    starttime=startdt,
                                    endtime=enddt)

        # Define aggregation groups
        agg = groupby.copy()
        agg.insert(0, 'month')


        # get agg statistics
        df = self.aggregate_df(df=obsdf, agg=agg, method=[agg_method, 'std'])
        # multiindex to string format
        df = df.reset_index()
        for idx_col in agg:
            df[idx_col] = df[idx_col].astype(str)
        df = df.set_index(agg)

        # unstack multiindex and unstack column levels if multiple groupby labels
        aggdf = df.unstack(groupby).droplevel(0, axis='columns')
        values_df = aggdf[[agg_method]].droplevel(0, axis='columns') #df with values
        std_df = aggdf[['std']].droplevel(0, axis='columns') #df with errors

        def _format_for_plotting(df):

            if len(groupby) > 1:
                df.columns = [' ,'.join(col).strip() for col in df.columns.values]

            # sort months order
            months = ['January', 'February', 'March', 'April', 'May', 'June',
                      'July', 'August', 'September', 'October', 'November',
                      'December']
            df = df.reindex(months)
            return df

        plotdf = _format_for_plotting(values_df)
        if errorbands:
            errordf= _format_for_plotting(std_df)
        else:
            errordf=None

        # title
        desc_dict = self.data_template[obstype].to_dict()

        if title is None:
            title=f'Anual {desc_dict["description"]} cycle plot per {groupby}.'
        else:
            title = str(title)

        # ylabel
        if y_label is None:

            y_label = f'{desc_dict["description"]} ({desc_dict["units"]})'
        else:
            y_label = str(y_label)



        if plot:
            ax = diurnal_plot(diurnaldf = plotdf,
                              errorbandsdf = errordf,
                              title = title,
                              plot_settings = self.settings.app['plot_settings']['anual'],
                              colorby = 'name',
                              lcz_dict = None,
                              data_template=self.data_template,
                              obstype=obstype,
                              y_label = y_label,
                              legend=legend)
        return plotdf


    def get_diurnal_statistics(self, obstype='temp', stations=None,
                               startdt=None, enddt=None, plot=True,
                               title=None, y_label=None, legend=True,
                               colorby='name', errorbands=False,
                               verbose=False):
        """
        Create an average diurnal cycle for the observations.

        (In the plot, each station is represed by a line.)


        Parameters
        ----------
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
        colorby : 'name' or 'lcz', optional
            If 'name' the plotted lines will be colored per station, if 'lcz' the colors represent the stations lcz. The default is 'name'.
        errorbands : bool, optional
            If True, the std is representd in the plot by colored bands. The default is False.
        verbose : True, optional
            If True, the dataframse with aggregation information are returned . The default is False.

        Returns
        -------
        tuple (if verbose)
            A tuple of dataframes is returned when verbose is True.

        Note
        --------
        If a timezone unaware datetime is given as an argument, it is interpreted
        as if it has the same timezone as the observations.

        """

        obsdf = self.df

        # Filter stations
        if not stations is None:
            if isinstance(stations, str):
                stations = [stations]

            obsdf = subset_stations(obsdf, stations)

        # Filter datetimes
        obsdf = datetime_subsetting(df=obsdf,
                                    starttime=startdt,
                                    endtime=enddt)

        # check if lcz is available if required
        if colorby == 'lcz':
            if self.metadf['lcz'].isnull().any():
                print("ERROR: Not all stations have a LCZ. Update the LCZ's first or use colorby='name'. ")
                return None



        # Get hours for all records
        obsdf = obsdf.reset_index()
        obsdf['hour'] = obsdf['datetime'].dt.hour


        agg_column_name = obstype #aggregate the measured obstypes
        startdt = obsdf['datetime'].min()
        enddt = obsdf['datetime'].max()


        # groupby and take the mean per station per hour.
        stats = obsdf.groupby(['name', 'hour'])[agg_column_name].agg(['mean', 'std', 'median'])

        hourly_avg = stats['mean'].unstack().transpose()


        if plot:
            # get lcz groups if needed
            if ((colorby == 'lcz') & (not 'lcz' in self.metadf.columns)):
                print('Warning: No LCZ information, thus colorby will be set to name.')
                colorby = 'name'
            if colorby =='lcz':
                lcz_dict = self.metadf['lcz'][hourly_avg.columns.to_list()].to_dict()
            else:
                lcz_dict = None

            # generate title
            if title is None:
                startdtstr = datetime.strftime(startdt, format=self.settings.app["print_fmt_datetime"])
                enddtstr = datetime.strftime(enddt, format=self.settings.app["print_fmt_datetime"])

                title=f'Hourly average {obstype} diurnal cycle for period {startdtstr} - {enddtstr}'


            # generate errorbands df
            if errorbands:
                stddf = stats['std'].unstack().transpose()
            else:
                stddf = None

            # extract timezone
            tzstring = self.df.index.get_level_values('datetime').tz.zone


            # Make plot
            ax = diurnal_plot(diurnaldf = hourly_avg,
                         errorbandsdf = stddf,
                         title = title,
                         plot_settings = self.settings.app['plot_settings']['diurnal'],
                         colorby = colorby,
                         lcz_dict = lcz_dict,
                         data_template=self.data_template,
                         obstype=obstype,
                         y_label = y_label,
                         legend=legend)

            ax.xaxis.set_major_formatter('{x:.0f} h')
            ax.set_xlabel(f'Hours (timezone: {tzstring})')


        if verbose:
            return hourly_avg, stats

        return hourly_avg

    def get_diurnal_statistics_with_reference(self, refstation, obstype='temp',
                                              tollerance='30T', stations=None,
                                              startdt=None, enddt=None,
                                              plot=True, title=None,
                                              y_label=None, legend=True,
                                              colorby='name', errorbands=False,
                                              verbose=False):
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
        obstype : str, optional
            Element of the metobs_toolkit.observation_types The default is 'temp'.
        tollerance : Timedelta or str, optional
            The tollerance string or object representing the maximum translation in time to find a reference
            observation for each observation. Ex: '5T' is 5 minuts, '1H', is one hour. The default is '30T'.
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
        colorby : 'name' or 'lcz', optional
            If 'name' the plotted lines will be colored per station, if 'lcz' the colors represent the stations lcz. The default is 'name'.
        errorbands : bool, optional
            If True, the std is representd in the plot by colored bands. The upper bound represents +1 x std, the lower bound -1 x std. The default is False.
        verbose : True, optional
            If True, the dataframse with aggregation information are returned . The default is False.

        Returns
        -------
        tuple (if verbose)
            A tuple of dataframes is returned when verbose is True.

        Note
        --------
        If a timezone unaware datetime is given as an argument, it is interpreted
        as if it has the same timezone as the observations.

        """




        obsdf = self.df

        # Filter stations
        if not stations is None:
            if isinstance(stations, str):
                stations = [stations]
                stations.append(refstation)
            obsdf = subset_stations(obsdf, stations)


        # Filter datetimes
        obsdf = datetime_subsetting(df=obsdf,
                                    starttime=startdt,
                                    endtime=enddt)

        # check if lcz is available if required
        if colorby == 'lcz':
            if self.metadf['lcz'].isnull().any():
                print("ERROR: Not all stations have a LCZ. Update the LCZ's first or use colorby='name'. ")
                return None

        obsdf = obsdf[obstype].reset_index()



        # extract refernce from observations
        refdf = obsdf[obsdf['name'] == refstation]
        obsdf = obsdf[obsdf['name']!= refstation]

        # Syncronize observations with the reference observations
        refdf = refdf.rename(columns={obstype: 'ref_'+obstype, 'datetime': 'ref_datetime'})
        mergedf = pd.merge_asof(left=obsdf.sort_values('datetime'),
                                right=refdf[['ref_datetime', 'ref_'+obstype]].sort_values('ref_datetime'),
                                right_on="ref_datetime",
                                left_on="datetime",
                                direction="nearest",
                                tolerance=pd.Timedelta(tollerance),
                                )

        startdt = refdf['ref_datetime'].min()
        enddt = refdf['ref_datetime'].max()

        # Compute difference
        agg_column_name = 'difference'
        mergedf[agg_column_name] = mergedf[obstype] - mergedf['ref_'+obstype]

        # Get hour column
        mergedf['hour'] = mergedf['datetime'].dt.hour

        # overwrite the obsdf
        obsdf = mergedf
        # groupby and take the mean per station per hour.
        stats = obsdf.groupby(['name', 'hour'])[agg_column_name].agg(['mean', 'std', 'median'])

        hourly_avg = stats['mean'].unstack().transpose()


        if plot:
            # get lcz groups if needed
            if ((colorby == 'lcz') & (not 'lcz' in self.metadf.columns)):
                print('Warning: No LCZ information, thus colorby will be set to name.')
                colorby = 'name'

            if colorby =='lcz':
                lcz_dict = self.metadf['lcz'][hourly_avg.columns.to_list()].to_dict()
            else:
                lcz_dict = None

            # generate title
            if title is None:
                startdtstr = datetime.strftime(startdt, format=self.settings.app["print_fmt_datetime"])
                enddtstr = datetime.strftime(enddt, format=self.settings.app["print_fmt_datetime"])
                title=f'Hourly average {obstype} diurnal cycle, with {refstation} as reference, for period {startdtstr} - {enddtstr}'


            # generate errorbands df
            if errorbands:
                stddf = stats['std'].unstack().transpose()
            else:
                stddf = None

            # extract timezone
            tzstring = self.df.index.get_level_values('datetime').tz.zone


            # Make plot
            ax = diurnal_plot(diurnaldf = hourly_avg,
                             errorbandsdf = stddf,
                             title = title,
                             plot_settings = self.settings.app['plot_settings']['diurnal'],
                             colorby = colorby,
                             lcz_dict = lcz_dict,
                             data_template=self.data_template,
                             obstype = obstype,
                             y_label = y_label,
                             legend=legend,
                             show_zero_horizontal = True)
            ax.xaxis.set_major_formatter('{x:.0f} h')
            ax.set_xlabel(f'Hours (timezone: {tzstring})')


        if verbose:
            return hourly_avg, stats, obsdf

        return hourly_avg


    def get_aggregated_diurnal_statistics(self, obstype='temp', stations=None,
                                          aggregation=['lcz', 'datetime'],
                                          aggregation_method='mean',
                                          startdt=None, enddt=None, plot=True,
                                          title=None, y_label=None, legend=True,
                                          errorbands=False, verbose=False):

        """
        Create an average diurnal cycle for an aggregated categorie. A commen
        example is to aggregate to the LCZ's, so to get the diurnal cycle per LCZ
        rather than per station.

        (In the plot, each aggregated category different from datetime, is represed by a line.)

        Parameters
        ----------
        obstype : str, optional
            Element of the metobs_toolkit.observation_types The default is 'temp'.
        stations : list, optional
            List of station names to use. If None, all present stations will be used. The default is None.
        aggregation : list, optional
            List of variables to aggregate to. These variables should either a
            categorical observation type, a categorical column in the metadf or
            a time aggregation. All possible time aggreagetions are: ['minute',
            'hour', 'month', 'year', 'day_of_year',
            'week_of_year', 'season']. The default is ['lcz', 'datetime'].
        aggregation_method : str, optional
            Which (numpy) function is used to aggregate the observations. The default is 'mean'.
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
            If True, the dataframse with aggregation information are returned . The default is False.

        Returns
        -------
        tuple (if verbose)
            A tuple of dataframes is returned when verbose is True.

        Note
        -------
        If a timezone unaware datetime is given as an argument, it is interpreted
        as if it has the same timezone as the observations.

        """

        obsdf = self.df

        # Filter stations
        if not stations is None:
            if isinstance(stations, str):
                stations = [stations]

            obsdf = subset_stations(obsdf, stations)


        # Filter datetimes
        obsdf = datetime_subsetting(df=obsdf,
                                    starttime=startdt,
                                    endtime=enddt)


        if bool(aggregation):
            # check if datetime is in the aggreagation, otherwise no time component is left
            if not 'datetime' in aggregation:
                print(f'WARNING: To make a diurnal cycle with aggregation, the "datetime" must be in the aggregation list. "datetime" is added to it.')
                aggregation.append('datetime')
            obsdf = self.aggregate_df(df = obsdf, agg=aggregation,
                                      method=aggregation_method)

        obsdf = obsdf.reset_index()
        # Create identifiers to form unique hours
        obsdf['hour'] = obsdf['datetime'].dt.hour


        # aggregation scheme setup
        agg_column_name = obstype #aggregate the measured obstypes
        groupby_list = aggregation
        groupby_list.append('hour')
        groupby_list.remove('datetime')

        # for plot titles
        startdt = obsdf['datetime'].dropna().min()
        enddt = obsdf['datetime'].dropna().max()


        # groupby and take the mean per station per hour.
        stats = obsdf.groupby(groupby_list)[agg_column_name].agg(['mean', 'std', 'median'])

        hourly_avg = stats['mean'].unstack().transpose()

        if plot:

            # generate title
            if title is None:
                startdtstr = datetime.strftime(startdt, format=self.settings.app["print_fmt_datetime"])
                enddtstr = datetime.strftime(enddt, format=self.settings.app["print_fmt_datetime"])
                title=f'Hourly average {obstype} diurnal cycle for period {startdtstr} - {enddtstr} grouped by {groupby_list}'



            # generate errorbands df
            if errorbands:
                stddf = stats['std'].unstack().transpose()
            else:
                stddf = None

            # extract timezone
            tzstring = self.df.index.get_level_values('datetime').tz.zone


            # Make plot
            ax = diurnal_plot(diurnaldf = hourly_avg,
                             errorbandsdf = stddf,
                             title = title,
                             plot_settings = self.settings.app['plot_settings']['diurnal'],
                             colorby = 'name',
                             lcz_dict = None,
                             data_template=self.data_template,
                             obstype=obstype,
                             y_label = y_label,
                             legend=legend)

            ax.xaxis.set_major_formatter('{x:.0f} h')
            ax.set_xlabel(f'Hours (timezone: {tzstring})')

        if verbose:
            return hourly_avg, stats

        return hourly_avg

    # =============================================================================
    # Correlations analysis
    # =============================================================================

    def get_lc_correlation_matrices(self, obstype=['temp'], groupby_labels=['hour']):
        """
        A method to compute the Pearson correlation between an obervation type
        and present landcover fractions in the metadf.

        The correlations are computed per group as defined by unique combinations
        of the groupby_labels.

        A dictionary is returnd where each key represents a unique combination of
        the groupby_labels. The value is a dictionary with the following keys
        and values:
        * cor matrix: the Pearson correlation matrix
        * significance matrix: the significance (p-)values of the correlations.
        * combined matrix: A human readable combination of the correlations and their p values. Indicate by *, ** or *** representing p-values < 0.05, 0.01 and 0.001 respectively.

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
            'week_of_year', 'season']. The default is ['lcz', 'datetime'].. The default is ['hour'].

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
                df = df.merge(self.metadf[[group_lab]],
                              how='left',
                              left_on='name',
                              right_index=True)

        for group_lab in groupby_labels:
            assert group_lab in df.columns, f'"{group_lab}" is found in the observations of possible groupby_labels.'


        # subset columns
        relev_columns = [label for label in groupby_labels] #to avoid deep copy import
        relev_columns.append('name')
        relev_columns.extend(obstype)
        df = df[relev_columns]

        # find landcover columnnames in the metadf
        lc_columns = [col for col in self.metadf.columns if (('_' in col ) & (col.endswith('m')))]

        # get landcover data
        lc_df = self.metadf[lc_columns]

        if lc_df.empty:
            print('WARNING: No landcover columns found in the metadf. Landcover correlations cannot be computed.')
            return None


        # merge together
        df = df.merge(lc_df, how='left', left_on='name', right_index = True)

        # remove name column if it is not explicit in the groupby labels
        if 'name' not in groupby_labels:
            df = df.drop(columns=['name'])



        # create return
        cor_dict = {}

        #Iterate over all groups

        # avoid futur pandas warning for groupby labels of len==1
        if len(groupby_labels) == 1:
            groups = df.groupby(groupby_labels[0])
        else:
            groups = df.groupby(groupby_labels)


        for group_lab, groupdf in groups:
            # No correlations can be computed when no variance is found
            if groupdf.shape[0] <=1:
                print(f'No variance found in correlationd group {group_lab}. Correlation thus not be computed for this group: {groupdf}.')
                continue
            # drop groupby labels
            groupdf = groupdf.drop(columns=groupby_labels, errors='ignore')

            rho = groupdf.corr(method='pearson')
            pval = groupdf.corr(method=lambda x, y: pearsonr(x, y)[1]) - np.eye(*rho.shape)
            # represent p values by stars
            p_stars = pval.applymap(lambda x: ''.join(['*' for t in [.05, .01, .001] if x<=t]))

            # combined human readable df
            comb_df = pd.DataFrame(index=rho.index)
            for col in rho.columns:
                comb_df[col] = rho[col].apply(lambda x: f"{x:.02f}") + ' ' + p_stars[col]

            cor_dict[group_lab] = {'cor matrix': rho,
                                   'significance matrix': pval,
                                   'combined matrix': comb_df}


        # Update attribute
        self.lc_cor_dict = cor_dict
        self._lc_cor_obstype = obstype
        self._lc_groupby_labels = groupby_labels

        return cor_dict


    def plot_correlation_heatmap(self, groupby_value=None, title=None):
        """
        Make a heatmap plot af a correaltion matrix. To specify which correlation
        matrix to plot, specify the group value using the groupby_value argument.

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
        assert bool(self.lc_cor_dict), 'No correlation matrices found, use the metod get_lc_correlation_matrices first.'

        if groupby_value is None:
            groupby_value = list(self.lc_cor_dict.keys())[0]
            print('WARNING: No groupby_value is given, so the first groupby value (={groupby_value}) will be used!')
            print(f'INFO: The correlations are computed over {self._lc_groupby_labels} with the following unique values: {list(self.lc_cor_dict.keys())}')

        # check if groupby value exists
        assert groupby_value in self.lc_cor_dict.keys(), f'{groupby_value} not found as a groupby value. These are all the possible values: {self.lc_cor_dict.keys()}'


        if title is None:
            title = f'Correlation heatmap for group: {self._lc_groupby_labels} = {groupby_value}'

        heatmap_plot(cor_dict = self.lc_cor_dict[groupby_value],
                     title=title,
                     heatmap_settings = self.settings.app['plot_settings']['correlation_heatmap'])


    def plot_correlation_variation(self, title=None):
        """
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
        assert bool(self.lc_cor_dict), 'No correlation matrices found, use the metod get_lc_correlation_matrices first.'

        # check if correlation evolution information is available
        if len(self.lc_cor_dict.keys()) <= 1:
            print(f'WARNING: Only one correlation group is found: {self.lc_cor_dict.keys()}')
            print('The variance plot can not be made.')
            return



        if title is None:
            title = f'Correlation scatter for group: {self._lc_groupby_labels}'


        correlation_scatter(full_cor_dict = self.lc_cor_dict,
                            groupby_labels = self._lc_groupby_labels,
                            obstypes =self._lc_cor_obstype,
                            title=title,
                            cor_scatter_settings = self.settings.app['plot_settings']['correlation_scatter'])



def _make_time_derivatives(df, required, get_all=False):
    """ construct time derivated columns if required.
        datetime must be a column."""

    if ('minute' in required) | (get_all):
        df['minute'] = df['datetime'].dt.minute
    if ('hour' in required) | (get_all):
        df['hour'] = df['datetime'].dt.hour
    if ('month' in required) | (get_all):
        df['month'] = df['datetime'].dt.month_name()
    if ('year' in required) | (get_all):
        df['year'] = df['datetime'].dt.year
    if ('day_of_year' in required) | (get_all):
        df['day_of_year'] = df['datetime'].dt.day_of_year
    if ('week_of_year' in required) | (get_all):
        df['week_of_year'] = df['datetime'].dt.isocalendar()['week']
    if ('season' in required) | (get_all):
        df['season'] = get_seasons(df['datetime'])

    return df



def get_seasons(datetimeseries,
                start_day_spring = '01/03' ,
                start_day_summer = '01/06',
                start_day_autumn = '01/09',
                start_day_winter = '01/12'):

    """ Convert a datetimeseries to a season label (i.g. categorical). """


    spring_startday = datetime.strptime(start_day_spring, '%d/%m')
    summer_startday = datetime.strptime(start_day_summer, '%d/%m')
    autumn_startday = datetime.strptime(start_day_autumn, '%d/%m')
    winter_startday = datetime.strptime(start_day_winter, '%d/%m')


    seasons = pd.Series(index=['spring', 'summer', 'autumn', 'winter'],
                        data=[spring_startday, summer_startday, autumn_startday, winter_startday],
                        name='startdt').to_frame()
    seasons['day_of_year'] = seasons['startdt'].dt.day_of_year - 1

    bins = [0]
    bins.extend(seasons['day_of_year'].to_list())
    bins.append(366)

    labels = ['winter', 'spring', 'summer', 'autumn', 'winter']



    return pd.cut(x = datetimeseries.dt.day_of_year,
                  bins = bins,
                  labels=labels,
                  ordered=False,
                  )







def filter_data(df, metadf, quarry_str):
    """
    Function to filter a dataframe by a user definde string expression. This
    can be used to filter the observation to specific meteorological conditions
    (i.e. low windspeeds, high humidity, cold temperatures, ...)

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
        contain number and expressions like <, >, ==, >=, *, +, .... Multiple filters
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
    df = _make_time_derivatives(df, required=' ', get_all=True)

    # merge together on name
    mergedf = df.merge(metadf, how='left', on='name')

    #apply filter
    filtered = mergedf.query(expr=quarry_str)

    # split to df and metadf
    filter_df = filtered[df_init_cols]
    filter_metadf = filtered[metadf_init_cols]

    # set indexes
    filter_df = filter_df.set_index(df_init_idx)
    filter_metadf = filter_metadf.set_index(metadf_init_idx)

    return filter_df, filter_metadf