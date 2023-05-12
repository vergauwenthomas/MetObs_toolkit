#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 15:35:07 2023

@author: thoverga
"""
from datetime import datetime
import pandas as pd

from metobs_toolkit.plotting_functions import diurnal_plot

from metobs_toolkit.df_helpers import (init_multiindexdf,
                                        datetime_subsetting,
                                        subset_stations)

class Analysis():
    """ The Analysis class contains methods for analysing diurnal cycles and landcover effects"""
    def __init__(self, obsdf, metadf, settings, data_template):
        self.df = obsdf
        self.metadf = metadf
        self.settings = settings
        self.data_template = data_template



    # =============================================================================
    #     Setters
    # =============================================================================


    def subset_period(self, startdt, enddt):
        """
        Subset the observations of the Analysis to a specific period.

        Parameters
        ----------
        startdt : datetime.datetime
            The start datetime to filter the observations to.
        enddt : datetime.datetime
            The end datetime to filter the observations to.

        Returns
        -------
        None.

        """
        if not isinstance(startdt, type(datetime(2020,1,1))):
            print(f' {startdt} not a datetime type. Ignore subsetting!')
            return
        if not isinstance(enddt, type(datetime(2020,1,1))):
            print(f' {enddt} not a datetime type. Ignore subsetting!')
            return

        self.df = datetime_subsetting(self.df, startdt, enddt)

    # =============================================================================
    #   Helpers
    # =============================================================================


    def aggregate_df(self, df, agg=['lcz', 'datetime'], method='mean'):
        """
        Aggregate observations to a (list of) categories.

        The output will be a dataframe that is aggregated to one, or more categories.
        A commen example is aggregating to LCZ's.


        Parameters
        ----------
        df : pandas.DataFrame
            The observations to aggregate.
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
        Present columns that ar non-numeric and are not in the agg list are not present in the return,
        since these values cannot be aggregated.

        """
        df = df.reset_index()

        # merge relevant info to the df for aggregation

        if 'lcz' in agg:
            if not 'lcz' in self.metadf:
                print('Warning: Aggregation to LCZ not possible because no LCZ information found.')
                return df
            else:
                df = pd.merge(df, self.metadf[['lcz']],
                                  how='left', left_on='name',
                                  right_index=True)

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
    def get_diurnal_statistics(self, obstype='temp', stations=None,
                               startdt=None, enddt=None, plot=True, colorby='name',
                               errorbands=False, verbose=False):
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
            If True, a diurnal plot is made. The default is True.
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
            diurnal_plot(diurnaldf = hourly_avg,
                         errorbandsdf = stddf,
                         title = title,
                         tzstr=tzstring,
                         plot_settings = self.settings.app['plot_settings']['diurnal'],
                         colorby = colorby,
                         lcz_dict = lcz_dict,
                         data_template=self.data_template,
                         obstype=obstype)


        if verbose:
            return hourly_avg, stats

        return hourly_avg

    def get_diurnal_statistics_with_reference(self, obstype='temp', refstation=None, tollerance='30T', stations=None,
                               startdt=None, enddt=None, plot=True, colorby='name', errorbands=False, verbose=False):
        """
        Create an average diurnal cycle for the observation differences of a reference station.

        All observational values are converted to differences with the closest
        (in time) reference observation. No reference observation is found when
        the time difference is larger than the tollerance.

        (In the plot, each station is represed by a line.)

        Parameters
        ----------
        obstype : str, optional
            Element of the metobs_toolkit.observation_types The default is 'temp'.
        refstation : str, optional
            Name of the station to use as a reference. If None, or not in the observations, no reference is
            used the 'get_diurnal_statistics()' is called with the same arguments. The default is None.
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

        """




        obsdf = self.df
        # Check if refstation is a valid station
        if not refstation in obsdf.index.get_level_values('name').unique():
            print(f'WARNING: refstation {refstation} is not found in the dataframe. Continue diurnal statistics without reference.')
            self.get_diurnal_statistic(obstype=obstype,
                                       stations=stations,
                                       startdt=startdt,
                                       enddt=enddt,
                                       colorby=colorby,
                                       plot=plot,
                                       errorbands=errorbands,
                                       verbose=verbose,
                                       data_template=self.data_template)

        # Filter stations
        if not stations is None:
            if isinstance(stations, str):
                stations = [stations]

            if not refstation is None:
                stations.append(refstation)

            obsdf = subset_stations(obsdf, stations)

        # Filter datetimes
        obsdf = datetime_subsetting(df=obsdf,
                                    starttime=startdt,
                                    endtime=enddt)


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
            diurnal_plot(diurnaldf = hourly_avg,
                         errorbandsdf = stddf,
                         title = title,
                         tzstr=tzstring,
                         plot_settings = self.settings.app['plot_settings']['diurnal'],
                         colorby = colorby,
                         lcz_dict = lcz_dict,
                         data_template=self.data_template,
                         obstype = obstype,
                         show_zero_horizontal = True)


        if verbose:
            return hourly_avg, stats, obsdf

        return hourly_avg


    def get_aggregated_diurnal_statistics(self, obstype='temp', stations=None, aggregation=['lcz', 'datetime'], aggregation_method='mean',
                               startdt=None, enddt=None, plot=True, errorbands=False, verbose=False):

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
            List of variables to aggregate to. "datetime" is added to this list if it is not present,
            becaus else there is no time evolution.
        aggregation_method : str, optional
            Which (numpy) function is used to aggregate the observations. The default is 'mean'.
        startdt : datetime.datetime, optional
            The start datetime of the observations to use. If None, all timestamps will be used. The default is None.
        enddt : datetime.datetime, optional
            The end datetime of the observations to use. If None, all timestamps will be used. The default is None.
        plot : bool, optional
            If True, a diurnal plot is made. The default is True.
        errorbands : bool, optional
            If True, the std is representd in the plot by colored bands. The default is False.
        verbose : True, optional
            If True, the dataframse with aggregation information are returned . The default is False.

        Returns
        -------
        tuple (if verbose)
            A tuple of dataframes is returned when verbose is True.

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
        return obsdf
        startdt = obsdf['datetime'].dropna().min()
        enddt = obsdf['datetime'].dropna().max()


        # groupby and take the mean per station per hour.
        stats = obsdf.groupby(groupby_list)[agg_column_name].agg(['mean', 'std', 'median'])

        hourly_avg = stats['mean'].unstack().transpose()

        if plot:

            # generate title
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
            diurnal_plot(diurnaldf = hourly_avg,
                         errorbandsdf = stddf,
                         title = title,
                         tzstr=tzstring,
                         plot_settings = self.settings.app['plot_settings']['diurnal'],
                         colorby = 'name',
                         lcz_dict = None,
                         data_template=self.data_template,
                         obstype=obstype)


        if verbose:
            return hourly_avg, stats

        return hourly_avg

