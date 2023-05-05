#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 15:35:07 2023

@author: thoverga
"""
from datetime import datetime
import pandas as pd


from metobs_toolkit.df_helpers import (init_multiindexdf,
                                        datetime_subsetting,
                                        subset_stations)

class Analysis():
    def __init__(self, obsdf, metadf, settings):
        self.df = obsdf
        self.metadf = metadf
        self.settings = settings



    # =============================================================================
    #     Setters
    # =============================================================================


    def subset_period(self, startdt, enddt):
       if not isinstance(startdt, type(datetime)):
           print(f' {startdt} not a datetime type. Ignore subsetting!')
           return
       if not isinstance(enddt, type(datetime)):
           print(f' {enddt} not a datetime type. Ignore subsetting!')
           return

       self.df = datetime_subsetting(self.df, startdt, enddt)

    # =============================================================================
    #   Helpers
    # =============================================================================
    def _subset_stations(self, stationslist):
        df = self.df.loc[self.df.index.get_level_values(
                    'name').isin(stationslist)]

        present_stations = df.index.get_level_values('name')
        not_present_stations = list(set(stationlist) - set(present_stations))
        if len(not_present_stations)!=0:
            print(f'WARNING: The stations: {not_present_stations} not found in the dataframe.')

        return df


    # =============================================================================
    #   Analyse method
    # =============================================================================

    def get_diurnal_statistics(self, obstype='temp',refstation=None, stations=None, verbose=False,
                               startdt=None, enddt=None):

        obsdf = self.df

        # Filter stations
        if not stations is None:
            if isinstance(stations, str):
                stations = [stations]

            if not refstation is None:
                stations.append(refstation)
            obsdf = subset_stations(obsdf, stations)

        # Check if refstation is a valid station
        if not refstation is None:
            if not refstation in obsdf.index.get_level_values('name').unique():
                print(f'WARNING: refstation {refstation} is not found in the dataframe. Continue diurnal statistics without reference.')
                refstation=None



        # Filter datetimes
        if not startdt is None:
            assert isinstance(startdt, type(datetime)), f'{startdt} is not a datetime.datetime.'
        if not enddt is None:
            assert isinstance(enddt, type(datetime)), f'{enddt} is not a datetime.datetime.'

        obsdf = datetime_subsetting(df=obsdf,
                                    starttime=startdt,
                                    endtime=enddt)

        obsdf = obsdf[obstype].reset_index()
        # Get hours for all records
        obsdf['hour'] = obsdf.reset_index()['datetime'].dt.hour

        if not refstation is None:
            # Create extra identifiers to form unique hours
            obsdf['day'] = obsdf['datetime'].dt.day
            obsdf['month'] = obsdf['datetime'].dt.month
            obsdf['year'] = obsdf['datetime'].dt.year

            # extract refernce from observations
            refdf = obsdf[obsdf['name'] == refstation]
            obsdf = obsdf[obsdf['name']!= refstation]

            # Merge refdf to obsdf using unique hours
            refdf = refdf.rename(columns={obstype: 'ref_'+obstype})
            mergedf = pd.merge(left=obsdf, right=refdf[['ref_'+obstype, 'year', 'month', 'day', 'hour']],
                               how='left', on=['year', 'month', 'day', 'hour'])
            # Compute difference
            agg_column_name = 'difference'
            mergedf[agg_column_name] = mergedf[obstype] - mergedf['ref_'+obstype]

            # overwrite the obsdf
            obsdf = mergedf

        else:
            agg_column_name = obstype #aggregate the measured obstypes


        # groupby and take the mean per station per hour.
        stats = obsdf.groupby(['name', 'hour'])[agg_column_name].agg(['mean', 'std', 'median'])

        hourly_avg = stats['mean'].unstack().transpose()

        if verbose:
            return hourly_avg, stats

        return hourly_avg


# TODO:
    1. dirunal cycle: Make plot function with options to add std as bands, add baseline as black,
    colorby attribute (station/lcz)

    2. add filter to only use ok observations

    3. Setup interaction with Dataset (creation of analysis instance)







    # def make_diurnal_cycle(self, stations=None, obstype='temp',
    #                        relative=False, refstation=None):

    #     # filter to stations and obstype
    #     if stations==None:
    #         stations = self.df.reset_index().name.unique()
    #     stations = set(stations)
    #     stations = list(stations)
    #     for stat in stations:
    #         if not isinstance(stat, str):
    #             print(f'{stat} is not in the station list, therefore the first station is chosen')
    #             return
    #     df = self._subset_stations(stations)

    #     if obstype not in df.columns.values:
    #         print('Input a valid variable for the plotting')
    #         return

    #     # if relative is true we calculate the relative values compared to the reff station
    #     if relative == True:
    #         if isinstance(refstation, str):
    #             if refstation not in stations:
    #                 df_ref=self._subset_stations([refstation])
    #             else:
    #                 print('the refstation need to be different from the other stations')
    #                 return
    #         else:
    #             print('chose one valid reff station')
    #             return
    #         df = df[obstype]
    #         df = df.reset_index()
    #         df_ref = df_ref[obstype]
    #         df_ref = df_ref.reset_index()
    #         df_ref = df_ref.rename(columns={obstype:refstation})
    #         df_ref['datetime'] = pd.to_datetime(df_ref['datetime'])
    #         df_ref = df_ref.set_index('datetime')
    #         df_ref = df_ref.resample('60min').first()
    #         df_ref = df_ref.reset_index()
    #         df_ref.pop('name')
    #         for stat in stations:
    #             df_temp = df[df['name']==stat]
    #             df_temp = df_temp.rename(columns={obstype:stat})
    #             df_temp['datetime'] = pd.to_datetime(df_temp['datetime'])
    #             df_temp = df_temp.set_index('datetime')
    #             df_temp = df_temp.resample('60min').first()
    #             df_temp = df_temp.reset_index()
    #             df_temp.pop('name')
    #             df_ref = df_ref.merge(df_temp,how='left',left_on='datetime',  right_on='datetime')
    #         df_total = df_ref
    #         df_total['hour'] = pd.to_datetime(df_total['datetime'])
    #         df_total['hour'] = df_total['hour'].dt.hour
    #         df_total.pop('datetime')
    #         for stat in stations:
    #             df_total[stat] = df_total[stat]-df_total[refstation]
    #         df_total.pop(refstation)
    #         df_final = pd.DataFrame(columns=df_total.columns, index=range(0,24))
    #         for h in range(0,24):
    #             df_final.loc[h] = df_total[df_total['hour']==h].mean()
    #         df_final=df_final.set_index('hour')
    #         return df_final
    #     #if relative is false we only plot the desired stations
    #     elif relative==False:
    #         df = df[obstype]
    #         df = df.reset_index()
    #         stat_init = stations[0]
    #         print(stat_init)
    #         df_total = df[df['name']==stat_init]
    #         df_total = df_total.rename(columns={obstype:stat_init})
    #         df_total['datetime'] = pd.to_datetime(df_total['datetime'])
    #         df_total = df_total.set_index('datetime')
    #         df_total = df_total.resample('60min').first()
    #         df_total = df_total.reset_index()
    #         df_total.pop('name')
    #         if len(stations)!=1:
    #             for stat in stations[1:]:
    #                 df_temp = df[df['name']==stat]
    #                 df_temp = df_temp.rename(columns={obstype:stat})
    #                 df_temp['datetime'] = pd.to_datetime(df_temp['datetime'])
    #                 df_temp = df_temp.set_index('datetime')
    #                 df_temp = df_temp.resample('60min').first()
    #                 df_temp = df_temp.reset_index()
    #                 df_temp.pop('name')
    #                 df_total = df_total.merge(df_temp,how='left',left_on='datetime',  right_on='datetime')
    #         df_total['hour'] = pd.to_datetime(df_total['datetime'])
    #         df_total['hour'] = df_total['hour'].dt.hour
    #         df_total.pop('datetime')
    #         df_final = pd.DataFrame(columns=df_total.columns, index=range(0,24))
    #         for h in range(0,24):
    #             df_final.loc[h] = df_total[df_total['hour']==h].mean()
    #         df_final=df_final.set_index('hour')
    #         return df_final
    #     else:
    #         print('Input needs to be True or False for relative')
    #         return



