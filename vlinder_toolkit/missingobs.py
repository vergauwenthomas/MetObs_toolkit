#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 16:12:41 2023

@author: thoverga
"""



import pandas as pd
import numpy as np
import logging




logger = logging.getLogger(__name__)

# =============================================================================
# Missing observation

# a missing observation is a missing timestamp
# =============================================================================




class Missingob_collection:
    def __init__(self, missing_obs_series):
        self.series = missing_obs_series


    def get_station_missingobs(self, name):
        if name in self.series.index:
            return Missingob_collection(self.series.loc[[name]])
        else:
            # return empty collection
            series = pd.Series(data=[], name='datetime', dtype=object)
            series.index.name = 'name'
            return Missingob_collection(series)


    def remove_missing_from_obs(self, obsdf):

        #Normally there are no missing records in the obsdf
        missing_multiidx = pd.MultiIndex.from_arrays(arrays=[self.series.index.to_list(),
                                                             self.series.to_list()],
                                                  names=[u'name', u'datetime'])

        obsdf = obsdf.drop(index=missing_multiidx, errors='ignore')

        return obsdf


    def get_missing_indx_in_obs_space(self, obsdf, resolutionseries):
        """
        Function to found which missing timestamps are expected in the observation space.
        Because of time coarsening not all missing timestamps are expected in observation space.

        This function handles each station seperatly because stations can have differnent resolution/timerange.



        Parameters
        ----------
        obsdf : pandas.DataFrame()
            Dataset.df.
        resolutionseries : pandas.Series() or Timedelta
            Dataset.metadf['dataset_resolution'].

        Returns
        -------
        missing_obsspace : pandas.MultiIndex
            The multiindex (name - datetime) is returned with the missing timestamps that are expexted in the observation space.

        """

        missing_obsspace = pd.MultiIndex(levels=[['name'],['datetime']],
                                 codes=[[],[]],
                                 names=[u'name', u'datetime'])

        # per stationtion because stations can have different resolutions/timerange
        for sta in self.series.index.unique():

            # Get missing observations in IO space
            sta_missing = self.series.loc[sta]
            if not isinstance(sta_missing, type(pd.Series(dtype=object))):
                sta_missing = pd.Series(data=[sta_missing], index=[sta],
                                        dtype=object)


            # Get start, end and frequency of the observation in obs space
            startdt = obsdf.xs(sta, level='name').index.min()
            enddt = obsdf.xs(sta, level='name').index.max()
            obs_freq = resolutionseries.loc[sta]

            # Make datetimerange
            obsrange = pd.date_range(start=startdt,
                                     end=enddt,
                                     freq = obs_freq,
                                     inclusive="both")

            # # Look which missing timestamps appears obsspace
            sta_missing =sta_missing[sta_missing.isin(obsrange)]


            #Convert to multiindex
            if sta_missing.empty:
                continue

            sta_missing_idx =  pd.MultiIndex.from_arrays(arrays=[[sta]*len(sta_missing),
                                                                 sta_missing.to_numpy()],
                                                      names=[u'name', u'datetime'])

            missing_obsspace = missing_obsspace.append(sta_missing_idx)


        return missing_obsspace