#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 15:35:07 2023

@author: thoverga
"""
from datetime import datetime

# from vlinder_toolkit.dataset import Dataset
from vlinder_toolkit.df_helpers import (init_multiindexdf,
                                        datetime_subsetting)

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
        return df


    # =============================================================================
    #   Analyse method
    # =============================================================================
    def make_diurnal_cycle(self, stations=None, obstype='temp',
                           relative=False, refstation=None):

        # filter to stations and obstype
        if isinstance(stations, str):
            df = self._subset_stations([stations])

