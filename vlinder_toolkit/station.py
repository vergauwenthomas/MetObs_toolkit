#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 24 08:42:56 2022

@author: thoverga
"""

import pandas as pd
import logging

from .plotting_functions import timeseries_plot
from .settings import Settings
logger = logging.getLogger(__name__)


#%%
class Station:
    def __init__(self, name, df, outliersdf, gapsdf, meta_series, data_template):
        self.name = name
        self.df = df
        self.outliersdf = outliersdf
        self.gapsdf = gapsdf
        self.meta_series = meta_series
        self.data_template = data_template
        
        
    def make_plot(self, variable='temp', title=None):
         
        """
        Make a timeseries plot of one attribute.
        Parameters
        ----------
        variable : str, optional
            Name of attribute to plot. Must be one of [temp, radiation_temp, humidity, precip, wind_speed wind_gust, wind_direction, pressure, pressure_at_sea_level].
            The default is 'temp'.
        **kwargs : 
            named-arguments that are passed to matplolib.pyplot.plot()
        Returns
        -------
        ax : AxesSubplot
            AxesSubplot is returned so layer can be added to it.
        """
        logger.info(f'Make {variable} plot for Station {self.name}.')
        #Load default plot settings
        default_settings=Settings.plot_settings['time_series']
      
        
        #Make title
        if isinstance(title, type(None)):
            title = self.name + ': ' + self.data_template.loc['description',variable]
    
        
        #make figure
        ax = timeseries_plot(dtseries=self.df[variable],
                             title=title,
                             xlabel='',
                             ylabel=self.data_template.loc['units',variable],
                             figsize=default_settings['figsize'],
                             )
                        
        
        
        return ax
        
        
        
    