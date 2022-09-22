#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The class object for a Vlinder/mocca station
@author: thoverga
"""

import pandas as pd

class Station:
    def __init__(self, station_name, network_name):
        self.network = network_name
        self.name = station_name
        
        #Meta data without processing
        self.lat = []
        self.lon = []
        self.location = None #ex. Boerenkreek
        
        
        #Observations
        self.temp = pd.Series()
        self.radiation_temp = pd.Series() 
        
        self.humidity = pd.Series()
        
        self.precip = pd.Series()
        self.precip_sum = pd.Series()
        
        self.wind_speed = pd.Series()
        self.wind_gust = pd.Series()
        self.wind_direction = pd.Series()
        
        self.pressure = pd.Series()
        self.pressure_at_sea_level = pd.Series()
        
        
    



