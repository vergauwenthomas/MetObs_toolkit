#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module contains all relevant data for each observation type.

Created on Fri Sep 23 09:48:46 2022

@author: thoverga
"""

datetime_settings = {
    'input_date_format': '%Y-%m-%d',
    'input_time_format': '%H:%M:%S',
    
    'string_representation_format': '%Y-%m-%d %H:%M:%S',
    
    }


observational_info = {
     '_date': {'units': 'days', 'format': 'object'},
     '_time': {'units': 'seconds', 'format': 'object'},
     'temp': {'units': r'$^o$C', 'format': 'float64', 'description': 'temperature'},
     'humidity': {'units': '%', 'format': 'float64', 'description': 'relative humidity'},
     'pressure': {'units': 'pa', 'format': 'float64', 'description': 'airpressure'},
     'precip': {'units': r'l/m$^2$', 'format': 'float64', 'description': 'precipitation intensity'},
     'precip_sum': {'units': r'l/m^2', 'format': 'float64', 'description': 'precipitation cumulated from midnight'},
     'wind_direction': {'units': r'° from North (CW)', 'format': 'float64', 'description': 'Wind direction'},
     'wind_speed': {'units': r'm/s', 'format': 'float64', 'description': 'windspeed'},
     'wind_gust': {'units': r'm/s', 'format': 'float64', 'description': 'windgust'},
     'pressure_at_sea_level': {'units': 'pa', 'format': 'float64', 'description': 'pressure at sea level'},
     'radiation_temp': {'units': r'??', 'format': 'float64', 'description': 'Radiative temperature'},
     'name': {'format': 'object'},
    }


#%% restructure_func
def compress_dict(valuesname):
    returndict = {}
    for key, item in observational_info.items():
        if valuesname in item:
            returndict[key] = item[valuesname]
    return returndict



#%% Derived variables
dtypedict = compress_dict('format')
description_dict = compress_dict('description')
unit_dict = compress_dict('units')




read_datetime_format = datetime_settings['input_date_format'] + ' ' + datetime_settings['input_time_format']