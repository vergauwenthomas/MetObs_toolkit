#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The class object for a Vlinder/mocca station
@author: thoverga
"""



plot_settings = {}
# =============================================================================
# General plot settings
# =============================================================================


# =============================================================================
# Timeseries plots
# =============================================================================
plot_settings['time_series'] = {
    
    #shape
    'figsize' : (10,5),
    
    }
# =============================================================================
# Spatial plot settings
# =============================================================================

plot_settings['spatial_geo'] = {
    #projection
    'proj' : 'Orthographic', #Orthographic or AlbersEqualArea
    'extent' : [ 2.260609, 49.25,  6.118359, 52.350618], #used if observatioons are within
    
    #colors
    # 'cmap' : 'Set1',
    'cmap' : 'inferno_r',
    'n_for_categorical' : 5, #number of quantiles for cat data (not for LCZ) 
    
    #shape
    'figsize': (10,5),
    
    #datetime
    'fmt': "%d/%m/%Y %H:%M:%S UTC" 
    }



print_settings = {
    "fmt_datetime":"%d/%m/%Y %H:%M:%S",
    "max_print_per_line":"40"
    }

# =============================================================================
# variables display strings
# =============================================================================
vars_display = {
    'network': 'network',
    'name': 'station name',
    'call_name': 'pseudo name',
    'location': 'region',
    
    
    'lat': 'latitude',
    'lon': 'longtitude',
    
    'temp': 'temperature',
    'radiation_temp': 'radiation temperature',
    'humidity': 'humidity',
    'precip': 'precipitation intensity',
    'precip_sum': 'cummulated precipitation',
    'wind_speed': 'wind speed',
    'wind_gust': 'wind gust speed',
    'wind_direction': 'wind direction',
    'pressure': 'air pressure',
    'pressure_at_sea_level': 'corrected pressure at sea level',
    
    'lcz':'LCZ'
    
    }