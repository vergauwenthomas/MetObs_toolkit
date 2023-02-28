#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 15:59:40 2023

@author: thoverga
"""


import os, sys
from pathlib import Path
toolkit_folder = Path(__file__).resolve().parents[2]
sys.path.append(str(toolkit_folder)) #turn of in operation mode


import vlinder_toolkit

# =============================================================================
# General settings
# =============================================================================



# =============================================================================
# IO page settings 
# =============================================================================

meta_names = {'name':{'dtype':'object', 'description': 'ID or name of a station.'},
              'lat':{'dtype':'float64', 'description': 'Latitude of the station'},
              'lon':{'dtype':'float64', 'description': 'Longitude of the station.'},
              'location':{'dtype':'object', 'description': 'City or region of this station.'},
              'call_name':{'dtype':'object', 'description': 'Nickname of the station'},
              'network': {'dtype':'object', 'description': 'Network of the station'}}



dt_names = {'datetime':{'dtype':'object', 'format':'%Y-%m-%d %H:%M:%S', 'description': 'Map only af a datetime column is available.'},
            '_date':{'dtype':'object', 'format':'%Y-%m-%d', 'description': 'Map only af a date column is available.'},
            '_time':{'dtype':'object', 'format':'%H:%M:%S', 'description': 'Map only af a time column is available.'}}



obs_names = {'temp':{'dtype':'float64', 'description': '2m-temperature', 'units':['Celcius', 'K']},
             'humidity':{'dtype':'float64', 'description': 'Relative humidity', 'units':['%']},
             'pressure':{'dtype':'float64', 'description': 'Observed atmospheric pressure', 'units':['pa']},
             'precip':{'dtype':'float64', 'description': 'Precipitation intensity', 'units':['l/m²']},
             #TODO aanvullen
             }