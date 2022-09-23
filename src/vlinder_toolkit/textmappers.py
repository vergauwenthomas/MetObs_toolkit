#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This script contians standardized mapper dicts that are used in IO. 


Created on Fri Sep 23 09:30:38 2022

@author: thoverga
"""

download_cols_to_class_cols_map = {
     'Datum': '_date',
     'Tijd (UTC)': '_time',
     'Temperatuur': 'temp',
     'Vochtigheid': 'humidity',
     'Luchtdruk': 'pressure',
     'Neerslagintensiteit': 'precip',
     'Neerslagsom': 'precip_sum',
     'Windrichting': 'wind_direction',
     'Windsnelheid': 'wind_speed',
     'Rukwind': 'wind_gust',
     'Luchtdruk_Zeeniveau': 'pressure_at_sea_level',
     'Globe Temperatuur': 'radiation_temp',
     'Vlinder': 'name'
    }

