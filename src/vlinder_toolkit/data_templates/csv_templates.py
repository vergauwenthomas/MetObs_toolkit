#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 08:12:02 2022

@author: thoverga
"""


# %%
# All templates or combined in a list, so if the template is not specified, the corresponding template can be found by iteration.

csv_templates_list = [] #Note that the order of elements is of importance.

# =============================================================================
# VLINDER CSV templates 
# =============================================================================


#templates have nested dict structure where the keys are the column names in the csv file, and the 
# values contain the mapping information to the toolkit classes and names. 



vlinder_brian_csv_template = {
        'Vlinder': {'varname': 'name',
                    'dtype': 'object'},
        
         'Datum': {'varname': '_date',
                   'fmt':'%Y-%m-%d',
                   'dtype': 'object' },
         'Tijd (UTC)': {'varname':'_time',
                        'fmt': '%H:%M:%S',
                        'dtype': 'object',
                        'timezone': 'UTC'},
         'Temperatuur': {'varname':'temp',
                        'units': r'$^o$C',
                        'dtype': 'float64',
                        'description': 'temperature' },
         'Vochtigheid': {'varname':'humidity',
                         'units': '%',
                         'dtype': 'float64',
                         'description': 'relative humidity'},
         'Luchtdruk': {'varname':'pressure', 
                       'units': 'pa',
                       'dtype': 'float64',
                       'description': 'airpressure'},
         'Neerslagintensiteit': {'varname':'precip',
                                 'units': r'l/m$^2$',
                                 'dtype': 'float64',
                                 'description': 'precipitation intensity'},
         'Neerslagsom': {'varname':'precip_sum', 
                         'units': r'l/m^2',
                         'dtype': 'float64',
                         'description': 'precipitation cumulated from midnight'},
         'Windrichting': {'varname': 'wind_direction',
                          'units': r'° from North (CW)',
                          'dtype': 'float64',
                          'description': 'Wind direction'},
         'Windsnelheid': {'varname':'wind_speed',
                          'units': r'm/s',
                          'dtype': 'float64',
                          'description': 'windspeed'},
         'Rukwind': {'varname': 'wind_gust',
                     'units': r'm/s',
                     'dtype': 'float64',
                     'description': 'windgust'},
         'Luchtdruk_Zeeniveau': {'varname': 'pressure_at_sea_level',
                                 'units': 'pa',
                                 'dtype': 'float64',
                                 'description': 'pressure at sea level'},
         'Globe Temperatuur': {'varname': 'radiation_temp',
                               'units': r'celscius denk ik??',
                               'dtype': 'float64',
                               'description': 'Radiative temperature'},
         
        }

csv_templates_list.append(vlinder_brian_csv_template)


vlinder_static_meta_data = {
    'ID': {'varname': '_ID',
                'dtype': 'object'},
    'VLINDER': {'varname': 'name',
                'dtype': 'object'},
    'lat': {'varname': 'lat',
                'dtype': 'float64'},
    'lon': {'varname': 'lon',
                'dtype': 'float64'},
    'stad': {'varname': 'location',
                'dtype': 'object'},
    'benaming': {'varname': 'call_name',
                'dtype': 'object'},
    'Network': {'varname': 'network',
                'dtype': 'object'},
    }

csv_templates_list.append(vlinder_static_meta_data)

        
    
    
    