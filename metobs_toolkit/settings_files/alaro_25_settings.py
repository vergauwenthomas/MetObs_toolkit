#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 10:56:38 2023

@author: thoverga
"""

al25_mapinfo = {
    'ALARO_2.5' : {
        "band_of_use": {"temp": {"name": "T2M", "units": "Celcius"},
                        "temp_ISBA": {"name": "T2M_ISBA", "units": "Celcius"},
                        "temp_TEB": {"name": "T2M_TEB", "units": "Celcius"},
                         },
        'other_mapping': {'datetime' : {'name' : 'date',
                                        'fmt': '%Y-%m-%d %H:%M:%S',
                                        'tz': 'UTC'},
                          'name': {'name': 'locationname'}},

    }
}
