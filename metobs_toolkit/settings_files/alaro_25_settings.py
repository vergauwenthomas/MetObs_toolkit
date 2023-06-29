#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 27 10:56:38 2023

@author: thoverga
"""

al25_mapinfo = {
    'ALARO_2.5' : {
        "band_of_use":{
                        # Temperatures
                        "temp": {"name": "SFX.T2M", "units": "Celcius"},
                        # "temp_ISBA": {"name": "T2M_ISBA", "units": "Celcius"},
                        # "temp_TEB": {"name": "T2M_TEB", "units": "Celcius"},

                        # Humidity
                        "humidity": {"name": "SFX.HU2M", "units": "percentage"}, #??

                        # Wind
                        "windspeed": {"name": "SFX.W10M", "units": "m/s"}, #??

                        # radiation
                        "net_radiation": {"name": "SFX.RN", "units": "W/m²"}, #??

                        # Fluxes
                        "heat_flux": {"name": "SFX.H", "units": "W/m²"}, #??"
                        "lat_heat_flux": {"name": "SFX.LE", "units": "W/m²"}, #??"
                        "ground_heat_flux": {"name": "SFX.GFLUX", "units": "W/m²"}, #??"






                         },
        'other_mapping': {'datetime' : {'name' : 'date',
                                        'fmt': '%Y-%m-%d %H:%M:%S',
                                        'tz': 'UTC'},
                          'name': {'name': 'name'}},

        'conversions': {'humidity' : 100.0} #multiply by

    }
}
