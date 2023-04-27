#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  4 14:03:58 2022

@author: thoverga
"""

vlinder_metadata_db_template = {
    "VLINDER": {"varname": "name", "dtype": "object"},
    "ID": {"varname": "id", "dtype": "object"},  # for merging
    "Location": {"varname": "call_name", "dtype": "object"},
    "stad": {"varname": "location", "dtype": "object"},
    "Latitude": {"varname": "lat", "dtype": "float"},
    "Longitude": {"varname": "lon", "dtype": "float"},
}


vlinder_observations_db_template = {
    "StationID": {"varname": "id", "dtype": "object"},  # for merging
    "datetime": {
        "varname": "datetime",
        "fmt": "%Y-%m-%d %H:%M:%S",
        "dtype": "object",
        "timezone": "UTC",
    },
    "temperature": {
        "varname": "temp",
        "units": r"$^o$C",
        "dtype": "float64",
        "description": "temperature",
    },
    "humidity": {
        "varname": "humidity",
        "units": "%",
        "dtype": "float64",
        "description": "relative humidity",
    },
    "pressure": {
        "varname": "pressure",
        "units": "pa",
        "dtype": "float64",
        "description": "airpressure",
    },
    "RainIntensity": {
        "varname": "precip",
        "units": r"l/m$^2 per ?? tijdseenheid$",
        "dtype": "float64",
        "description": "precipitation intensity",
    },
    "RainVolume": {
        "varname": "precip_sum",
        "units": r"l/m^2",
        "dtype": "float64",
        "description": "precipitation cumulated from midnight",
    },
    "WindDirection": {
        "varname": "wind_direction",
        "units": r"° from North (CW)",
        "dtype": "float64",
        "description": "Wind direction",
    },
    "WindSpeed": {
        "varname": "wind_speed",
        "units": r"m/s",
        "dtype": "float64",
        "description": "windspeed",
    },
    "WindGust": {
        "varname": "wind_gust",
        "units": r"m/s",
        "dtype": "float64",
        "description": "windgust",
    },
    "pressure_0": {
        "varname": "pressure_at_sea_level",
        "units": "pa",
        "dtype": "float64",
        "description": "pressure at sea level",
    },
    "BlackGlobeTemp": {
        "varname": "radiation_temp",
        "units": r"celscius denk ik??",
        "dtype": "float64",
        "description": "Radiative temperature",
    },
}
