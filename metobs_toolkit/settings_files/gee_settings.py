#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Google earth engine dataset settings

@author: thoverga
"""


gee_datasets = {
    "global_lcz_map": {
        "location": "RUB/RUBCLIM/LCZ/global_lcz_map/v1",  # GEE location
        "usage": "LCZ",  # Human readable application domain
        "band_of_use": "LCZ_Filter",  # band to use for imagecollections (or None if no band available)
        "value_type": "categorical",  # categorical or numeric
        "dynamical": False,  # time evolution? To be used for timeseries
        "scale": 100,
        "is_image": False,
        "is_imagecollection": True,
        "categorical_mapper": {
            1: "Compact highrise",  # mapvalue: (color, human class)
            2: "Compact midrise",
            3: "Compact lowrise",
            4: "Open highrise",
            5: "Open midrise",
            6: "Open lowrise",
            7: "Lightweight lowrise",
            8: "Large lowrise",
            9: "Sparsely built",
            10: "Heavy industry",
            11: "Dense Trees (LCZ A)",
            12: "Scattered Trees (LCZ B)",
            13: "Bush, scrub (LCZ C)",
            14: "Low plants (LCZ D)",
            15: "Bare rock or paved (LCZ E)",
            16: "Bare soil or sand (LCZ F)",
            17: "Water (LCZ G)",
        },
        "credentials": "Demuzere M.; Kittner J.; Martilli A.; Mills, G.; Moede, C.; Stewart, I.D.; van Vliet, J.; Bechtel, B. A global map of local climate zones to support earth system modelling and urban-scale environmental science. Earth System Science Data 2022, 14 Volume 8: 3835-3873. doi:10.5194/essd-14-3835-2022",
    },
    "DEM": {
        "location": "CGIAR/SRTM90_V4",  # GEE location
        "usage": "elevation",  # Human readable application domain
        "band_of_use": "elevation",  # band to use for imagecollections (or None if no band available)
        "value_type": "numeric",  # categorical or numeric
        "dynamical": False,  # time evolution? To be used for timeseries
        "scale": 100,
        "is_image": True,
        "is_imagecollection": False,
        "credentials": "SRTM Digital Elevation Data Version 4",
    },
    "ERA5_hourly": {
        "location": "ECMWF/ERA5_LAND/HOURLY",  # GEE location
        "usage": "ERA5",  # Human readable application domain
        "band_of_use": {"temp": {"name": "temperature_2m", "units": "K"}},
        # band mapper to use for imagecollections (or None if no band available)
        "value_type": "numeric",  # categorical or numeric
        "dynamical": True,  # time evolution? To be used for timeseries
        "scale": 2500,
        "is_image": False,
        "is_imagecollection": True,
        "time_res": "1H",
        "credentials": "",
    },
    "worldcover": {
        "location": "ESA/WorldCover/v200",  # GEE location
        "usage": "landcover",  # Human readable application domain
        "band_of_use": "Map",  # band to use for imagecollections (or None if no band available)
        "value_type": "categorical",  # categorical or numeric
        "dynamical": False,  # time evolution? To be used for timeseries
        "scale": 10,
        "is_image": False,
        "is_imagecollection": True,
        "categorical_mapper": {
            10: "Tree cover",  # mapvalue: (color, human class)
            20: "Shrubland",
            30: "Grassland",
            40: "Cropland",
            50: "Built-up",
            60: "Bare / sparse vegetation",
            70: "Snow and ice",
            80: "Permanent water bodies",
            90: "Herbaceous wetland",
            95: "Mangroves",
            100: "Moss and lichen",
        },
        "aggregation": {
            "water": [70, 80, 90, 95],
            "pervious": [10, 20, 30, 40, 60, 100],
            "impervious": [50],
        },
        "credentials": "https://spdx.org/licenses/CC-BY-4.0.html",
    },
}
