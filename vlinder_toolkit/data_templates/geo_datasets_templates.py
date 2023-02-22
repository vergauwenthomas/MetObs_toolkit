#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 11:40:47 2022

@author: thoverga
"""

geo_datasets = []

lcz_dataset = {
    "usage": 'LCZ',
    "format": 'geotiff',
    "epsg": "epsg:4326",
    "version": 'v1 - 2022, lcz for 2018',
    "reference": "Matthias Demuzere, Jonas Kittner, Alberto Martilli, Gerald Mills, Christian Moede, Iain D. Stewart, Jasper van Vliet, & Benjamin Bechtel. (2022). Global map of Local Climate Zones (1.0.0) [Data set]. Zenodo. https://doi.org/10.5281/zenodo.6364594",
    "covers": {
        1:{'cover_name': 'LCZ-1',
           'class_description': 'Compact highrise',
           'hex_color': '#910613'},
        2:{'cover_name': 'LCZ-2',
           'class_description': 'Compact midrise',
           'hex_color': '#D9081C'},
        3:{'cover_name': 'LCZ-3',
           'class_description': 'Compact lowrise',
           'hex_color': '#FF0A22'},
        4:{'cover_name': 'LCZ-4',
           'class_description': 'Open highrise',
           'hex_color': '#C54F1E'},
        5:{'cover_name': 'LCZ-5',
           'class_description': 'Open midrise',
           'hex_color': '#FF6628'},
        6:{'cover_name': 'LCZ-6',
           'class_description': 'Open lowrise',
           'hex_color': '#FF985E'},
        7:{'cover_name': 'LCZ-7',
           'class_description': 'Lightweight low-rise',
           'hex_color': '#FDED3F'},
        8:{'cover_name': 'LCZ-8',
           'class_description': 'Large lowrise',
           'hex_color': '#BBBBBB'},
        9:{'cover_name': 'LCZ-9',
           'class_description': 'Sparsely built',
           'hex_color': '#FFCBAB'},
        10:{'cover_name': 'LCZ-10',
           'class_description': 'Heavy Industry',
           'hex_color': '#565656'},
        11:{'cover_name': 'LCZ-11 (A)',
           'class_description': 'Dense trees',
           'hex_color': '#006A18'},
        12:{'cover_name': 'LCZ-12 (B)',
           'class_description': 'Scattered trees',
           'hex_color': '#00A926'},
        13:{'cover_name': 'LCZ-13 (C)',
           'class_description': 'Bush, scrub',
           'hex_color': '#628432'},
        14:{'cover_name': 'LCZ-14 (D)',
           'class_description': 'Low plants',
           'hex_color': '#B5DA7F'},
        15:{'cover_name': 'LCZ-15 (E)',
           'class_description': 'Bare rock or paved',
           'hex_color': '#000000'},
        16:{'cover_name': 'LCZ-16 (F)',
           'class_description': 'Bare soil or sand',
           'hex_color': '#FCF7B1'},
        17:{'cover_name': 'LCZ-17 (G)',
           'class_description': 'Water',
           'hex_color': '#656BFA'}
        }
    
    }


geo_datasets.append(lcz_dataset)
