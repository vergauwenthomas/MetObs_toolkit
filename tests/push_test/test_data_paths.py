#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 09:53:55 2023

@author: thoverga
"""

import sys, os

from pathlib import Path

from os.path import join

lib_folder = Path(__file__).resolve().parents[2]
# sys.path.append(str(lib_folder))
# print(str(lib_folder))

sys.path.insert(0, str(lib_folder))
import metobs_toolkit

test_data_dir = os.path.join(str(lib_folder), 'tests', 'test_data')


testdata = {
    # demo
    # 'demo' : {
    #         'datafile': metobs_toolkit.demo_datafile,
    #         'metadatafile':metobs_toolkit.demo_metadatafile,
    #         'template': metobs_toolkit.demo_template,
    #         'kwargs':{},
    #         'coarsen': '20T',
    #         },

    # wide test dataset
    'debug_wide' : {
            'datafile': join(test_data_dir, 'debug_wide.csv'),
            'metadatafile':join(test_data_dir, 'debug_wide_metadata.csv'),
            'template': join(test_data_dir, 'debug_wide_template.csv'),
            'kwargs':{'long_format' : False,
                      'obstype' : 'temp'},
            'coarsen': '20T',
            },

      # Single station dataset
      'single_station' : {
              'datafile': join(test_data_dir, 'single_station.csv'),
              'metadatafile':join(test_data_dir, 'single_station_metadata.csv'),
              'template': join(test_data_dir, 'single_station_template.csv'),
              'kwargs':{'long_format' : False,
                        'obstype' : 'temp'},
              'coarsen': '20T',
              },


    # breaking
    'breaking data' : {
            'datafile': join(test_data_dir, 'testdata_breaking.csv'),
            'metadatafile': None,
            'template': join(test_data_dir, 'template_breaking.csv'),
            'kwargs':{},
            'coarsen': '60T',
            },

    # Kobe congo (single station)
    'Congo_single_station' : {
            'datafile': join(test_data_dir,'testdata_testday', 'Kobe','meteo_soil_clean_2023-01-19.csv'),
            'metadatafile':join(test_data_dir,'testdata_testday', 'Kobe','CONGO_meta.csv'),
            'template': join(test_data_dir,'testdata_testday', 'Kobe','CONGO_template.csv'),
            'kwargs':{},
            'coarsen': '60T',
            },

      # Single Netatmo station Sara
      'single_netatmo_sara_station' : {
              'datafile': join(test_data_dir,'testdata_testday', 'Sara','Outdoor_module_Netatmo_Sara_small.csv'),
              'metadatafile': join(test_data_dir,'testdata_testday', 'Sara','metadata_Outdoor_module_Netatmo_Sara_new.csv'),
              'template': join(test_data_dir,'testdata_testday', 'Sara','template_sara.csv'),
              'kwargs':{'freq_estimation_method' : 'median'},
              'coarsen': '60T',
              },
      # Vlinders 2022
      'vlindergent2022':{
              'datafile': join(test_data_dir,'testdata_testday', 'Sara','Vlinder_gent_2022.csv'),
              'metadatafile': join(test_data_dir,'testdata_testday', 'Sara','all_vlinders_metadata.csv'),
              'template': join(test_data_dir,'testdata_testday', 'Sara','bigvlinder_templatefile.csv'),
              'kwargs':{'freq_estimation_method' : 'median'},
              'coarsen': '60T',

                    }

    }

