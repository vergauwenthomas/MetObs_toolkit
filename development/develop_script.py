#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:25:02 2022
@author: thoverga
"""

#%%

# import metobs_toolkit

import os
import sys
from pathlib import Path


lib_folder = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(lib_folder))

import metobs_toolkit


#%%


file = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/testdata_testday/Amsterdam_D2222z6together.csv'
metafile ='/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/testdata_testday/Latlon_stations_Amsterdam.csv'



# metobs_toolkit.build_template_prompt()


#%%
# Make an empty dataset
# dataset = metobs_toolkit.Dataset()

# # Add the demo data files to the dataset settings
# dataset.update_settings(input_data_file = metobs_toolkit.demo_datafile,
#                         input_metadata_file = metobs_toolkit.demo_metadatafile,
#                         data_template_file = metobs_toolkit.demo_template,
#                         metadata_template_file = metobs_toolkit.demo_template, # Contains also the metadata mapping
#                         output_folder = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/development'
#                         )

# dataset.import_data_from_file()



#%%

import pandas as pd



data = pd.read_csv(file)
columnnames = data.columns.to_list()


metadata = pd.read_csv(metafile)
metacolumnnames = metadata.columns.to_list()



data_test = data.iloc[0].to_dict()


template_dict = {'datetime': {'orig_name': 'DateTime', 'format': '%Y-%m-%d %H:%M:%S'}, 'temp': {'units': 'Celcius', 'description': 'feiemoij'}}
metatemplate_dict = {'name': {'orig_name': 'Station'}, 'lat': {'orig_name': 'Latitude'}, 'lon': {'orig_name': 'Longitude'}}


stationnames = [
 'D2194',
 'D2195',
 'D2198',
 'D2199',
 'D2221',
 'D2222',
 'D2223',
 'D2225',
 'D2226',
 'D2227',
 'D2228',
 'D2229',
 'D2230',
 'D2231',
 'D2235',
 'D2236',
 'D2237',
 'D2239',
 'D2240',
 'D2241',
 'D2244',
 'D2245',
 'D2246',
 'D2247']


