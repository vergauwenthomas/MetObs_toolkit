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


# file = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/testdata_testday/Amsterdam_D2222z6together.csv'
# metafile ='/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/testdata_testday/Latlon_stations_Amsterdam.csv'



# # metobs_toolkit.build_template_prompt()


#%%
# Make an empty dataset
dataset = metobs_toolkit.Dataset()

# Add the demo data files to the dataset settings
dataset.update_settings(input_data_file = metobs_toolkit.demo_datafile,
                        input_metadata_file = metobs_toolkit.demo_metadatafile,
                        data_template_file = metobs_toolkit.demo_template,
                        metadata_template_file = metobs_toolkit.demo_template, # Contains also the metadata mapping
                        output_folder = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/development'
                        )

dataset.import_data_from_file()



#%%

dataset.get_lcz()
# sta = dataset.get_station('vlinder01')
# sta.make_plot()

dataset.make_geo_plot(obstype='temp')

#%%


analy = dataset.get_analysis()
analy.get_aggregated_diurnal_statistics( title='hhoi', legend=False)
