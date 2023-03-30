#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:25:02 2022
@author: thoverga
"""

# import vlinder_toolkit
import os, sys
import pandas as pd
import numpy as np
from pathlib import Path


lib_folder = Path(__file__).resolve().parents[1]
sys.path.append(str(lib_folder))
print(str(lib_folder))

import vlinder_toolkit

# % Import

testdatafile = os.path.join(str(lib_folder), 'tests', 'test_data',  'testdata_okt_small.csv')

static_data =  os.path.join(str(lib_folder), 'static_data', 'vlinder_metadata.csv')

# lcz_map = os.path.join(str(lib_folder), 'physiograpy', 'lcz_filter_v1.tif')


#% Setup dataset
settings = vlinder_toolkit.Settings()
settings.update_settings(input_data_file=testdatafile,
                          input_metadata_file=static_data,
                          # geotiff_lcz_file=lcz_map
                          output_folder='/home/thoverga/Documents/VLINDER_github/vlinder_toolkit'
                          )


dataset = vlinder_toolkit.Dataset()
dataset.import_data_from_file(coarsen_timeres=False)



# dataset.apply_quality_control()

# test = dataset.get_qc_stats()

# dataset.write_to_csv(filename='remove_me', add_final_labels=True)





#%%

obs = {'temp': {'map_column': 'Temperatuur', 'map_unit': 'Celcius', 'map_desc': '2mT'}, 'radiation_temp': {'map_column': np.nan, 'map_unit': np.nan, 'map_desc': np.nan}, 'humidity': {'map_column': 'Vochtigheid', 'map_unit': '%', 'map_desc': 'Relative humidity'}, 'precip': {'map_column': np.nan, 'map_unit': np.nan, 'map_desc': np.nan}, 'precip_sum': {'map_column': np.nan, 'map_unit': np.nan, 'map_desc': np.nan}, 'wind_speed': {'map_column': np.nan, 'map_unit': np.nan, 'map_desc': np.nan}, 'wind_gust': {'map_column': np.nan, 'map_unit': np.nan, 'map_desc': np.nan}, 'wind_direction': {'map_column': np.nan, 'map_unit': np.nan, 'map_desc': np.nan}, 'pressure': {'map_column': np.nan, 'map_unit': np.nan, 'map_desc': np.nan}, 'pressure_at_sea_level': {'map_column': np.nan, 'map_unit': np.nan, 'map_desc': np.nan}}

dt = {'datetime': {'map_column': np.nan, 'map_fmt': np.nan}, 'date': {'map_column': 'Datum', 'map_fmt': '%Y-%m-%d'}, 'time': {'map_column': 'Tijd (UTC)', 'map_fmt': '%H:%M:%S'}}

meta = {'name': {'map_column': 'Vlinder'}, 'lat': {'map_column': 'Windrichting'}, 'lon': {'map_column': 'Rukwind'}, 'location': {'map_column': np.nan}, 'call_name': {'map_column': np.nan}, 'network': {'map_column': np.nan}}

#%%

df = pd.DataFrame(obs).transpose()
df.index.name = 'varname'
df = df.dropna(axis = 0, how = 'all')
df = df.reset_index()
df = df.rename(columns={'map_column': 'template column name',
                        'map_unit': 'units',
                        'map_desc': 'description'})
df['dtype'] = 'float64'






