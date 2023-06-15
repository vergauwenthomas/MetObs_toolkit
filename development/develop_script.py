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

tmp_pickle=os.path.join(lib_folder, 'development', 'tmp', 'dev_pickle.pkl')

import metobs_toolkit


#%%

# # data
# era5_congo_file = '/home/thoverga/Downloads/era5_data_kongo.csv'
# data_file = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/testdata_testday/Kobe/meteo_soil_clean_2023-01-19.csv'
# metadata_file = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/testdata_testday/Kobe/CONGO_meta.csv'
# template_file ='/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/testdata_testday/Kobe/CONGO_template.csv'
#%%

dataset = metobs_toolkit.Dataset()

dataset.update_settings(output_folder=None,
                        input_data_file=metobs_toolkit.demo_datafile,
                        input_metadata_file=metobs_toolkit.demo_metadatafile,
                        data_template_file=metobs_toolkit.demo_template,
                        metadata_template_file=metobs_toolkit.demo_template)


dataset.import_data_from_file()

dataset.coarsen_time_resolution()

dataset.apply_quality_control()


dataset.update_gaps_and_missing_from_outliers(n_gapsize=10)



#%%


erafile = "/home/thoverga/Downloads/era5_data_3.csv"

model = metobs_toolkit.Modeldata('era')
# dataset.get_modeldata()


model.set_model_from_csv(erafile)



#%%
for gap in dataset.gaps:
    print(gap)


#%%

sta = dataset.get_station('vlinder01')
sta.gaps[0].get_info()
sta.fill_gaps_automatic(model, max_interpolate_duration_str='10H', overwrite_fill=True)



#%%
# sta.fill_gaps_linear()
# sta.fill_gaps_era5(modeldata=model)

# stanames = dataset.df.index.get_level_values('name').unique().to_list()

# for stan in stanames:
#     sta = dataset.get_station(stan)
#     print(f' {stan} has these gaps: {sta.gaps}')
#     sta.fill_gaps_automatic(model, max_interpolate_duration_str='10H')



#%%





