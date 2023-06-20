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




# # data
# era5_congo_file = '/home/thoverga/Downloads/era5_data_kongo.csv'
data_file = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/testdata_testday/Ian/ATHTS01_all.csv'
# metadata_file = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/testdata_testday/Kobe/CONGO_meta.csv'
template_file ='/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/testdata_testday/Ian/template.csv'




#%%

# dataset = metobs_toolkit.Dataset()

# dataset.update_settings(output_folder=None,
#                         input_data_file=metobs_toolkit.demo_datafile,
#                         # input_data_file = data_file,
#                         # input_metadata_file=metobs_toolkit.demo_metadatafile,
#                         data_template_file=metobs_toolkit.demo_template,
#                         # data_template_file = template_file,
#                         # metadata_template_file=metobs_toolkit.demo_template,
#                         )


# dataset.import_data_from_file()
# dataset.coarsen_time_resolution()




#%%


model_data = metobs_toolkit.Modeldata("ERA5")
# model_data.get_ERA5_data(metadf=dataset.metadf, startdt=datetime(2022, 9, 1), enddt=datetime(2022, 9, 2))

model_data.add_band_to_gee_dataset(mapname='ERA5_hourly',
                                   bandname = 'hophop',
                                   obstype='temp',
                                   units='idk',
                                   overwrite=True)



# model_data.list_gee_datasets()


# model_data.set_model_from_csv('/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/era5_data.csv')










