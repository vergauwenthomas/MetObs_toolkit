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

dataset = metobs_toolkit.Dataset()



#%%

dataset.update_settings(output_folder=None,
                        input_data_file=metobs_toolkit.demo_datafile,
                        # input_data_file = data_file,
                        input_metadata_file=metobs_toolkit.demo_metadatafile,
                        data_template_file=metobs_toolkit.demo_template,
                        # data_template_file = template_file,
                        metadata_template_file=metobs_toolkit.demo_template,
                        )


dataset.import_data_from_file()
dataset.coarsen_time_resolution()

dataset.apply_quality_control()

dataset.update_gaps_and_missing_from_outliers()

#%%
dataset.get_missing_obs_info()
dataset.fill_missing_obs_linear()


#%%







#%% test




# model_data.convert_units_to_tlk(obstype='temp',
#                                 target_unit_name='Celcius',
#                                 conv_expr=None)

# print(model_data)
# model_data.convert_units_to_tlk(obstype='temp',
#                                 target_unit_name='hot_Celcius',
#                                 conv_expr="x + 10")

# model_data.get_ERA5_data(metadf=dataset.metadf, startdt=datetime(2022, 9, 1), enddt=datetime(2022, 9, 2))

# print(model_data)
# model_data.add_gee_dataset(mapname='worldcover',
#                            gee_location='location/loc/dfmijfe',
#                            obstype='temp',
#                            bandname='temp 2m passive field',
#                            units ='C',
#                            scale = 100,
#                            time_res='1H',
#                            is_image=False,
#                            is_numeric=True,
#                            credentials='bladiblie')


# model_data.list_gee_datasets()


# model_data.set_model_from_csv('/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/era5_data.csv')










