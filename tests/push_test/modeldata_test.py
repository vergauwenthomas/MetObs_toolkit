#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 20:32:45 2023

@author: thoverga
"""


import sys, os

from pathlib import Path


lib_folder = Path(__file__).resolve().parents[2]

import metobs_toolkit
# print(metobs_toolkit.__version__)

#%% Import dataset

dataset = metobs_toolkit.Dataset()

dataset.update_settings(output_folder=None,
                        input_data_file=metobs_toolkit.demo_datafile,
                        template_file=metobs_toolkit.demo_template,
                        input_metadata_file=metobs_toolkit.demo_metadatafile
                        )


dataset.import_data_from_file()
dataset.coarsen_time_resolution()

# dataset.get_modeldata()

#%% test adding gee information
model_data = metobs_toolkit.Modeldata("ERA5_hourly")
model_data.add_band_to_gee_dataset(bandname='surface_pressure',
                                   obstype='pressure',
                                   units='pa')

model_data.add_gee_dataset(mapname='new dataset name',
                            gee_location='location/loc/dfmijfe',
                            obstype='temp',
                            bandname='temp 2m passive field',
                            units ='C',
                            scale = 100,
                            time_res='1H',
                            is_image=False,
                            is_numeric=True,
                            credentials='bladiblie')
#%% Import modeldata
model_data = metobs_toolkit.Modeldata("ERA5_hourly")

csv_file = os.path.join(lib_folder, 'tests', 'test_data', 'era5_modeldata_test.csv')

model_data.set_model_from_csv(csv_file)

#%% Test repr

print(model_data)

#%% test saving and importing
outfolder = os.path.join(lib_folder, 'tests', 'test_data')
pkl_file = 'delete_me_if_you_see_me'
# save
model_data.save_modeldata(outputfolder=outfolder, filename=pkl_file)

# read it again
newmod = metobs_toolkit.Modeldata('ERA5_hourly')
newmod2 = newmod.import_dataset(folder_path=outfolder, filename=pkl_file+'.pkl')

# delete file
fullpath = os.path.join(outfolder, pkl_file+'.pkl')
if os.path.exists(fullpath):
    os.remove(fullpath)


#%% test interpolation
interpdf = model_data.interpolate_modeldata(dataset.df.index)

assert interpdf[interpdf['temp'].isnull()].shape == (28, 1), 'Error in modeldata interpolation'




#%% Test plotting

a = model_data.df.shape

model_data.make_plot(stationnames=['vlinder01', 'vlinder02'])


assert model_data.df.shape == (10052, 1), 'Shape of modeldata df changed after plotting.'


model_data.make_plot(dataset=dataset, show_outliers=False)

