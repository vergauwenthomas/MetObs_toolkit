#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:25:02 2022
@author: thoverga
"""

#%%

# import metobs_toolkit
#
import os
import sys
from pathlib import Path
import pandas as pd


lib_folder = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(lib_folder))

# add testdata paths
from tests.push_test.test_data_paths import testdata




tmp_pickle=os.path.join(lib_folder, 'development', 'tmp', 'dev_pickle.pkl')

import metobs_toolkit




#%%

# outfile = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/testdata_testday/Sara/Vlinder_gent_2022.csv'
# infile = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/testdata_testday/Sara/Vlinder_2022.csv'



# def get_small_subset(datafile, nrows=10, **kwargs):
#     return pd.read_csv(datafile, chunksize=nrows+2, **kwargs).get_chunk(nrows)
# test = get_small_subset(infile)


# gentlijst = ['vlinder02', 'vlinder01', 'vlinder05', 'vlinder27']

# df = pd.read_csv(infile)
# subdf = df[df['name'].isin(gentlijst)]


# subdf.to_csv(outfile)




#%%

# # use_dataset = 'debug_wide'
use_dataset = 'demo'

dataset = metobs_toolkit.Dataset()


dataset.update_settings(output_folder=None,
                        input_data_file=testdata[use_dataset]['datafile'],
                        input_metadata_file=testdata[use_dataset]['metadatafile'],
                        data_template_file=testdata[use_dataset]['template'],
                        metadata_template_file=testdata[use_dataset]['template'],
                        )



dataset.import_data_from_file(**testdata[use_dataset]['kwargs'])

# dataset.coarsen_time_resolution(freq = testdata[use_dataset]['coarsen'])



#%%


analysis = dataset.get_analysis()
stats = analysis.get_diurnal_statistics_with_reference(obstype='temp',
                                                       refstation='vlinder01', # define a (rural) reference station of your dataset, insert the name here
                                                       stations=['vlinder02','vlinder27','vlinder28'], # here you can select the stations you want to include, for example: stations=['vlinder01','vlinder02','vlinder25','vlinder27','vlinder28'],
                                                       #if None then all stations are selected
                                                       startdt=None,
                                                       enddt=None,
                                                       plot=True,
                                                       colorby='name',
                                                       errorbands=False, # standard deviation of both reference station and station included
                                                       verbose=False)



#%%

#%%

