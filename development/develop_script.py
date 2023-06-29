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
use_dataset = 'vlindergent2022'

dataset = metobs_toolkit.Dataset()


dataset.update_settings(output_folder=None,
                        input_data_file=testdata[use_dataset]['datafile'],
                        input_metadata_file=testdata[use_dataset]['metadatafile'],
                        data_template_file=testdata[use_dataset]['template'],
                        metadata_template_file=testdata[use_dataset]['template'],
                        )
tz ='UTC'

dataset.update_timezone(timezonestr=tz)


dataset.import_data_from_file(**testdata[use_dataset]['kwargs'])

dataset.coarsen_time_resolution(freq = testdata[use_dataset]['coarsen'])




#%%



ann= dataset.get_analysis()

plotdf = ann.get_anual_statistics(groupby=['name', 'hour'])




#%%

