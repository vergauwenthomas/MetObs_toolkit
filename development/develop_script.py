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


lib_folder = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(lib_folder))

# add testdata paths
from tests.push_test.test_data_paths import testdata




tmp_pickle=os.path.join(lib_folder, 'development', 'tmp', 'dev_pickle.pkl')

import metobs_toolkit



#%%
use_dataset = 'demo'
dataset = metobs_toolkit.Dataset()


dataset.update_settings(output_folder=None,
                        input_data_file=testdata[use_dataset]['datafile'],
                        input_metadata_file=testdata[use_dataset]['metadatafile'],
                        data_template_file=testdata[use_dataset]['template'],
                        metadata_template_file=testdata[use_dataset]['template'],
                        )


dataset.import_data_from_file(**testdata[use_dataset]['kwargs'])

dataset.coarsen_time_resolution(freq ='60T')
#%%

outputfolder = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit'


dataset.update_settings(output_folder=outputfolder)

dataset.fairness_coordinates_for_modeldata_csv_creator()







#%%
# import pandas as pd
# import copy


# metadf = dataset.metadf.copy()
# metadf= metadf[metadf['lat'].notna()]
# metadf= metadf[metadf['lon'].notna()]


# bounds = tuple(metadf.total_bounds)


# # add bounds as a column (avoid creating two files with data, and readin in problems in R)
# # savedf.loc[:,'bbox'] = [bounds] * savedf.shape[0]
# metadf['bbox'] = [bounds for _ in range(len(metadf))]


# # reset index so no problems in R
# metadf = metadf.reset_index()


# savedf = metadf[['name', 'lat', 'lon' , 'bbox']]

#%%



