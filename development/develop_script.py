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
import pandas as pd
import datetime


# Make an empty dataset
dataset = metobs_toolkit.Dataset()

# Add the demo data files to the dataset settings
dataset.update_settings(input_data_file = metobs_toolkit.demo_datafile,
                        input_metadata_file = metobs_toolkit.demo_metadatafile,
                        data_template_file = metobs_toolkit.demo_template,
                        metadata_template_file = metobs_toolkit.demo_template, # Contains also the metadata mapping
                        output_folder='/home/thoverga/Documents/VLINDER_github/MetObs_toolkitss'
                        )



# Load the data from the demo data files
dataset.import_data_from_file()

# dataset.coarsen_time_resolution()
dataset.get_altitude()
dataset.get_landcover()
dataset.get_lcz()



dataset.save_dataset(outputfolder='/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data',
                     filename='tests_dataset')



#%%
# outputfolder = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkitss'
dataset.save_dataset()

anal = dataset.get_analysis()


# anal = anal.apply_filter('lcz == "Compact midrise"')

anal.get_lc_correlation_matrices(obstype=['temp', 'humidity'], groupby_labels=['hour'])


anal.plot_correlation_variation()




#%%





