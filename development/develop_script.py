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
# data Ian
# datafile = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/WBGTdata_totaal.csv'

# data Wout
datafile = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/Outdoor_module_Netatmo_Sara_new.csv'
metafile ='/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/metadata_Outdoor_module_Netatmo_Sara_new.csv'



# metobs_toolkit.build_template_prompt()




# Make an empty dataset
dataset = metobs_toolkit.Dataset()

# Add the demo data files to the dataset settings
dataset.update_settings(input_data_file = metobs_toolkit.demo_datafile,
                        input_metadata_file = metobs_toolkit.demo_metadatafile,
                        data_template_file = metobs_toolkit.demo_template,
                        metadata_template_file = metobs_toolkit.demo_template, # Contains also the metadata mapping
                        output_folder = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/development'
                        )

