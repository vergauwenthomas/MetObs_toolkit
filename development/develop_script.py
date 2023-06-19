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

# single station
file = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/single_station.csv'
metadfile = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/single_station_metadata.csv'


# wide
file = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/debug_wide.csv'
metadfile = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/debug_wide_metadata.csv'



#metobs_toolkit.build_template_prompt(debug=True)

from datetime import datetime

demo_metadata="/home/mivieijra/Documents/CLIMPACTH/toolkit/vlinder_toolkit/metobs_toolkit/datafiles/demo_metadatafile.csv"
demo_data="/home/mivieijra/Documents/CLIMPACTH/toolkit/vlinder_toolkit/metobs_toolkit/datafiles/demo_datafile.csv"
demo_template="/home/mivieijra/Documents/CLIMPACTH/toolkit/vlinder_toolkit/metobs_toolkit/datafiles/demo_templatefile.csv"

dataset = metobs_toolkit.Dataset()
dataset.update_settings(input_data_file=demo_data,  
    input_metadata_file=demo_metadata,
    data_template_file=demo_template,
    metadata_template_file=demo_template)

dataset.import_data_from_file()


#%%


model_data = metobs_toolkit.Modeldata("ERA5")
model_data.get_ERA5_data(metadf=dataset.metadf, startdt=datetime(2022, 9, 1), enddt=datetime(2022, 9, 2))


#%%
test = model_data.make_plot(dataset=dataset, starttime=datetime(2022, 9, 1), endtime=datetime(2022, 9, 2))
