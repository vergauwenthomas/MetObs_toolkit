#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:25:02 2022
@author: thoverga
"""

#%%
import vlinder_toolkit
import os
import sys
import pandas as pd
import numpy as np
from pathlib import Path


lib_folder = Path(__file__).resolve().parents[1]
sys.path.append(str(lib_folder))




#%% % Import


testdatafile = os.path.join(
    str(lib_folder), 'tests', 'test_data',  'testdata_okt_small.csv')

static_data = os.path.join(
    str(lib_folder), 'static_data', 'vlinder_metadata.csv')

#%%
dataset = vlinder_toolkit.Dataset()


dataset.update_settings(input_data_file=testdatafile,
                        input_metadata_file=static_data,
                        output_folder='/home/thoverga/Documents/VLINDER_github/vlinder_toolkit'
                        )


from datetime import datetime
startdt = datetime(2020, 8,1)
enddt = datetime(2020,8,2)
test = dataset.import_data_from_database(start_datetime=startdt, end_datetime=enddt)


# dataset.import_data_from_file(coarsen_timeres=True)

# dataset.apply_quality_control()
# dataset.get_qc_stats()
# dataset.make_geo_plot()
# dataset.make_plot()


#%%
print('done')

