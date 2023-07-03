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
import time
import math


lib_folder = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(lib_folder))

# add testdata paths
from tests.push_test.test_data_paths import testdata




tmp_pickle=os.path.join(lib_folder, 'development', 'tmp', 'dev_pickle.pkl')

import metobs_toolkit


#%%


datafile ='/home/thoverga/Downloads/sub_set.csv'

template_file = '/home/thoverga/Downloads/sub_set_template.csv'

dataset = metobs_toolkit.Dataset()

dataset.update_settings(input_data_file=datafile,
                        data_template_file=template_file)

dataset.import_data_from_file()
# print(dataset)

#%%
# # use_dataset = 'debug_wide'
# use_dataset = 'single_netatmo_sara_station'
use_dataset = 'vlindergent2022'
# use_dataset = 'demo'
dataset = metobs_toolkit.Dataset()


dataset.update_settings(output_folder=None,
                        input_data_file=testdata[use_dataset]['datafile'],
                        input_metadata_file=testdata[use_dataset]['metadatafile'],
                        data_template_file=testdata[use_dataset]['template'],
                        metadata_template_file=testdata[use_dataset]['template'],
                        )



dataset.import_data_from_file(**testdata[use_dataset]['kwargs'])

dataset.coarsen_time_resolution(freq = testdata[use_dataset]['coarsen'])
# dataset.apply_quality_control()
# dataset.get_lcz()
# dataset.update_gaps_and_missing_from_outliers()
#%%

ann = dataset.get_analysis()

#%%

from datetime import datetime
# Compute mean annual cycle for each station + plot
# Create diurnal cycle for the meteorological summer of 2022 (June, July, August)
stats = ann.get_diurnal_statistics(obstype='temp', # here you can change the varible for which you want to plot the diurnal cycle
                                    stations=None, # here you can select the stations you want to include, for example: stations=['vlinder01','vlinder02','vlinder25','vlinder27','vlinder28']
                                    startdt= datetime(2022,6,1), # here you can change the start date and time
                                    enddt= datetime(2022,8,31), # here you can change the end date and time
                                    plot=True, # create immediatly a plot, if false, then no plot is created
                                    colorby='lcz', # here you can change the color of the lines in the graph
                                    errorbands=False, # when you set this to True, then error bands are created around the curves based on the standard deviation
                                    verbose=False) #if True, an extra dataframe with the std is returned aswell.




# ann.get_diurnal_statistics(colorby='lcz',
#                             errorbands=True)


# ann.get_diurnal_statistics_with_reference(refstation='vlinder01')

# stats = ann.get_aggregated_diurnal_statistics(obstype='temp', # here you can change the varible for which you want to plot the diurnal cycle
#                                               stations=None, # here you can select the stations you want to include, for example: stations=['vlinder01','vlinder02','vlinder25','vlinder27','vlinder28']
#                                               aggregation=['name','hour'],
#                                               aggregation_method='mean',
#                                               horizontal_axis='season',
#                                               startdt=None, # here you can change the start date and time
#                                               enddt=None, # here you can change the end date and time
#                                               plot=True, # create immediatly a plot, if false, then no plot is created
#                                               errorbands=True, # when you set this to True, then error bands are created around the curves based on the standard deviation
#                                               verbose=True) #if True, an extra dataframe with the std is returned aswell.


#%%

