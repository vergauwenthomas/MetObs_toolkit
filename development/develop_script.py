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
import logging
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

# handle logs to spyder console
toolkit_logger = metobs_toolkit.loggers

toolkit_logger.addHandler(logging.StreamHandler())


#%%
import os


newmod = metobs_toolkit.Modeldata('ERA5_hourly')
newmod = newmod.import_modeldata(folder_path = '/home/thoverga/Documents/VLINDER_github/MetObs_GUI/metobs_gui/cache/modeldata',
                               filename='my_modeldata.pkl')

newmod.make_plot()






#%%
# use_dataset = 'debug_wide'
# use_dataset = 'single_netatmo_sara_station'
use_dataset = 'vlindergent2022'
use_dataset = 'demo'

dataset = metobs_toolkit.Dataset()


dataset.update_settings(output_folder=None,
                        input_data_file=testdata[use_dataset]['datafile'],
                        input_metadata_file=testdata[use_dataset]['metadatafile'],
                        template_file=testdata[use_dataset]['template'],
                        )




dataset.import_data_from_file(**testdata[use_dataset]['kwargs'])

dataset.coarsen_time_resolution(freq = testdata[use_dataset]['coarsen'])
dataset.apply_quality_control()
# dataset.get_lcz()
dataset.update_gaps_and_missing_from_outliers()
#%%

dataset.fill_missing_obs_linear()

dataset.fill_gaps_linear()



#%%

start_gaps = dataset.get_gaps_df()['start_gap'].min()
end_gaps =  dataset.get_gaps_df()['end_gap'].max()




#%%
dataset.update_qc_settings(obstype='radiation_temp', rep_max_valid_repetitions=3)
dataset.apply_quality_control(obstype='radiation_temp')

#%%
dataset.make_plot(obstype='radiation_temp', colorby='label')




#%%
dataset.get_qc_stats(obstype='radiation_temp')


df = dataset.df[['radiation_temp']]
outliers = dataset.outliersdf.xs('radiation_temp', level='obstype')




comb_df = dataset.combine_all_to_obsspace()
comb_df = comb_df.xs('radiation_temp', level='obstype')
#%%





#%%
dataset.update_gaps_and_missing_from_outliers()
dataset.fill_gaps_linear()


#%%
# use_dataset = 'debug_wide'
# dataset = metobs_toolkit.Dataset()


# dataset.update_settings(output_folder=None,
#                         input_data_file=testdata[use_dataset]['datafile'],
#                         input_metadata_file=testdata[use_dataset]['metadatafile'],
#                         template_file=testdata[use_dataset]['template'],
#                         )



# dataset.import_data_from_file()


#%%
# # use_dataset = 'debug_wide'
# use_dataset = 'single_netatmo_sara_station'
# use_dataset = 'vlindergent2022'
# use_dataset = 'demo'
# dataset = metobs_toolkit.Dataset()


# dataset.update_settings(output_folder=None,
#                         input_data_file=testdata[use_dataset]['datafile'],
#                         input_metadata_file=testdata[use_dataset]['metadatafile'],
#                         data_template_file=testdata[use_dataset]['template'],
#                         metadata_template_file=testdata[use_dataset]['template'],
#                         )



# dataset.import_data_from_file(**testdata[use_dataset]['kwargs'])

# dataset.coarsen_time_resolution(freq = testdata[use_dataset]['coarsen'])
# dataset.apply_quality_control()
# dataset.get_lcz()
# dataset.update_gaps_and_missing_from_outliers()
#%%
# dataset.apply_quality_control()
# dataset.update_gaps_and_missing_from_outliers()
# dataset.fill_gaps_linear()

#%%
