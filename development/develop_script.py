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
from pathlib import Path


lib_folder = Path(__file__).resolve().parents[1]
sys.path.append(str(lib_folder))




#%% % Import


testdatafile = os.path.join(
    str(lib_folder), 'tests', 'test_data',  'testdata_okt_small.csv')

static_data = os.path.join(
    str(lib_folder), 'static_data', 'vlinder_metadata.csv')



# #% Setup dataset



dataset = vlinder_toolkit.Dataset()
dataset.update_settings(input_data_file=testdatafile,
                        input_metadata_file=static_data,
                        output_folder='/home/thoverga/Documents/VLINDER_github/vlinder_toolkit'
                        )

# dataset.apply_quality_control()


dataset.import_data_from_file(coarsen_timeres=True)
# dataset.apply_quality_control()


<<<<<<< HEAD
=======

dataset.import_data_from_file(coarsen_timeres=True)


# dataset.set_timezone()

# # dataset.import_data_from_file(coarsen_timeres=True)

# # dataset.apply_quality_control()
# # dataset.get_qc_stats()
dataset.make_geo_plot()
dataset.make_plot()
>>>>>>> master

dataset.make_plot(stationnames=['vlinder01', 'vlinder02'])


<<<<<<< HEAD
# era_modeldata = dataset.get_modeldata(stations='vlinder01')

# dataset.fill_gaps_era5(era_modeldata)








#%%
=======
#%%
from datetime import datetime
import pandas as pd

startdt = datetime(2023, 3,24)
enddt = datetime(2023, 3, 28)


df = pd.DataFrame()




df['datetime'] = pd.date_range(startdt, enddt, freq='60T')


df['date_lt'] = pd.to_datetime(df['datetime']).dt.tz_localize('Europe/Brussels',
                                                                  ambiguous='infer',
                                                                  nonexistent='shift_forward')

df['date_utc'] =pd.to_datetime(df['date_lt'].dt.tz_convert('UTC'))



import pytz
test = pytz.common_timezones


#%%
test = dataset.df.copy()
st# test.index. = test.index.get_level_values('datetime').tz_localize('Europe/Brussels',
#                                                                   ambiguous='infer',
#                                                                   nonexistent='shift_forward')

test.index = test.index.set_levels(test.index.levels[1].tz_localize('Europe/Brussels',
                                                                  ambiguous='infer',
                                                                  nonexistent='shift_forward'), level=1)

>>>>>>> master
