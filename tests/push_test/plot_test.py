# -*- coding: utf-8 -*-

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 08:29:37 2022

@author: thoverga
"""

import sys, os

from pathlib import Path


lib_folder = Path(__file__).resolve().parents[2]
sys.path.append(str(lib_folder))
print(str(lib_folder))

from src import vlinder_toolkit

#%% IO testdata

testdatafile = os.path.join(str(lib_folder), 'tests', 'test_data',  'vlinderdata_small.csv')
metadatafile = os.path.join(str(lib_folder), 'static_data', 'vlinder_metadata.csv')



settings = vlinder_toolkit.Settings()
settings.update_settings(input_data_file=testdatafile,
                         input_metadata_file=metadatafile)
settings.check_settings()

dataset = vlinder_toolkit.Dataset()
dataset.import_data_from_file(coarsen_timeres=True)

#%% timeseries plots of dataset
dataset.make_plot()
dataset.make_plot(variable='humidity',stationnames=['vlinder02', 'vlinder17'], legend=False)

from datetime import datetime
dataset.make_plot(starttime=datetime(2022,9,4),
                  endtime=datetime(2022,9,6), title='test')

#%% timeseries plot of station

dataset.get_station('vlinder05').make_plot()

#%% Make spatial plot
dataset.make_geo_plot()
dataset.make_geo_plot(variable='wind_direction', timeinstance=datetime(2022,9,5,12,0))