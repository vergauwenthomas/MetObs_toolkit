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

import matplotlib.pyplot as plt
import matplotlib.dates as mdates


ax = dataset.make_plot(legend=True)
ax = dataset.make_plot(colorby='label')


# from matplotlib.dates import AutoDateFormatter, AutoDateLocator

# xtick_locator = AutoDateLocator()
# xtick_formatter = AutoDateFormatter(xtick_locator)

# # ax = plt.axes()
# ax.xaxis.set_major_locator(xtick_locator)
# ax.xaxis.set_major_formatter(xtick_formatter)





# print(ax.xaxis.get_major_formatter())

# # locator = mdates.AutoDateLocator()
# formatter = mdates.AutoDateFormatter(mdates.AutoDateLocator())
# # formatter = mdates.DateFormatter(fmt='%Y/%m/%d %H:%M:%S')
# ax.xaxis.set_major_formatter(formatter)
# # ax.xaxis.set_minor_formatter(formatter)
# print(ax.xaxis.get_major_formatter())

#%%



