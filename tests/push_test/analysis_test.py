#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  8 15:55:13 2023

@author: thoverga
"""


import sys, os

from pathlib import Path
import pandas as pd

import metobs_toolkit

lib_folder = Path(__file__).resolve().parents[2]
# sys.path.append(str(lib_folder))
print(str(lib_folder))


testdatafile = os.path.join(str(lib_folder), "tests", "test_data", 'testdata_okt_small.csv')
static_data = os.path.join(str(lib_folder),'static_data','vlinder_metadata.csv')

lcz_dict = {'vlinder01': 'Low plants (LCZ D)',
 'vlinder02': 'Open midrise',
 'vlinder03': 'Open midrise',
 'vlinder04': 'Sparsely built',
 'vlinder05': 'Water (LCZ G)',
 'vlinder06': 'Scattered Trees (LCZ B)',
 'vlinder07': 'Compact midrise',
 'vlinder08': 'Compact midrise',
 'vlinder09': 'Scattered Trees (LCZ B)',
 'vlinder10': 'Compact midrise',
 'vlinder11': 'Open lowrise',
 'vlinder12': 'Open highrise',
 'vlinder13': 'Compact midrise',
 'vlinder14': 'Low plants (LCZ D)',
 'vlinder15': 'Sparsely built',
 'vlinder16': 'Water (LCZ G)',
 'vlinder17': 'Scattered Trees (LCZ B)',
 'vlinder18': 'Low plants (LCZ D)',
 'vlinder19': 'Compact midrise',
 'vlinder20': 'Compact midrise',
 'vlinder21': 'Sparsely built',
 'vlinder22': 'Low plants (LCZ D)',
 'vlinder23': 'Low plants (LCZ D)',
 'vlinder24': 'Dense Trees (LCZ A)',
 'vlinder25': 'Water (LCZ G)',
 'vlinder26': 'Open midrise',
 'vlinder27': 'Compact midrise',
 'vlinder28': 'Open lowrise'}



dataset = metobs_toolkit.Dataset()
dataset.update_settings(input_data_file=testdatafile,
                        input_metadata_file=static_data,
                        # data_template_file= template,
                        output_folder='/home/thoverga/Documents'
                        )



dataset.import_data_from_file()

dataset.metadf['lcz'] = pd.Series(lcz_dict) #to avoid gee interaction
#%%



an = dataset.get_analysis()

#%%
teststa =  ['vlinder01', 'vlinder02', 'vlinder03']

from datetime import datetime
startdt = datetime(2022,10,6)

#%%
test1 = an.get_diurnal_statistics(colorby='name', stations=teststa, startdt=startdt, verbose=True)

#%%

test2 = an.get_diurnal_statistics_with_reference(refstation='vlinder08',colorby='name', verbose=True)


#%%
test3 = an.get_aggregated_diurnal_statistics(aggregation=['lcz'], verbose=True)


