#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:25:02 2022
@author: thoverga
"""

#%%

import metobs_toolkit

import os
import sys
from pathlib import Path


lib_folder = Path(__file__).resolve().parents[1]
sys.path.append(str(lib_folder))




#%% % Import


testdatafile = os.path.join(str(lib_folder), 'tests', 'test_data',  'testdata_okt_small.csv')

static_data = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/static_data/vlinder_metadata.csv'
# template = os.path.join(str(lib_folder), 'tests', 'test_data',  'wide_test_template.csv')
# template = os.path.join(str(lib_folder), 'tests', 'test_data',  'wide_test_template.csv')




# #% Setup dataset

dataset = metobs_toolkit.Dataset()
dataset.update_settings(input_data_file=testdatafile,
                        input_metadata_file=static_data,
                        # data_template_file= template,
                        output_folder='/home/thoverga/Documents'
                        )



dataset.import_data_from_file()
# dataset.get_lcz()





#%%

an = dataset.get_analysis()


teststa =  ['vlinder01', 'vlinder02']

test2 = an.get_diurnal_statistics_with_reference(refstation='vlinder08',colorby='name', verbose=True, errorbands=True)



from datetime import datetime
startdt = datetime(2022,10,6)

#%%
# test1 = an.get_diurnal_statistics(colorby='name', stations=teststa, startdt=startdt, verbose=True)
# an.get_diurnal_statistics(colorby='name', stations=teststa, errorbands=True)
#%%

# test2 = an.get_diurnal_statistics_with_reference(refstation='vlinder08',colorby='name', verbose=True)


#%%
# test3 = an.get_aggregated_diurnal_statistics(aggregation=['lcz'], verbose=True)

obstype = 'temp'
template = dataset.data_template






