#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 08:59:40 2022

@author: thoverga
"""

import sys
print('SYS path: ', list(sys.path))
print('Python version: ', sys.version)
#TO show installed package location and info:
    
# pip show vlinder_toolkit
import vlinder_toolkit


print('Succesfull imported!')


#%% Test
import os
testdatafile = os.path.join('/home/thoverga/Documents/VLINDER_github/vlinder_toolkit/tests/test_data',  'vlinderdata.csv')


settings = vlinder_toolkit.Settings()
settings.show()


settings.update_settings(input_file=testdatafile)
settings.check_settings()


et = vlinder_toolkit.Dataset()

dataset.import_data_from_file()


station = dataset.get_station('vlinder02')
stationdf = station.df()
print(stationdf.head())




try:
    ax =station.make_plot()
except:
    print('coulnd not make plot ...')
    
    
print("FINISH")