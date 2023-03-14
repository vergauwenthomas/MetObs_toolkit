#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:25:02 2022
@author: thoverga
"""

# import vlinder_toolkit
import os, sys
import pandas as pd
import numpy as np
from pathlib import Path


lib_folder = Path(__file__).resolve().parents[1]
sys.path.append(str(lib_folder))
print(str(lib_folder))

import vlinder_toolkit

# % Import

testdatafile = os.path.join(str(lib_folder), 'tests', 'test_data',  'testdata_breaking.csv')

##tatic_data =  os.path.join(str(lib_folder), 'static_data', 'vlinder_metadata.csv')




# #% Setup dataset
settings = vlinder_toolkit.Settings()
settings.update_settings(input_data_file=testdatafile,
                          input_metadata_file=None,
                          # geotiff_lcz_file=lcz_map
                          output_folder='/home/thoverga/Documents/VLINDER_github/vlinder_toolkit'
                          )

template_file = os.path.join(str(lib_folder), 'tests', 'test_data',  'template_breaking.xls')
settings.add_excel_template(template_file)

dataset = vlinder_toolkit.Dataset()

dataset.import_data_from_file(coarsen_timeres=True)


dataset.apply_quality_control()

dataset.write_to_csv()






# test = dataset.get_qc_stats()



#%%



