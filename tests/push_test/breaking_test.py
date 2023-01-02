#!/usr/bin/env python3
# -*- coding: utf-8 -*-



# PROGRAM FOR TESTING THE BREAKING DATAFILE
"""
Created on Tue Nov 29 12:19:03 2022

@author: mivieijra
"""

import sys, os

from pathlib import Path


lib_folder = Path(__file__).resolve().parents[2]
sys.path.append(str(lib_folder))
print(str(lib_folder))

from src import vlinder_toolkit

testdata = os.path.join(str(lib_folder), 'tests', 'test_data',  'testdata_breaking.csv')

settings = vlinder_toolkit.Settings()
settings.update_settings(input_data_file=testdata)

#% add template
template_file = os.path.join(str(lib_folder), 'tests', 'test_data',  'template_breaking.xls')
settings.add_excel_template(template_file)

dataset_coarsened = vlinder_toolkit.Dataset()
dataset_coarsened.import_data_from_file(coarsen_timeres=True)

try:
    dataset_coarsened.apply_quality_control()
except:
    print("The quality control doesn't work when coarsening the dataset")

#%%
dataset = vlinder_toolkit.Dataset()
dataset.import_data_from_file(coarsen_timeres=False)

dataset.apply_quality_control()
outliersdf = dataset.get_final_qc_labels()
df = dataset.input_df
df = df.merge(outliersdf['temp_final_label'], how='outer', left_index=True, right_index=True)
df['temp_final_label'] = df['temp_final_label'].fillna(value='ok')

indices_missing_timestamp = df[df['temp_final_label'] == 'missing timestamp'].index
df.loc[indices_missing_timestamp,'flags'] = 'missing timestamp'

indices_gap_timestamp = df[df['temp_final_label'] == 'missing timestamp (gap)'].index
df.loc[indices_gap_timestamp,'flags'] = 'missing timestamp (gap)'

#dataset.get_qc_stats()

if not df['flags'].equals(df['temp_final_label']):
    print('Timestamps with wrong label are: ', list(df.index[df['flags'] != df['temp_final_label']]))
    #sys.exit('There is a problem with the quality control')
#else:
    #sys.exit('The quality control is performing as expected')
    


