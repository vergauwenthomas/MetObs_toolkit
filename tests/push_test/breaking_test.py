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

#####################################################################
# Set settings for QC
minimal_gapsize = 40    #gaps defined as n times the highest frequency on IO.
dupl_dropping = False #method used to drop duplicated timestamps

persistance_time_window_to_check = '1h'   # Use this format as example: "1h20min50s"
min_num_obs = 5   #Minimum number of records in window to perform persistance check

max_valid_repetitions = 5 # Maximal number of repetitions that is allowed

min_value = -15.0  # Minimal allowed value
max_value = 50.0 # Maximal allowed value

max_increase_per_second = 8.0/3600.0   # Maximal allowed increase per second (for window variation check)
max_decrease_per_second = 10.0/3600.0   # Maximal allowed decrease per second (for window variation check)
time_window_to_check = '1h'   # Use this format as example: "1h20min50s"
min_window_members = 3 # Minimal number of records in window to perform check

max_increase_per_second_step = 8.0/3600.0   # Maximal allowed increase per second (for step check)
max_decrease_per_second_step = -10.0/3600.0   # Maximal allowed increase per second (for step check)


settings.qc_check_settings['gaps_finder']['gapsize_n'] = minimal_gapsize

settings.qc_check_settings['duplicated_timestamp']['keep'] = dupl_dropping

settings.qc_check_settings['persistance']['temp']['time_window_to_check'] = persistance_time_window_to_check
settings.qc_check_settings['persistance']['temp']['min_num_obs'] = min_num_obs

settings.qc_check_settings['repetitions']['temp']['max_valid_repetitions'] = max_valid_repetitions

settings.qc_check_settings['gross_value']['temp']['min_value'] = min_value
settings.qc_check_settings['gross_value']['temp']['max_value'] = max_value

settings.qc_check_settings['window_variation']['temp']['max_increase_per_second'] = max_increase_per_second
settings.qc_check_settings['window_variation']['temp']['max_decrease_per_second'] = max_decrease_per_second
settings.qc_check_settings['window_variation']['temp']['time_window_to_check'] = time_window_to_check
settings.qc_check_settings['window_variation']['temp']['min_window_members'] = min_window_members

settings.qc_check_settings['step']['temp']['max_increase_per_second'] = max_increase_per_second_step
settings.qc_check_settings['step']['temp']['max_decrease_per_second'] = max_decrease_per_second_step
#####################################################################


#% add template
template_file = os.path.join(str(lib_folder), 'tests', 'test_data',  'template_breaking.xls')
settings.add_excel_template(template_file)

dataset_coarsened = vlinder_toolkit.Dataset()
dataset_coarsened.import_data_from_file(coarsen_timeres=True)
dataset_coarsened.apply_quality_control()


#%%
dataset = vlinder_toolkit.Dataset()
dataset.import_data_from_file(coarsen_timeres=False)
dataset.apply_quality_control()

outliersdf = dataset.get_final_qc_labels()
df = dataset.input_df

indices_dupl_df = df[df['flags'] =='duplicated timestamp outlier'].index.sortlevel()[0]
df = df[df['flags'] != 'duplicated timestamp outlier']

indices_dupl_outliers = outliersdf[outliersdf['temp_final_label'] =='duplicated timestamp outlier'].index.sortlevel()[0]
outliersdf = outliersdf[outliersdf['temp_final_label'] != 'duplicated timestamp outlier']

if not indices_dupl_df.equals(indices_dupl_outliers):
    if len(indices_dupl_outliers.difference(indices_dupl_df)) > 0:
        print('Timestamps with wrong duplicate label are: ', indices_dupl_outliers.difference(indices_dupl_df))
    else:
        print('Timestamps with missing duplicate label are: ', indices_dupl_df.difference(indices_dupl_outliers))
    sys.exit('There is a problem with the duplicates')

df = df.merge(outliersdf['temp_final_label'], how='outer', left_index=True, right_index=True)
df['temp_final_label'] = df['temp_final_label'].fillna(value='ok')

indices_missing_timestamp = df[df['temp_final_label'] == 'missing timestamp'].index
df.loc[indices_missing_timestamp,'flags'] = 'missing timestamp'

indices_gap_timestamp = df[df['temp_final_label'] == 'missing timestamp (gap)'].index
df.loc[indices_gap_timestamp,'flags'] = 'missing timestamp (gap)'

dataset.get_qc_stats()
dataset_coarsened.get_qc_stats()

if not df['flags'].equals(df['temp_final_label']):
    print('Timestamps with wrong label are: ', list(df.index[df['flags'] != df['temp_final_label']]))
    sys.exit('There is a problem with the quality control')
else:
    print('The quality control is performing as expected')
    


