#!/usr/bin/env python3
# -*- coding: utf-8 -*-



# PROGRAM FOR TESTING THE BREAKING DATAFILE
"""
Created on Tue Nov 29 12:19:03 2022

@author: mivieijra
"""

import sys, os

from pathlib import Path

import vlinder_toolkit
lib_folder = Path(__file__).resolve().parents[2]
# sys.path.append(str(lib_folder))
print(str(lib_folder))




#x = all(keys in ['a', 'b', 'c'] for keys in ['c', 'b', 'a'])

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


settings.gaps_settings['gaps_finder']['gapsize_n'] = minimal_gapsize

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
template_file = os.path.join(str(lib_folder), 'tests', 'test_data',  'template_breaking.csv')
settings.add_csv_template(template_file)

dataset_coarsened = vlinder_toolkit.Dataset()
dataset_coarsened.import_data_from_file(coarsen_timeres=True)
dataset_coarsened.apply_quality_control()

_ = dataset_coarsened.get_qc_stats()
#%%
dataset = vlinder_toolkit.Dataset()
dataset.import_data_from_file(coarsen_timeres=False)
dataset.apply_quality_control()

dataset.make_plot(stationnames=['Fictional'],colorby='label', show_outliers=True)


#%% Compare manual and toolkit labeling


man_df = dataset.input_df #manual label

tlk_df = dataset.combine_all_to_obsspace()



#%%
all_manual_labels = list(man_df['flags'].unique())
manual_to_tlkit_label_map = {
     'ok': 'ok',
     'in step outlier group': settings.qc_checks_info['step']['outlier_flag'],
     'repetitions outlier': settings.qc_checks_info['repetitions']['outlier_flag'],
     'duplicated timestamp outlier': settings.qc_checks_info['duplicated_timestamp']['outlier_flag'],
     'gross value outlier': settings.qc_checks_info['gross_value']['outlier_flag'],
     'in window variation outlier group': settings.qc_checks_info['window_variation']['outlier_flag'],
     'persistance outlier': settings.qc_checks_info['persistance']['outlier_flag']
    }

#check if the mapper is still up to date
assert all([True for label in all_manual_labels if label in manual_to_tlkit_label_map.keys()]), 'Update the manual to toolkit mapper'



# =============================================================================
# iterate over all labels and validate if the indices are equal between manual and toolkit
# =============================================================================

for man_label, tlk_label in manual_to_tlkit_label_map.items():
    print(f' Testing equality of the {tlk_label} with the manual labeling ({man_label}).')
    
    man_idx = man_df[man_df['flags'] == man_label].index.sort_values()
    tlk_idx = tlk_df[tlk_df['temp_final_label'] == tlk_label].index.sort_values()

    
    if not tlk_idx.equals(man_idx):
        print(f'ERROR: wrong labels for {tlk_label}')
        
        print(f'differences tlkit --> manual: { tlk_idx.difference(man_idx)}')
        print(f'differences manual --> tlkit: {man_idx.difference(tlk_idx)}')
        sys.exit(1)
    
    else:
        print('OK!')


# =============================================================================
# test missing Gaps
# =============================================================================



from datetime import datetime
import pandas as pd

manual_missing_gaps = [{'name': 'Fictional', 'start_gap': datetime(2020,9,14,1,30), 'end_gap': datetime(2020,9,14,23,55)}] #UPDATE MANUALLY !!!!!!!!!!

print('Testing the gaps')

man_gapsdf = pd.DataFrame().from_records(manual_missing_gaps)
man_gapsdf = man_gapsdf.set_index('name')

tlk_gapsdf = dataset.gaps.to_df()
tlk_gapsdf = tlk_gapsdf[list(man_gapsdf.columns)]



if not tlk_gapsdf.equals(man_gapsdf):
    print(f'ERROR: wrong gaps detection')
    
    print(f'differences tlkit --> manual: { tlk_gapsdf.difference(man_gapsdf)}')
    print(f'differences manual --> tlkit: {man_gapsdf.difference(tlk_gapsdf)}')
    sys.exit(1)

else:
    print('OK!')



# =============================================================================
# test missing missing timestamps
# =============================================================================


#This has to be done properly, there are too much missing timestamps for manual labelling. 
# as a dirty fix, now only count the number of missing timestamps and see it that is oke

number_missing_timestamps = {'1': 1,
                             'Fictional' : 307}

print('Testing the missing obs')
tlk_missing_series = dataset.missing_obs.series

for station in tlk_missing_series.index.unique():
    sta_n_missing = tlk_missing_series[[station]].shape[0]
    
    if sta_n_missing != number_missing_timestamps[station]:
        print(f'ERROR: wrong number of missing obs for station: {station}')
        
        print(f'number of missing obs for {station} by manual work: {number_missing_timestamps[station]}')
        print(f'number of missing obs for {station} by toolkit: {tlk_missing_series}')
        sys.exit(1)
    
    
print('OK!')
    









# # indices_dupl_df = df[df['flags'] =='duplicated timestamp outlier'].index.sortlevel()[0]
# # df = df[df['flags'] != 'duplicated timestamp outlier']

# # indices_dupl_outliers = outliersdf[outliersdf['temp_final_label'] =='duplicated timestamp outlier'].index.sortlevel()[0]
# # outliersdf = outliersdf[outliersdf['temp_final_label'] != 'duplicated timestamp outlier']


# df = df[df['flags'] != 'duplicated timestamp outlier']
# outliersdf = outliersdf[outliersdf['temp_final_label'] != 'duplicated timestamp outlier']






# if not indices_dupl_df.equals(indices_dupl_outliers):
#     if len(indices_dupl_outliers.difference(indices_dupl_df)) > 0:
#         print('Timestamps with wrong duplicate label are: ', indices_dupl_outliers.difference(indices_dupl_df))
#     else:
#         print('Timestamps with missing duplicate label are: ', indices_dupl_df.difference(indices_dupl_outliers))
#     sys.exit('There is a problem with the duplicates')

# df = df.merge(outliersdf['temp_final_label'], how='outer', left_index=True, right_index=True)
# df['temp_final_label'] = df['temp_final_label'].fillna(value='ok')

# indices_missing_timestamp = df[df['temp_final_label'] == 'missing timestamp'].index
# df.loc[indices_missing_timestamp,'flags'] = 'missing timestamp'

# indices_gap_timestamp = df[df['temp_final_label'] == 'missing timestamp (gap)'].index
# df.loc[indices_gap_timestamp,'flags'] = 'missing timestamp (gap)'


# dataset.get_qc_stats()



# if not df['flags'].equals(df['temp_final_label']):
#     print('Timestamps with wrong label are: ', list(df.index[df['flags'] != df['temp_final_label']]))
#     sys.exit('There is a problem with the quality control')
# else:
#     print('The quality control is performing as expected')
    


