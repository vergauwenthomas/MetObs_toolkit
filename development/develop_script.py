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


# testdatafile = os.path.join(str(lib_folder), 'tests', 'test_data',  'testdata_okt_small.csv')
testdatafile = os.path.join(str(lib_folder), 'tests', 'test_data',  'testdata_breaking.csv')


template = os.path.join(str(lib_folder), 'tests', 'test_data',  'template_breaking.csv')

static_data = os.path.join(
    str(lib_folder), 'static_data', 'vlinder_metadata.csv')



# #% Setup dataset

dataset = metobs_toolkit.Dataset()
dataset.update_settings(input_data_file=testdatafile,
                        # input_metadata_file=static_data,
                        data_template_file= template,
                        output_folder='/home/thoverga/Documents'
                        )



dataset.import_data_from_file()
dataset.apply_quality_control(gross_value=True,
                              persistance=False)


dataset.get_qc_stats('humidity')
#%%
applieddf = dataset._applied_qc
print(applieddf)



#%%

obstype='temp'


# cobmine all and get final label
comb_df = dataset.combine_all_to_obsspace()

# subset to relevant columnt
relev_columns = [obstype, obstype+'_final_label']
comb_df = comb_df[relev_columns]


#%%

import pandas as pd

comb_df = comb_df
obstype = obstype
checks_info =dataset.settings.qc['qc_checks_info']
gaps_info =dataset.settings.gap['gaps_info']


# ----------------------------------------------------------------------



outlier_labels = [qc['outlier_flag'] for qc in checks_info.values()]

if not obstype+'_final_label' in comb_df.columns:
    print(f'Final observation label for {obstype} is not computed!')
    # return (None, None, None)

final_counts = comb_df[obstype+'_final_label'].value_counts()

#add missing labels
# QC labels
non_triggered_labels_dict = {}
#fill with zeros for non-triggered checks
for outl_label in outlier_labels:
    if not outl_label in final_counts.index:
        non_triggered_labels_dict[outl_label] = 0

#gaps
if not gaps_info['gap']['outlier_flag'] in final_counts.index:
    non_triggered_labels_dict[gaps_info['gap']['outlier_flag']] = 0

#missing timestamps
if not gaps_info['missing_timestamp']['outlier_flag'] in final_counts.index:
    non_triggered_labels_dict[gaps_info['missing_timestamp']['outlier_flag']] = 0


non_triggered_labels = pd.Series(non_triggered_labels_dict)
final_counts = pd.concat([final_counts, non_triggered_labels])
tot_n_obs= final_counts.sum()

# to percentages
final_counts = (final_counts/tot_n_obs)*100.0



# ------- aggregate outliers ----------



# 1 agg to ok - outlier - gap - missing


try:
    agg_ok = final_counts["ok"].squeeze()
except KeyError:
    agg_ok = 0.


agg_dict = {
    'ok': agg_ok,
    'QC outliers': final_counts.loc[final_counts.index.isin(outlier_labels)].sum(),
    'missing (gaps)': final_counts[gaps_info['gap']['outlier_flag']].squeeze(),
    'missing (individual)': final_counts[gaps_info['missing_timestamp']['outlier_flag']].squeeze()
    }


#2 indevidual outliers
outl_dict =  final_counts.loc[final_counts.index.isin(outlier_labels)].to_dict()


# 3 Effectivenes per check

specific_counts = {}
# Note: some complexity because observations can be removed by privious executed checsk,
# so construct the counts in the order of the applied checks

applied_checks = dataset._applied_qc.loc[dataset._applied_qc['obstype'] == obstype]['checkname'].to_list()

percent_rejected_before = 0.
for checkname in applied_checks:
    try:
        specific_outliers = final_counts.loc[checks_info[checkname]['outlier_flag']]
    except KeyError:
        specific_outliers = 0.

    not_checked = percent_rejected_before
    ok = 100. - specific_outliers - not_checked

    specific_counts[checkname] = {'not checked': not_checked,
                                  'ok': ok,
                                  'outlier': specific_outliers}

    percent_rejected_before += specific_outliers


# add checks that are not performed
not_perf_checknames = [check for check in checks_info.keys() if not check in applied_checks]
for checkname in not_perf_checknames:
    specific_counts[checkname] = {'not checked': 100.,
                                  'ok': 0.,
                                  'outlier': 0.}


# add Gaps
gap_specific_counts = {
    'not checked': 0, #all obs are always checked
    'ok': 100.0 - final_counts[gaps_info['gap']['outlier_flag']],
    'outlier': final_counts[gaps_info['gap']['outlier_flag']]
    }
specific_counts[gaps_info['gap']['label_columnname']] = gap_specific_counts


#misssing timestamps
missing_specific_counts = {
    'not checked': 0, #all obs are always checked
    'ok': 100.0 - final_counts[gaps_info['missing_timestamp']['outlier_flag']],
    'outlier': final_counts[gaps_info['missing_timestamp']['outlier_flag']]
    }
specific_counts[gaps_info['missing_timestamp']['label_columnname']] = missing_specific_counts

# return (agg_dict, outl_dict, specific_counts)


# # compute freq statistics
# final_freq, outl_freq, specific_freq = get_freq_statistics(
#     comb_df = comb_df,
#     obstype=obstype,
#     checks_info=dataset.settings.qc['qc_checks_info'],
#     gaps_info =dataset.settings.gap['gaps_info'],
#     )


test = pd.DataFrame(specific_counts)
print(test)
