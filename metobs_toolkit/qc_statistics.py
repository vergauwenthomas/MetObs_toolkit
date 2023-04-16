# -*- coding: utf-8 -*-

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  25 13:44:54 2022

@author: thoverga
"""


import pandas as pd
import logging



logger = logging.getLogger(__name__)

def get_freq_statistics(comb_df, obstype, checks_info, gaps_info):

    outlier_labels = [qc['outlier_flag'] for qc in checks_info.values()]

    if not obstype+'_final_label' in comb_df.columns:
        print(f'Final observation label for {obstype} is not computed!')
        return (None, None, None)

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

    agg_dict = {
        'ok': final_counts["ok"].squeeze(),
        'QC outliers': final_counts.loc[final_counts.index.isin(outlier_labels)].sum(),
        'missing (gaps)': final_counts[gaps_info['gap']['outlier_flag']].squeeze(),
        'missing (individual)': final_counts[gaps_info['missing_timestamp']['outlier_flag']].squeeze()
        }


    #2 indevidual outliers
    outl_dict =  final_counts.loc[final_counts.index.isin(outlier_labels)].to_dict()


    # 3 Effectivenes per check

    qc_label_columns = [col for col in comb_df.columns if not col in [obstype, obstype+'_final_label' ] ]

    mapper = {qc_type['outlier_flag']: 'outlier' for qc_type in checks_info.values()}
    specific_counts = comb_df.replace(mapper).reset_index()[qc_label_columns] \
                            .transpose().apply(pd.Series.value_counts, axis=1).fillna(0) #\
                            # .to_dict(orient='index')

    # convert to percentages and dict
    specific_counts = ((specific_counts/tot_n_obs) * 100).to_dict(orient='index')

    # specific_counts[Settings.gaps_info['gap']['label_columnname']] = {}
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

    return (agg_dict, outl_dict, specific_counts)




