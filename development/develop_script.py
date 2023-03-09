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

testdatafile = os.path.join(str(lib_folder), 'tests', 'test_data',  'testdata_okt_small.csv')

static_data =  os.path.join(str(lib_folder), 'static_data', 'vlinder_metadata.csv')




# #% Setup dataset
settings = vlinder_toolkit.Settings()
settings.update_settings(input_data_file=testdatafile,
                          input_metadata_file=static_data,
                          output_folder='/home/thoverga/Documents/VLINDER_github/vlinder_toolkit'
                          )


dataset = vlinder_toolkit.Dataset()
dataset.import_data_from_file(coarsen_timeres=True)


dataset.apply_quality_control()

dataset.write_to_csv()






# test = dataset.get_qc_stats()



#%%



#%% 
sta = dataset.get_station('vlinder01')


def add_final_label_to_outliersdf(outliersdf, data_res_series):
    """
    V3
        This function creates a final label based on de individual qc labels. The final label will be that of the individual qc-label
        which rejected the obseration.
       
        This functions converts labels to numeric values, algebra to get final label, and inversly 
        convert to labels. This is faster than looping over the rows.

        Parameters
        ----------
        outliersdf : pandas.DataFrame 
            The dataset outliers dataframe containing the observations and QC labels. 
       
        data_res_series : Pandas.Series
            The series that contain the dataset resolution (values) per station (index). This 
            is stored in the dataset.metadf as column 'dataset_resolution'. These are used to explode the gaps.

        Returns
        -------
        outliersdf : pd.DataFrame
            The outliersdf with extra columns indicated by example 'temp_final_label' and 'humid_final_label'.

        """ 
       

      
    
    # order columns
    labels_columns = [column for column in outliersdf.columns if not column in settings.observation_types]
    #drop final columns if they are in the outliersdf
    labels_columns = [column for column in labels_columns if not column.endswith('_final_label')]
    
    checked_obstypes = [obstype for obstype in settings.observation_types if any([qc_column.startswith(obstype+'_') for qc_column in labels_columns])]
    columns_on_record_lvl = [info['label_columnname'] for checkname, info in settings.qc_checks_info.items() if info['apply_on'] == 'record']
   
    
    # Construct numeric mapper
    labels_to_numeric_mapper = {info['outlier_flag']:info['numeric_flag'] for info in settings.qc_checks_info.values()}
    # add 'ok' and 'not checked' labels
    labels_to_numeric_mapper['ok'] = 0
    labels_to_numeric_mapper['not checked'] = np.nan
    #invert numeric mapper
    inv_label_to_num = {v: k for k, v in labels_to_numeric_mapper.items()}
    
    
    #generete final label per obstype
    for obstype in checked_obstypes:
        # logger.debug(f'Generating final QC labels for {obstype}.')
        #Get qc column namse specific for this obstype
        specific_columns = [col for col in labels_columns if col.startswith(obstype+'_')]
        #add qc labels that are applicable on all obstypes
        specific_columns.extend(columns_on_record_lvl)
        
        #Drop columns that are not present
        specific_columns = [colmn for colmn in specific_columns if colmn in outliersdf.columns]
        
        
        
        #get labels dataframe
        qc_df = outliersdf[specific_columns]
        num_qc_df = pd.DataFrame()
        num_qc_df = qc_df.applymap(labels_to_numeric_mapper.get )
       
    
        outliersdf[obstype+'_final_label'] = num_qc_df.sum(axis=1, skipna=True).map(inv_label_to_num)
       
   
    return outliersdf

def remove_outliers_from_obs(obsdf, outliersdf):
    return obsdf.loc[~obsdf.index.isin(outliersdf.index)]


# def combine_all_to_obsspace(self):
#     """
#     Combine observations, outliers, gaps and missing timesteps to one dataframe in the resolution of the dataset.
#     Final quality labels are calculated for all checked obstypes.
    
#     If an observation value exist for an outlier, it will be used in the corresponding obstype column.

#     Returns
#     -------
#     comb_df : pandas.DataFrame()
#         Multi index dataframe with observations and labels.

#     """

# from .df_helpers import remove_outliers_from_obs
outliersdf = sta.outliersdf

# get final label
sta_resolution = pd.Series(index=[sta.name],
                           data=[sta.meta_series['dataset_resolution']],
                           name='dataset_resolution')

outliersdf = add_final_label_to_outliersdf(outliersdf, sta_resolution)

# add gaps observations and fill with default values
gapsidx = sta.gaps.get_gaps_indx_in_obs_space(sta.df, sta.outliersdf, sta_resolution)
gapsdf = gapsidx.to_frame()


# add missing observations if they occure in observation space
missingidx = sta.missing_obs.get_missing_indx_in_obs_space(sta.df, sta_resolution)
missingdf = missingidx.to_frame()


#get observations
df = sta.df
# remove outliers (represented by Nan's) from observations
df = remove_outliers_from_obs(df, outliersdf)
   

   
#initiate default values
for col in outliersdf.columns:
    if  col in settings.observation_types:
        default_value_gap = np.nan #nan for observations
        default_value_missing = np.nan
        default_value_obs = df[col]
    elif col.endswith('_final_label'):
        default_value_gap =  settings.gaps_info['gap']['outlier_flag'] #'gap' for final label
        default_value_missing =  settings.gaps_info['missing_timestamp']['outlier_flag'] #'is_missing_timestamp' for final label
        default_value_obs = 'ok'
    else: 
        default_value_gap = 'not checked'
        default_value_missing = 'not checked'
        default_value_obs = 'ok'
    
    gapsdf[col] = default_value_gap
    missingdf[col] = default_value_missing
    df[col] = default_value_obs


# sort columns
gapsdf = gapsdf[list(outliersdf.columns)]
missingdf = missingdf[list(outliersdf.columns)]
df = df[list(outliersdf.columns)]

#Merge all together
comb_df = pd.concat([df, outliersdf, gapsdf, missingdf]).sort_index()
# return comb_df