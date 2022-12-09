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

from src import vlinder_toolkit

# % Import

testdatafile = os.path.join(str(lib_folder), 'tests', 'test_data',  'testdata_okt_small.csv')

static_data =  os.path.join(str(lib_folder), 'static_data', 'vlinder_metadata.csv')

# lcz_map = os.path.join(str(lib_folder), 'physiograpy', 'lcz_filter_v1.tif')


#% Setup dataset
settings = vlinder_toolkit.Settings()
settings.update_settings(input_data_file=testdatafile,
                          input_metadata_file=static_data,
                          # geotiff_lcz_file=lcz_map
                          output_folder='/home/thoverga/Documents/VLINDER_github/vlinder_toolkit'
                          )


dataset = vlinder_toolkit.Dataset()
dataset.import_data_from_file(coarsen_timeres=False)



dataset.apply_quality_control()


dataset.write_to_csv(filename='remove_me', add_final_labels=True)




#%%


filename="remove_me"
include_outliers=True
add_final_labels=True


observation_types = ['temp', 'radiation_temp', 'humidity', 'precip',
                     'precip_sum', 'wind_speed', 'wind_gust', 'wind_direction',
                     'pressure', 'pressure_at_sea_level']

location_info = ['network', 'lat', 'lon', 'lcz', 'call_name', 'location' ]

check_settings = settings.qc_check_settings
checks_info = settings.qc_checks_info

def init_outlier_multiindexdf():
    my_index = pd.MultiIndex(levels=[['name'],['datetime']],
                             codes=[[],[]],
                             names=[u'name', u'datetime'])

   
    df = pd.DataFrame(index=my_index)
    return df



def gaps_to_outlier_format(gapsdf, dataset_res_series):

    checkname = 'gaps_finder'
    exploded_gaps_df = init_outlier_multiindexdf()
    for station, gapinfo in gapsdf.iterrows():
        gap_timestamps=pd.date_range(start=gapinfo['start_gap'],
                                     end=gapinfo['end_gap'],
                                     freq=dataset_res_series.loc[station])
        multi_idx = pd.MultiIndex.from_tuples(list(zip([station] * len(gap_timestamps),
                                                       gap_timestamps)),
                                              names=['name', 'datetime'])
        exploded_gaps_df = pd.concat([exploded_gaps_df,
                                      pd.DataFrame(data=checks_info[checkname]['outlier_flag'],
                                                   index=multi_idx,
                                                   columns=[checks_info[checkname]['label_columnname']])])
    return exploded_gaps_df

def add_final_label_to_outliersdf(outliersdf, gapsdf, data_res_series):
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
        gapsdf : pandas.Dataframe
            The dataset gaps dataframe. The gaps will be exploded and added to the outliersdf.
        data_res_series : Pandas.Series
            The series that contain the dataset resolution (values) per station (index). This 
            is stored in the dataset.metadf as column 'dataset_resolution'. These are used to explode the gaps.

        Returns
        -------
        outliersdf : pd.DataFrame
            The outliersdf with extra columns indicated by example 'temp_final_label' and 'humid_final_label'.

        """ 
       

    #combine gaps with outliers
    gaps_exploded = gaps_to_outlier_format(gapsdf, data_res_series)
    outliersdf = pd.concat([outliersdf, gaps_exploded])
    #fill Nan's by 'ok'
    outliersdf[settings.qc_checks_info['gaps_finder']['label_columnname']].fillna(value='ok',
                                                                      inplace=True)
    
    # order columns
    labels_columns = [column for column in outliersdf.columns if not column in observation_types]
    checked_obstypes = [obstype for obstype in observation_types if any([qc_column.startswith(obstype+'_') for qc_column in labels_columns])]
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











# logger.info('Writing the dataset to a csv file')
assert not isinstance(settings.output_folder, type(None)), 'Specify settings.output_folder in order to export a csv.'

#add final quality control labels per observation type
if add_final_labels:
    outliersdf = add_final_label_to_outliersdf(outliersdf=dataset.outliersdf,
                                      gapsdf= dataset.gapsdf,
                                      data_res_series = dataset.metadf['dataset_resolution'])
else:
    outliersdf = dataset.outliersdf
 

#Get observations and metadata columns in the right order
# logger.debug('Merging data and metadata')


#make column ordering
df_columns = observation_types.copy() #observations
df_columns.extend(location_info) #metadata


qc_columns = [col for col in outliersdf if col.endswith('_label')] #add qc labels
df_columns.extend(qc_columns)
df_columns.insert(0, 'datetime') # timestamp as first column
df_columns.insert(1, 'name') #station name as second column




#unstack observations and merge with metadf
df = dataset.df 
df[qc_columns] = 'ok'
df = pd.concat([df, outliersdf])
df = df.reset_index()

metadf = dataset.metadf.reset_index()
df = df.merge(metadf, how='left', on='name')

#sort and subset columns
df = df[df_columns]

        
#find observation type that are not present
ignore_obstypes = [col for col in observation_types if df[col].isnull().all()]

df = df.drop(columns=ignore_obstypes)

# logger.debug(f'Skip quality labels for obstypes: {ignore_obstypes}.')

df = df.sort_values(['name', 'datetime'])

#make filename
if isinstance(filename, type(None)):
    startstr = dataset.df.index.min().strftime('%Y%m%d') 
    endstr = dataset.df.index.max().strftime('%Y%m%d') 
    filename= 'dataset_' + startstr + '_' + endstr
else:
    if filename.endswith('.csv'):
        filename = filename[:-4] #to avoid two times .csv.csv
    
filepath = os.path.join(settings.output_folder, filename + '.csv')

#write to csv in output folder
# logger.info(f'write dataset to file: {filepath}')
# df.to_csv(path_or_buf=filepath,
#                sep=';',
#                na_rep='NaN',
#                index=True)        





#%%




