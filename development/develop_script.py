#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:25:02 2022
@author: thoverga
"""

#%%

# import metobs_toolkit

import os
import sys
from pathlib import Path


lib_folder = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(lib_folder))

import metobs_toolkit


#%%
import pandas as pd
import datetime


# Make an empty dataset
dataset = metobs_toolkit.Dataset()

# Add the demo data files to the dataset settings
dataset.update_settings(input_data_file = metobs_toolkit.demo_datafile,
                        input_metadata_file = metobs_toolkit.demo_metadatafile,
                        data_template_file = metobs_toolkit.demo_template,
                        metadata_template_file = metobs_toolkit.demo_template, # Contains also the metadata mapping
                        output_folder='/home/thoverga/Documents/VLINDER_github/MetObs_toolkitss'
                        )



# Load the data from the demo data files
dataset.import_data_from_file()

# dataset.coarsen_time_resolution()

# dataset.get_lcz()

#%%
# outputfolder = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkitss'
dataset.save_dataset()



# Make an empty dataset
dataset = metobs_toolkit.Dataset()

# dataset2 = dataset.import_dataset(folder_path=outputfolder)

# print(dataset2)


#%%



# anal = dataset.get_analysis()


# anal.get_aggregated_diurnal_statistics(aggregation=['season', 'name'])



#%%
# from datetime import datetime


# def get_seasons(datetimeseries,
#                 start_day_spring = '01/03' ,
#                 start_day_summer = '01/06',
#                 start_day_autumn = '01/09',
#                 start_day_winter = '01/12'):

#     spring_startday = datetime.strptime(start_day_spring, '%d/%m')
#     summer_startday = datetime.strptime(start_day_summer, '%d/%m')
#     autumn_startday = datetime.strptime(start_day_autumn, '%d/%m')
#     winter_startday = datetime.strptime(start_day_winter, '%d/%m')


#     seasons = pd.Series(index=['spring', 'summer', 'autumn', 'winter'],
#                         data=[spring_startday, summer_startday, autumn_startday, winter_startday],
#                         name='startdt').to_frame()
#     seasons['day_of_year'] = seasons['startdt'].dt.day_of_year - 1

#     bins = [0]
#     bins.extend(seasons['day_of_year'].to_list())
#     bins.append(366)

#     labels = ['winter', 'spring', 'summer', 'autumn', 'winter']



#     return pd.cut(x = datetimeseries.dt.day_of_year,
#                   bins = bins,
#                   labels=labels,
#                   ordered=False,
#                   )




# def aggregate_df(analysis, df, agg=['lcz', 'datetime'], method='mean'):
#     """
#     Aggregate observations to a (list of) categories.

#     The output will be a dataframe that is aggregated to one, or more categories.
#     A commen example is aggregating to LCZ's.


#     Parameters
#     ----------
#     df : pandas.DataFrame
#         The observations to aggregate.
#     agg : list, optional
#         The list of columnnames to aggregate to. If 'lcz' is included, the
#         lcz information is extracted from the Analysis.metadf. The default is ['lcz', 'datetime'].
#     method : str, optional
#         list of functions and/or function names, e.g. [np.sum, 'mean']. The default is 'mean'.

#     Returns
#     -------
#     pandas.DataFrame
#         A dataframe with the agg columns as an index. The values are the aggregated values.

#     Note
#     -------
#     Present columns that ar non-numeric and are not in the agg list are not present in the return,
#     since these values cannot be aggregated.

#     """
#     df = df.reset_index()

#     time_agg_keys = ['minute', 'hour', 'month', 'year', 'day_of_year',
#                      'week_of_year', 'season']

#     # scan trough the metadf for aggregation keys
#     for agg_key in agg:
#         if agg_key not in df.columns:
#             # look in metadf
#             if agg_key in analysis.metadf.columns:
#                 df = pd.merge(df, analysis.metadf[[agg_key]],
#                                   how='left', left_on='name',
#                                   right_index=True)




#     # Check if all agg keys are present or defined:
#     possible_agg_keys = time_agg_keys
#     possible_agg_keys.extend(list(df.columns))
#     unmapped = [agg_key for agg_key in agg if agg_key not in possible_agg_keys]
#     assert len(unmapped) == 0, f'cannot aggregate to unknown labels: {unmapped}.'
#     # define time aggregations
#     #1. minute
#     if 'minute' in agg:
#         df['minute'] = df['datetime'].dt.minute
#     if 'hour' in agg:
#         df['hour'] = df['datetime'].dt.hour
#     if 'month' in agg:
#         df['month'] = df['datetime'].dt.month_name()
#     if 'year' in agg:
#         df['year'] = df['datetime'].dt.year
#     if 'day_of_year' in agg:
#         df['day_of_year'] = df['datetime'].dt.day_of_year
#     if 'week_of_year' in agg:
#         df['week_of_year'] = df['datetime'].dt.week_of_year
#     if 'season' in agg:
#         df['season'] = get_seasons(df['datetime'])


#     # check if not all values are Nan
#     for agg_name in agg:
#         assert df[agg_name].isnull().all() == False, f'Aggregation to {agg_name} not possible because no valid values found for {agg_name}.'





#     # Aggregate the df
#     agg_df = df.groupby(agg).agg(method, numeric_only=True)
#     # sort index
#     agg_df = agg_df.reset_index()
#     agg_df = agg_df.set_index(agg)
#     return agg_df



# test = aggregate_df(analysis = anal,
#                     df = anal.df,
#                     agg=['lcz', 'month'] )




# #%%





# df = anal.df.reset_index()

# df['season'] = get_seasons(df['datetime'])




