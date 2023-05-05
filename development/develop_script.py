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


testdatafile = os.path.join(str(lib_folder), 'tests', 'test_data',  'testdata_okt_small.csv')


# template = os.path.join(str(lib_folder), 'tests', 'test_data',  'wide_test_template.csv')
# template = os.path.join(str(lib_folder), 'tests', 'test_data',  'wide_test_template.csv')




# #% Setup dataset

dataset = metobs_toolkit.Dataset()
dataset.update_settings(input_data_file=testdatafile,
                        # input_metadata_file=static_data,
                        # data_template_file= template,
                        output_folder='/home/thoverga/Documents'
                        )



dataset.import_data_from_file()


#%%
from metobs_toolkit.analysis import Analysis
an = Analysis(obsdf = dataset.df,
              metadf = dataset.metadf,
              settings = dataset.settings)


test = an.get_diurnal_statistics(refstation='vlinder08')


print(test)


#%%
# obstype='temp'
# stations = ['vlinder01', 'vlinder02', 'vlinder03']
# verbose=True




# obsdf = an.df[obstype]

# # Subset to selected stations if needed
# if not stations is None:
#     if isinstance(stations, str):
#         stations = [stations]
#     present_stations = [sta for sta in stations if sta in obsdf.index.get_level_values('name').unique()]
#     # warning if stationnames are not found
#     if len(list(set(stations) - set(present_stations))) != 0:
#         rint(f'WARNING: The stations: {list(set(stations) - set(present_stations))} not found in the Analysis.df.')

#     # use all stations when no station defined (with warning)
#     if len(present_stations) == 0:
#         print('WARNING: No available stations for subsetting, use all stations instead.')
#         present_stations = obsdf.index.get_level_values('name').unique().to_list()



# # Subset to stations
# obsdf = obsdf.reset_index()
# obsdf = obsdf[obsdf['name'].isin(present_stations)]

# # Get hours for all records
# obsdf['hour'] = obsdf.reset_index()['datetime'].dt.hour

# # groupby and take the mean per station per hour.
# # hourly_avg = obsdf.groupby(['name', 'hour'])[obstype].mean().unstack().transpose()
# hourly_avg = obsdf.groupby(['name', 'hour'])[obstype].agg(['mean', 'std', 'median'])
# # .unstack().transpose()







# hourly_avg.plot()





















#%%

def make_diurnal_cycle(analysis, stations=None, obstype='temp',
                       relative=False, refstation=None):

    # filter to stations and obstype
    if stations==None:
        stations = analysis.df.reset_index().name.unique()
    stations = set(stations)
    stations = list(stations)
    for stat in stations:
        if not isinstance(stat, str):
            print(f'{stat} is not in the station list, therefore the first station is chosen')
            return
    df = analysis._subset_stations(stations)

    if obstype not in df.columns.values:
        print('Input a valid variable for the plotting')
        return

    # if relative is true we calculate the relative values compared to the reff station
    if relative == True:
        if isinstance(refstation, str):
            if refstation not in stations:
                df_ref=analysis._subset_stations([refstation])
            else:
                print('the refstation need to be different from the other stations')
                return
        else:
            print('chose one valid reff station')
            return
        df = df[obstype]
        df = df.reset_index()
        df_ref = df_ref[obstype]
        df_ref = df_ref.reset_index()
        df_ref = df_ref.rename(columns={obstype:refstation})
        df_ref['datetime'] = pd.to_datetime(df_ref['datetime'])
        df_ref = df_ref.set_index('datetime')
        df_ref = df_ref.resample('60min').first()
        df_ref = df_ref.reset_index()
        df_ref.pop('name')
        for stat in stations:
            df_temp = df[df['name']==stat]
            df_temp = df_temp.rename(columns={obstype:stat})
            df_temp['datetime'] = pd.to_datetime(df_temp['datetime'])
            df_temp = df_temp.set_index('datetime')
            df_temp = df_temp.resample('60min').first()
            df_temp = df_temp.reset_index()
            df_temp.pop('name')
            df_ref = df_ref.merge(df_temp,how='left',left_on='datetime',  right_on='datetime')
        df_total = df_ref
        df_total['hour'] = pd.to_datetime(df_total['datetime'])
        df_total['hour'] = df_total['hour'].dt.hour
        df_total.pop('datetime')
        for stat in stations:
            df_total[stat] = df_total[stat]-df_total[refstation]
        df_total.pop(refstation)
        df_final = pd.DataFrame(columns=df_total.columns, index=range(0,24))
        for h in range(0,24):
            df_final.loc[h] = df_total[df_total['hour']==h].mean()
        df_final=df_final.set_index('hour')
        return df_final
    #if relative is false we only plot the desired stations
    elif relative==False:
        df = df[obstype]
        df = df.reset_index()
        stat_init = stations[0]
        print(stat_init)
        df_total = df[df['name']==stat_init]
        df_total = df_total.rename(columns={obstype:stat_init})
        df_total['datetime'] = pd.to_datetime(df_total['datetime'])
        df_total = df_total.set_index('datetime')
        df_total = df_total.resample('60min').first()
        df_total = df_total.reset_index()
        df_total.pop('name')
        if len(stations)!=1:
            for stat in stations[1:]:
                df_temp = df[df['name']==stat]
                df_temp = df_temp.rename(columns={obstype:stat})
                df_temp['datetime'] = pd.to_datetime(df_temp['datetime'])
                df_temp = df_temp.set_index('datetime')
                df_temp = df_temp.resample('60min').first()
                df_temp = df_temp.reset_index()
                df_temp.pop('name')
                df_total = df_total.merge(df_temp,how='left',left_on='datetime',  right_on='datetime')
        df_total['hour'] = pd.to_datetime(df_total['datetime'])
        df_total['hour'] = df_total['hour'].dt.hour
        df_total.pop('datetime')
        df_final = pd.DataFrame(columns=df_total.columns, index=range(0,24))
        for h in range(0,24):
            df_final.loc[h] = df_total[df_total['hour']==h].mean()
        df_final=df_final.set_index('hour')
        return df_final
    else:
        print('Input needs to be True or False for relative')
        return

