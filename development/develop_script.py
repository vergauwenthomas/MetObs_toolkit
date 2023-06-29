#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:25:02 2022
@author: thoverga
"""

#%%

# import metobs_toolkit
#
import os
import sys
from pathlib import Path
import pandas as pd
import time
import math


lib_folder = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(lib_folder))

# add testdata paths
from tests.push_test.test_data_paths import testdata




tmp_pickle=os.path.join(lib_folder, 'development', 'tmp', 'dev_pickle.pkl')

import metobs_toolkit


#%%
sara_folder = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/testdata_testday/Sara'

<<<<<<< Updated upstream

bigdatafile = os.path.join(sara_folder,'Vlinder_2022.csv')
template =os.path.join(sara_folder, 'bigvlinder_templatefile.csv')
metadata = os.path.join(sara_folder, 'all_vlinders_metadata.csv')


# dataset1:
datafile = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/testdata_testday/chunkdir/chunk_0.csv'


<<<<<<< HEAD

dataset = metobs_toolkit.Dataset()


dataset.update_settings(output_folder=None,
                        input_data_file=datafile,
                        # input_metadata_file=metadata,
                        data_template_file=template,
                        metadata_template_file=template,
                        )


dataset.import_data_from_file()
dataset.coarsen_time_resolution(freq='60T')

datafile2 = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/testdata_testday/chunkdir/chunk_2.csv'

dataset2 = metobs_toolkit.Dataset()


dataset2.update_settings(output_folder=None,
                        input_data_file=datafile2,
                        # input_metadata_file=metadata,
                        data_template_file=template,
                        metadata_template_file=template,
                        )


dataset2.import_data_from_file()

dataset2.coarsen_time_resolution(freq='60T')

# comb = dataset.combine_all_to_obsspace()

comb = dataset + dataset2



test = comb.metadf

#%%









# # # use_dataset = 'debug_wide'
# use_dataset = 'Congo_single_station'

# dataset = metobs_toolkit.Dataset()


# dataset.update_settings(output_folder=None,
#                         input_data_file=testdata[use_dataset]['datafile'],
#                         input_metadata_file=testdata[use_dataset]['metadatafile'],
#                         data_template_file=testdata[use_dataset]['template'],
#                         metadata_template_file=testdata[use_dataset]['template'],
#                         )


# dataset.import_data_from_file(**testdata[use_dataset]['kwargs'])

# dataset.coarsen_time_resolution(freq = testdata[use_dataset]['coarsen'])
# # dataset.apply_quality_control()
# #%%
# dataset.make_geo_plot(boundbox=[7.0, -14, 47.3, 14])
#%%
bigfile = '/home/thoverga/Documents/thesis_studenten/Amber/data/formatted_amsterdam.csv'
outputdir = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/testdata_testday/chunkdir'
from guppy import hpy
h = hpy() #create session

#%%

def n_lines_of_file(filepath):
    # read number of lines without ectually reading it
    with open(filepath, "rbU") as f:
        num_lines = sum(1 for _ in f)
    return num_lines
#%%
# nlines = n_lines_of_file(bigfile)

# print(nlines)



def get_small_subset(datafile, nrows=10, **kwargs):
    return pd.read_csv(datafile, chunksize=nrows+2, **kwargs).get_chunk(nrows)






def format_datetime_from_one_column(df, datetimecolumn, fmt):

    assert datetimecolumn in df.columns, f'{datetimecolumn} not in the df columns: {df.columns}'

    df['datetime'] = pd.to_datetime(df[datetimecolumn], format=fmt)
    return df




def format_datetime_from_two_columns(df, datecolumn, timecolumn, date_fmt, time_fmt):
    assert datecolumn in df.columns, f'{datecolumn} not in the df columns: {df.columns}'
    assert timecolumn in df.columns, f'{timecolumn} not in the df columns: {df.columns}'

    # add them stringwise together
    df['date'] = df[datecolumn].astype(str)
    df['time'] = df[timecolumn].astype(str)

    df['datetime'] = df['date'] + ' ' + df['time']
    fmt = f'{date_fmt} {time_fmt}'

    return format_datetime_from_one_column(df = df,
                                           datetimecolumn='datetime',
                                           fmt = fmt)




# test = get_small_subset(bigfile, nrows=10, sep=';')

# df = format_datetime_from_one_column(test, 'DateTime', '%Y-%m-%d %H:%M:%S')





#%%




def method1():
    # normal reading
    df = pd.read_csv(bigfile)
    print('\n **** method 1 **** \n')
    print(h.heap())
    del df


def method2(datafile, outputdir, chunksize=500000):

    # Calculate the number of chunks
    nlines = n_lines_of_file(datafile)
    n_chunks =math.ceil(nlines/chunksize)
    print(f' By estimate {n_chunks} chunks are needed! ')

    # chunked reading
    print('\n **** method 2 **** \n')
    for i,chunk in enumerate(pd.read_csv(datafile, chunksize=chunksize)):
        print(f' **** Memory usage after chunk {i}(/{n_chunks}) **** ')
        # print(h.heap())
        # print(chunk[0:10])
        outputfile = os.path.join(outputdir, f'chunk_{i}.csv')
=======
# outfile = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/testdata_testday/Sara/Vlinder_gent_2022.csv'
# infile = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/testdata_testday/Sara/Vlinder_2022.csv'



# def get_small_subset(datafile, nrows=10, **kwargs):
#     return pd.read_csv(datafile, chunksize=nrows+2, **kwargs).get_chunk(nrows)
# test = get_small_subset(infile)


# gentlijst = ['vlinder02', 'vlinder01', 'vlinder05', 'vlinder27']

# df = pd.read_csv(infile)
# subdf = df[df['name'].isin(gentlijst)]
>>>>>>> master

        chunk.to_csv(outputfile, index=False)
        # time.sleep(1)
        del chunk
        
    return n_chunks

<<<<<<< HEAD
# method1()
# method2('/home/thoverga/Downloads/Vlinder_2022.csv', 1000000)


# def format_datetime(df, datetimecolumn, fmt, datecolumn=None, timecolumn=None, )


#%% checking the propper functionality of the __add__ function (no accidental duplicates)


# template = metobs_toolkit.demo_template
# inputdata = metobs_toolkit.demo_datafile
# outputdir='D:/github/vlinder_toolkit/development/test'
# n=method2(inputdata, outputdir, 10000)

# chunk0 = 'D:/github/vlinder_toolkit/development/test/chunk_0.csv'
# chunk1 = 'D:/github/vlinder_toolkit/development/test/chunk_1.csv'


# dataset0 = metobs_toolkit.Dataset()
# dataset0.update_settings(output_folder=None,
#                             input_data_file=chunk0,
#                             #input_metadata_file=metadata,
#                             data_template_file=template,
#                             metadata_template_file=template,
#                             )
# dataset0.import_data_from_file()
# dataset0.coarsen_time_resolution(freq='60T')
=======
# subdf.to_csv(outfile)




#%%

# # use_dataset = 'debug_wide'
use_dataset = 'demo'

dataset = metobs_toolkit.Dataset()


dataset.update_settings(output_folder=None,
                        input_data_file=testdata[use_dataset]['datafile'],
                        input_metadata_file=testdata[use_dataset]['metadatafile'],
                        data_template_file=testdata[use_dataset]['template'],
                        metadata_template_file=testdata[use_dataset]['template'],
                        )
>>>>>>> master


# dataset1 = metobs_toolkit.Dataset()
# dataset1.update_settings(output_folder=None,
#                             input_data_file=chunk1,
#                             #input_metadata_file=metadata,
#                             data_template_file=template,
#                             metadata_template_file=template,
#                             )
# dataset1.import_data_from_file()
# dataset1.coarsen_time_resolution(freq='60T')

<<<<<<< HEAD
# dataset0.show()
# dataset1.show()
# dataset0.make_plot(obstype='temp',colorby='label')
# dataset1.make_plot(obstype='temp',colorby='label')
# bigdataset=dataset0.__add__(dataset1, gapsize=None)
# bigdataset.make_plot(obstype='temp',colorby='label')


#if Add worked then there should be no duplicate outliers in the final plot
#%%

# # # use_dataset = 'debug_wide'
# use_dataset = 'Congo_single_station'

# dataset = metobs_toolkit.Dataset()


# dataset.update_settings(output_folder=None,
#                         input_data_file=testdata[use_dataset]['datafile'],
#                         input_metadata_file=testdata[use_dataset]['metadatafile'],
#                         data_template_file=testdata[use_dataset]['template'],
#                         metadata_template_file=testdata[use_dataset]['template'],
#                         )


# dataset.import_data_from_file(**testdata[use_dataset]['kwargs'])
=======
dataset.import_data_from_file(**testdata[use_dataset]['kwargs'])

# dataset.coarsen_time_resolution(freq = testdata[use_dataset]['coarsen'])



#%%
>>>>>>> master

# dataset.coarsen_time_resolution(freq = testdata[use_dataset]['coarsen'])
# # dataset.apply_quality_control()
# #%%
# dataset.make_geo_plot(boundbox=[7.0, -14, 47.3, 14])

<<<<<<< HEAD
#%%
# from datetime import datetime
# lon_min, lat_min, lon_max, lat_max
# extentlist = [2.260609, 49.25, 6.118359, 52.350618]



=======
analysis = dataset.get_analysis()
stats = analysis.get_diurnal_statistics_with_reference(obstype='temp',
                                                       refstation='vlinder01', # define a (rural) reference station of your dataset, insert the name here
                                                       stations=['vlinder02','vlinder27','vlinder28'], # here you can select the stations you want to include, for example: stations=['vlinder01','vlinder02','vlinder25','vlinder27','vlinder28'],
                                                       #if None then all stations are selected
                                                       startdt=None,
                                                       enddt=None,
                                                       plot=True,
                                                       colorby='name',
                                                       errorbands=False, # standard deviation of both reference station and station included
                                                       verbose=False)



#%%
>>>>>>> master

#%%

<<<<<<< HEAD
=======
metobs_toolkit.build_template_prompt()
>>>>>>> Stashed changes
=======
>>>>>>> master
