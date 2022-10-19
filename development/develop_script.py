#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:25:02 2022

@author: thoverga
"""

# import vlinder_toolkit
import os, sys
import pandas as pd
from pathlib import Path


lib_folder = Path(__file__).resolve().parents[1]
sys.path.append(str(lib_folder))
print(str(lib_folder))

from src import vlinder_toolkit

#% Import

testdatafile = os.path.join(str(lib_folder), 'tests', 'test_data',  'testdata_okt.csv')

static_data =  os.path.join(str(lib_folder), 'static_data', 'vlinder_metadata.csv')



#% Setup dataset
settings = vlinder_toolkit.Settings()
settings.update_settings(input_data_file=testdatafile,
                          input_metadata_file=static_data,
                         )


dataset = vlinder_toolkit.Dataset()
dataset.import_data_from_file(coarsen_timeres=False)







sta = dataset.get_station('vlinder05')
#%%

import numpy as np


dataset = vlinder_toolkit.Dataset()
dataset.import_data_from_file(coarsen_timeres=False)
sta = dataset.get_station('vlinder02')






df = sta.df()



dropdt = df.index[5:11]

dropdt2 = df.index[16]

df = df.drop(index=dropdt)
df = df.drop(index=dropdt2)


#extrac observed frequencies
likely_freq = df.index.to_series().diff().value_counts().idxmax()


missing_datetimeindices = pd.date_range(start = df.index.min(),
                     end = df.index.max(),
                     freq=likely_freq).difference(df.index)

missing_df = pd.DataFrame(data=np.nan,
                          index=missing_datetimeindices,
                          columns=df.columns)

df = df.append(missing_df)


df = df.sort_index()




#%%
# =============================================================================
# checks
# =============================================================================

sta = dataset.get_station('vlinder02')

df_init = sta.df()
sta.make_plot(title='init temp')


sta = vlinder_toolkit.qc_checks.duplicate_timestamp(sta)
sta.make_plot(title='after timstamp dub qc')
sta = vlinder_toolkit.qc_checks.gross_value_check(sta)
sta.make_plot(title='after gross value qc')
sta = vlinder_toolkit.qc_checks.persistance(sta)
sta.make_plot(title='after persistance qc')

# df = sta.df()
# sta.make_plot()


# sta.make_plot()



