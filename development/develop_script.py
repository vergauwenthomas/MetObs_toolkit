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


lib_folder = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(lib_folder))

# add testdata paths
from tests.push_test.test_data_paths import testdata




tmp_pickle=os.path.join(lib_folder, 'development', 'tmp', 'dev_pickle.pkl')

import metobs_toolkit



#%%

# use_dataset = 'debug_wide'
use_dataset = 'single_netatmo_sara_station'

dataset = metobs_toolkit.Dataset()


dataset.update_settings(output_folder=None,
                        input_data_file=testdata[use_dataset]['datafile'],
                        input_metadata_file=testdata[use_dataset]['metadatafile'],
                        data_template_file=testdata[use_dataset]['template'],
                        metadata_template_file=testdata[use_dataset]['template'],
                        )


dataset.import_data_from_file(**testdata[use_dataset]['kwargs'])

# dataset.coarsen_time_resolution(freq ='60T')
# dataset.apply_quality_control()
#%%

print(dataset)
dataset.make_plot(title='Netatmo Sara, pure input freq')
print(f'frequencies: {dataset.metadf["dataset_resolution"]}')
# dataset.make_plot()
dataset.sync_observations(tollerance='2T')
dataset.make_plot(title='Netatmo Sara, syncronize with tollerance = 2min')


dataset.sync_observations(tollerance='3T')
dataset.make_plot(title='Netatmo Sara, syncronize with tollerance = 3min')



print(dataset)


#%%
import pandas as pd
from metobs_toolkit.data_import import import_data_from_csv, import_metadata_from_csv


use_dataset = 'single_netatmo_sara_station'

dataset = metobs_toolkit.Dataset()


dataset.update_settings(output_folder=None,
                        input_data_file=testdata[use_dataset]['datafile'],
                        input_metadata_file=testdata[use_dataset]['metadatafile'],
                        data_template_file=testdata[use_dataset]['template'],
                        metadata_template_file=testdata[use_dataset]['template'],
                        )


long_format=True
obstype=None,
obstype_unit=None
obstype_description=None
kwargs_data_read={}


df, template = import_data_from_csv(
    input_file=dataset.settings.IO["input_data_file"],
    template_file=dataset.settings.templates["data_template_file"],
    long_format=long_format,
    obstype=obstype,  # only relevant in wide format
    obstype_units = obstype_unit, # only relevant in wide format
    obstype_description = obstype_description, # only relevant in wide format
    kwargs_data_read = kwargs_data_read
)


# Set timezone information
df.index = df.index.tz_localize(
    tz=dataset.settings.time_settings["timezone"],
    ambiguous="infer",
    nonexistent="shift_forward",
)



# drop Nat datetimes if present
df = df.loc[pd.notnull(df.index)]









meta_df = import_metadata_from_csv(
    input_file=dataset.settings.IO["input_metadata_file"],
    template_file=dataset.settings.templates["metadata_template_file"],
    kwargs_metadata_read = {},
)



# in dataset of one station, the name is most often not present!
if not "name" in df.columns:
    # logger.warning(f'No station names find in the observations!' )

    # If there is ONE name in the metadf, than we use that name for
    # the df, else we use the default name
    if (('name' in meta_df.columns) & (meta_df.shape[0] == 1)):
        name = meta_df['name'].iloc[0]
        df['name'] = name
        # logger.warning(f'One stationname found in the metadata: {name}, this name is used for the data.')
    else:
        df["name"] = str(dataset.settings.app["default_name"])
        # for later merging, we add the name column with the default
        # also in the metadf
        meta_df['name'] =str(dataset.settings.app["default_name"])
        # logger.warning(
            # f'Assume the dataset is for ONE station with the \
            # default name: {self.settings.app["default_name"]}.')




# merge additional metadata to observations
meta_cols = [
    colname for colname in meta_df.columns if not colname.startswith("_")
]
additional_meta_cols = list(set(meta_cols).difference(df.columns))

if bool(additional_meta_cols):

    additional_meta_cols.append("name")  # merging on name
    # merge deletes datetime index somehow? so add it back.
    df_index = df.index
    df = df.merge(
        right=meta_df[additional_meta_cols], how="left", on="name"
    )
    df.index = df_index


# convert dataframe to multiindex (datetime - name)
df = df.set_index(["name", df.index])

#%%
from metobs_toolkit.df_helpers import get_freqency_series

freq_estimation_method='median'
freq_estimation_simplify='False'
freq_estimation_simplify_error='1T'





dataset._initiate_df_attribute(dataframe=df, update_metadf=True)
dataset._apply_qc_on_import()

fixed_freq_series = None
freq_series = get_freqency_series(
    df=dataset.df,
    method=freq_estimation_method,
    simplify=freq_estimation_simplify,
    max_simplify_error=freq_estimation_simplify_error,
)

freq_series_import = freq_series

#%% try timeseries cleanup







