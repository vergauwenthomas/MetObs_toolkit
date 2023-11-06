#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Introduction example script

@author: thoverga
"""

# This is an introduction to get started with the MetObs toolkit. Once the MetObs toolkit package is installed,
# you can import its functionality by:

import metobs_toolkit

# =============================================================================
# Dataset
# =============================================================================

# A dataset is a collection of all observational data. Most of the methods are
# applied directly on a dataset. Start by creating an empty dataset object:

your_dataset = metobs_toolkit.Dataset()

# The most relevant attributes of a Dataset are:
#    * .df --> a pandas DataFrame where all the observation data is stored
#    * .metadf --> a pandas DataFrame where all the metadata for each station is stored
#    * .settings --> a Settings object to store all specific settings.
#    * .missing_obs and .gaps --> here the missing records and gaps are stored if present.

# Note that each Dataset will be equipped with the default settings.


# We created a dataset and stored in under the variable 'aug_2020_all_vlinder'.
# The show function prints out an overview of data in the dataset:
your_dataset.show() # or .get_info()


# TIP: to get an extensive overview of an object, call the .show() method on it.


# =============================================================================
# Importing data
# =============================================================================


# To import your data into a Dataframe, the following files are required:

#    * data file: This is the csv file containing the observations
#    * (optional) metadata file: The csv file containing metadata for each station.
#    * template file: This is a csv file that is used to interpret your data file, and metadata file if present.

# In practice you need to start by creating a template file for your data. More
# information on the creation of the template can be found in documentation
# (under "Mapping to the toolkit").

# TIP: Use the template assistant of the toolkit for creating a template file by using:

# metobs_toolkit.build_template_prompt()



# To import data, you must specify the paths to your data, metadata and template are.
# For this example we use the demo data, metadata and template that comes with
# the toolkit.

your_dataset.update_settings(
    input_data_file=metobs_toolkit.demo_datafile, # path to the data file
    input_metadata_file=metobs_toolkit.demo_metadatafile,
    template_file=metobs_toolkit.demo_template,
)

# The settings of your dataset are updated with the required paths. Now the data can
# be imported into your empty Dataset.
your_dataset.import_data_from_file()

#%%
# =============================================================================
# Inspecting the Data
# =============================================================================

# To get an overview of the data stored in your dataset you can use
your_dataset.show()


# If you want to inspect the data in your dataset directly, you can take a look at the .df attribute
print(your_dataset.df.head())
# equivalent for the metadata
print(your_dataset.metadf.head())


# Generationg timeseries of a Dataset is easily done by:
your_dataset.make_plot(obstype='temp') #This is a basic timeseries plot of the temperatures


# =============================================================================
# Inspecting one station
# =============================================================================

# if you are interested in one station, you can extract all the info for that one station from the dataset by:

favorite_station = your_dataset.get_station(stationname="vlinder02")

#Favorite station now contains all the info of that one station. All methods
#that are applicable on a Dataset are also applicable on a Station. So to inspect your favorite station you can:

print(favorite_station.df.head())

# Or you can make a timeseries plot for a field of choice:
favorite_station.make_plot(obstype="humidity",
                           colorby='label') #colors indicate the label of an observation.
