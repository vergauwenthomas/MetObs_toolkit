#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Example script on GEE interactions.

@author: thoverga
"""
import metobs_toolkit



# =============================================================================
# Import data (See previous example script for explanations)
# =============================================================================
your_dataset = metobs_toolkit.Dataset()
your_dataset.update_settings(
    input_data_file=metobs_toolkit.demo_datafile, # path to the data file
    input_metadata_file=metobs_toolkit.demo_metadatafile,
    template_file=metobs_toolkit.demo_template,
)

your_dataset.import_data_from_file()


# =============================================================================
# Using the Google Earth Engine
# =============================================================================


# The Google Earth Engine (GEE) is used to extract geospatial information at the
# locations of the stations.
# To make use of the GEE, you need to setup a google developers account (free of charge).
# The details and steps are documented at the following page: https://vergauwenthomas.github.io/MetObs_toolkit/gee_authentication.html#

# Follow the steps on the 'Using Google Earth Engine' page and then continue.

# =============================================================================
# Extracting LCZ from GEE
# =============================================================================
# The metadata of your stations is stored in the .metadf attribute of your dataset.

print(your_dataset.metadf)

# In order to extract geospatial information for you stations, the lat and lon (latitude and longitude)
# of your stations must be present in the metadf. If so, than geospatial
# information will be extracted from GEE at these locations.

# To extract the Local Climate Zones (LCZ's) of your stations:

lcz_values = your_dataset.get_lcz()

# The first time, in each session, you are asked to authenticated by Google.
# Select your google account and billing project that you have set up and accept
# the terms of condition.
# NOTE: For small datarequest the read-only scopes are sufficient, for large
# datarequests this is insufficient because the data will be written directly to your google Drive.


# The LCZ's for all your stations are extracted
print(lcz_values)
# and the lcz column in your metadata of your dataset is updated
print(your_dataset.metadf['lcz'].head())

# To make a geospatial plot you can use the following method:
your_dataset.make_geo_plot(variable="lcz")

# =============================================================================
# Other geospatial info
# =============================================================================

# Similar as LCZ extraction you can extract the altitude of the stations (from
# a digital elevation model):

altitudes = your_dataset.get_altitude() #The altitudes are in meters above sea level.


# A more detailed description of the landcover/land use in the microenvironment
# can be extracted in the form of landcover fractions in a circular buffer for each station.

# You can select to aggregate the landcoverclasses to water - pervious and impervious,
# or set aggregation to false to extract the landcoverclasses as present in the worldcover_10m dataset.

aggregated_landcover = your_dataset.get_landcover(
                                        buffers=[100, 250], # a list of bufferradii in meters
                                        aggregate=True #if True, aggregate landcover classes to the water, pervious and impervious.
                                        )

print(aggregated_landcover)

# =============================================================================
# Interactive plotting a GEE dataset
# =============================================================================

# You can make an interactive spatial plot to visualize the stations spatially.
# This is done by creating an html file, that can be opened in your browser.

# First specify the location you want to save the plot
import os
spatial_plot_file= os.path.join(os.getcwd(), 'example_spatial_plot.html')
print(spatial_plot_file)


# To make an interactive plot of a GEE dataset one can use the following function:
your_dataset.make_gee_plot(gee_map='worldcover',
                           save=True,
                           outputfile =spatial_plot_file
                           )

