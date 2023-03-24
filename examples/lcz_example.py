#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 12:01:35 2022

@author: thoverga
"""
import os
from pathlib import Path
main_folder = Path(__file__).resolve().parents[1]
testdata_file = os.path.join(str(main_folder), 'tests', 'test_data',  'vlinderdata_small.csv' )
metadata = os.path.join(str(main_folder), 'static_data', 'vlinder_metadata.csv')


#import sys
#sys.path.append(str(main_folder))


import vlinder_toolkit


#%%

# =============================================================================
# Import data
# =============================================================================

                         
#1. Importing a dataset containing mulitple different stations is a function in the Dataset class. First we need to initiate a Dataset with a name of choise.
aug_2020_all_vlinders = vlinder_toolkit.Dataset()
aug_2020_all_vlinders.update_settings(input_data_file=testdata_file, #A demo data file, downloaded with brian tool: https://vlinder.ugent.be/vlinderdata/multiple_vlinders.php
                         input_metadata_file=metadata) 
#Coordinates needed to extract LCZ!
aug_2020_all_vlinders.import_data_from_file()


# The metadata is stored in the metadf attribute:
    
print(aug_2020_all_vlinders.metadf.head())

#you can see that there is no LCZ information yet. Als long as the coordinates are present, the lcz can be extracted.




# =============================================================================
#  Get LCZ
# =============================================================================

# 2. To use the LCZ functions, you need a google develeopers account to make use of google earth engine (gee).
# Creating first such an account, thans simply use this function to extract the LCZ for all stations in your metadata

aug_2020_all_vlinders.get_physiography_data(types=['lcz'])

#Now the metadata is updated with lcz information for each station in the 'lcz' column:
print(aug_2020_all_vlinders.metadf.head())
# or 

print(aug_2020_all_vlinders.metadf['lcz'])
    





# =============================================================================
# Analysing LCZ
# =============================================================================
    

#You can recompute the lcz for all stations by calling the get_lcz function on the metadata.


# To make a geospatial map of the LCZ of all stations:
aug_2020_all_vlinders.make_geo_plot(variable='lcz')











