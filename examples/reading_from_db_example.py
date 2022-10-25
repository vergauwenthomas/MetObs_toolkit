#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 25 13:00:32 2022

@author: thoverga
"""


import os
from pathlib import Path
main_folder = Path(__file__).resolve().parents[1]
testdata_file = os.path.join(str(main_folder), 'tests', 'test_data',  'vlinderdata.csv' )


import vlinder_toolkit

# =============================================================================
# Settings
# =============================================================================

# 1. Initiate settings object. This object contains all settings needed for furthur analysis
settings = vlinder_toolkit.Settings()




# 4. Updating DB settings
# To extract data directly from the database you need an active UGent VNP connection (or from inside the UGent)
# To connect to the database you need a specific login and password. (Contact Thomas (thomas.vergauwen@meteo.be) for this account.)

#The user and password are extracted from the environment-variables: VLINDER_DB_USER_NAME and VLINDER_DB_USER_PASW

#If this is done correctly than the user and password should appear when running settings.show()

settings.show()


# =============================================================================
# Importing from Database
# =============================================================================


# We have to start by first making an (empty) dataset
sept_2022_all_vlinders = vlinder_toolkit.Dataset()


#We need to specify the start-moment and end-moment for this period. To do this we need to import the datetime module (base python module)
from datetime import datetime

#To get the data it is a simple as:
sept_2022_all_vlinders.import_data_from_database(start_datetime=datetime(2022, 9,1), # 2022/09/01 00:00:00
                                                 end_datetime=datetime(2022,9,25,12,45)) #2022/09/25 12:45:00


