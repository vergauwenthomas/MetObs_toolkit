#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 12:01:35 2022

@author: thoverga
"""
import os
from pathlib import Path
main_folder = Path(__file__).resolve().parents[1]
testdata_file = os.path.join(str(main_folder), 'tests', 'test_data',  'testdata_okt_with_metadata.csv' )
metadata = os.path.join(str(main_folder), 'static_data', 'vlinder_metadata.csv')

import vlinder_toolkit



#%%


# =============================================================================
# Settings
# =============================================================================





# Settings contains all information about all the paths that will be used. Without proper initiating the settings, a script will not run.

settings = vlinder_toolkit.Settings()


settings.show()


# =============================================================================
# importing data (test)
# =============================================================================

settings.update_settings(input_data_file=testdata_file, #A demo data file, downloaded with brian tool: https://vlinder.ugent.be/vlinderdata/multiple_vlinders.php
                          input_metadata_file=metadata,
                          output_folder='/home/%s/output/' % os.getenv('USER')) #automatically fill username

your_dataset = vlinder_toolkit.Dataset()


#Run following line:
# your_dataset.import_data_from_file(settings) #Rember that you added the 



#You will get an error that no comatible template is found. This means that the package does not understand the column names of the data file.
#You can see this if we open the file and print the column names:
    
import pandas as pd

data = pd.read_csv(settings.input_data_file, sep=';')

print(list(data.columns))

#To import the data we need to add a template to the default templates.


# =============================================================================
#  The templates
# =============================================================================

#Because different networks store data in different formats, there is a need for templates. The templates are defenitions 
# on how to map the data from the csv files (or online databases) to the names that are used in the backhand of this package. 

#When reading a csv file, the package will scan all available templates, and see which one can be applied to the data. 


#To see all templates:

print(settings.template_list)

#This is a list with templates (stored in a dictionary.)

#This list is created from an csv file, where these templates are defined as different tabs in the csv file. To see the templates,
# you can make a copy to a folder of your choise using this function:
    
settings.copy_template_csv_file(target_folder='/home/%s/Desktop' % os.getenv('USER'))


#If you compair the template with the columnnames of the data file than you will see how this mapping is done. 



# =============================================================================
# adding a template
# =============================================================================

# If you want to analyse a new dataset by another network, than you need to add a template yourself. 
# Start by making a copy of the default templates (see above). Do  change the column names of this csv file, but only the cell values.

# If you dataset is one file (observations and meta data in one csv file), like in this example, than you can add the metadata to the end of the template:
    
#example
# varname	template column name	units	description	dtype	format
# name 	    Stationname			                        object	
# 					
# _date	    Datum_dummy			                        object	%Y-%m-%d
# _time	    Tijd (UTC)		 	                        object	%H:%M:%S
# 					
# temp	    Temperatur_test	            Celcius	2m-temperature	float64	       
# lat	    latitudecolumn                              float64                 #ADD metadata
# lon	    longtitudeclumn			                    float64                 #ADD metadat 




#save the file, and add the template to the settings object by specifying the path of the template you created.
settings.add_csv_template(csv_file='/home/%s/Desktop/default_templates.csv' % os.getenv('USER'))

#Now you can import your data from csv file, and when importing, the package will test your template if it can be applied.



# settings.update_settings(input_data_file=testdata_file, 
#                          input_metadata_file=metadata,
#                          output_folder='/home/%s/output/' % os.getenv('USER')) #automatically fill username

# your_dataset = vlinder_toolkit.Dataset()
# your_dataset.import_data_from_file(settings) #Rember that you added the input file in the settings object, this file will be used.









