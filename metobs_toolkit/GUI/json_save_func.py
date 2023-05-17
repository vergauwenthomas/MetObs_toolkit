#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 10:06:18 2023

@author: thoverga
"""

import os, sys
from pathlib import Path
import json


GUI_folder = Path(__file__).resolve().parents[0]
saved_values_file = os.path.join(str(GUI_folder),
                                 'save_gui_vals.json')




# =============================================================================
# Append to the jsonfile
# =============================================================================

# function to add to JSON
# def update_json_file(pydict, catkey=None):
#     with open(saved_values_file,'r+') as file:
#           # First we load existing data into a dict.
#         file_data = json.load(file)

#         if isinstance(catkey, type(None)):
#             file_data.update(pydict)
#         else:
#             # Join new_data with file_data inside catkey-dict
#             file_data[catkey].append(pydict)

#         # Sets file's current position at offset.
#         file.seek(0)
#         # convert back to json.
#         json.dump(file_data, file, indent = 4)

def update_json_file(update_dict, filepath=None):
    if isinstance(filepath, type(None)):
        filepath=saved_values_file
    # read the existing JSON data from the file or create an empty dict if it doesn't exist
    try:
        with open(filepath, 'r') as f:
            data = json.load(f)
    except FileNotFoundError:
        data = {}

    # update or add to the JSON data using the update_dict
    for key, value in update_dict.items():
        data[key] = value

    # write the updated JSON data to the file
    with open(filepath, 'w') as f:
        json.dump(data, f, indent = 4)



def read_json(jsonfilename):
    with open(jsonfilename,'r') as file:
          # First we load existing data into a dict.
        file_data = json.load(file)

    return file_data





def get_saved_vals():
    vals = read_json(saved_values_file)
    vals = {key: val for key, val in vals.items() if val != ""}
    return vals


#%%


