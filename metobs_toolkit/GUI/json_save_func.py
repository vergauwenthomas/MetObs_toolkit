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
def update_json_file(pydict, catkey=None):
    with open(saved_values_file,'r+') as file:
          # First we load existing data into a dict.
        file_data = json.load(file)

        if isinstance(catkey, type(None)):
            file_data.update(pydict)
        else:
            # Join new_data with file_data inside catkey-dict
            file_data[catkey].append(pydict)

        # Sets file's current position at offset.
        file.seek(0)
        # convert back to json.
        json.dump(file_data, file, indent = 4)


def read_json(jsonfilename):
    with open(jsonfilename,'r') as file:
          # First we load existing data into a dict.
        file_data = json.load(file)

    return file_data


def get_saved_vals():
    vals = read_json(saved_values_file)
    vals = {key: val for key, val in vals.items() if val != ""}
    return vals

