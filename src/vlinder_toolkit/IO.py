#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 16:24:06 2022

@author: thoverga
"""

import pandas as pd
import sys




def import_data(Settings):
    
    assert not isinstance(Settings.input_file, type(None)), "Specify input file in the settings!"
    
    df = pd.read_csv(Settings.input_file, sep=';')
    
    
    return df