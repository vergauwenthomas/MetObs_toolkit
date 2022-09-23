#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 16:24:06 2022

@author: thoverga
"""

import pandas as pd


from .textmappers import download_cols_to_class_cols_map
from .physical_info import dtypedict, read_datetime_format


def import_data_from_csv(Settings):
    
    assert not isinstance(Settings.input_file, type(None)), "Specify input file in the settings!"    
    df = pd.read_csv(Settings.input_file, sep=';')
    
    
    assert not df.empty, "Dataset is empty!"


    # rename columns to toolkit attriute names
    df = df.rename(columns=download_cols_to_class_cols_map)
    
    
    #format columns
    df = df.astype(dtype=dtypedict)
    
    #create datetime column
    df['datetime'] =pd.to_datetime(df['_date'] +' ' + df['_time'],
                                    format=read_datetime_format) 
    
    #Set datetime index
    df = df.set_index('datetime', drop=True, verify_integrity=False)
    
    
    #drop 'date' and 'time' columns
    df = df.drop(columns=['_date', '_time'])
    

    return df