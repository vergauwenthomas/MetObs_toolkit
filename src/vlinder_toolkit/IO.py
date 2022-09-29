#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 16:24:06 2022

@author: thoverga
"""

import pandas as pd

from .csv_templates import csv_templates


def import_data_from_csv(Settings):
    
    assert not isinstance(Settings.input_file, type(None)), "Specify input file in the settings!"    
    df = pd.read_csv(Settings.input_file, sep=';')
    
    
    assert not df.empty, "Dataset is empty!"


    # import template
    if isinstance(Settings.input_template, type(None)):
        templ = csv_templates.get_template_from_df_columns(df.columns)

    else:
        templ = Settings.input_template





    # rename columns to toolkit attriute names
    df = df.rename(columns=csv_templates.compress_dict(templ, 'varname'))
    
    
    #COnvert template to package-space
    template = csv_templates.template_to_package_space(templ)
    
    
    #format columns
    df = df.astype(dtype=csv_templates.compress_dict(template, 'dtype'))
    
    #create datetime column
    datetime_fmt = template['_date']['fmt'] + ' ' + template['_time']['fmt']
    df['datetime'] =pd.to_datetime(df['_date'] +' ' + df['_time'],
                                    format=datetime_fmt) 
    #TODO implement timezone settings
    
    
    #Set datetime index
    df = df.set_index('datetime', drop=True, verify_integrity=False)
    
    
    #drop 'date' and 'time' columns
    df = df.drop(columns=['_date', '_time'])
    
    
    #Keep only columns as defined in the template
    for column in df.columns:
        if not (column in template.keys()):
            df = df.drop(columns=[column])
    
  
    

    return df