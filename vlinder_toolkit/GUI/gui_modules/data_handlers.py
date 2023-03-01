#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 08:43:44 2023

@author: thoverga
"""

import pandas as pd
import numpy as np
from bokeh.models import ColumnDataSource



def read_csv(filepath):
    common_seperators = [';',',','    ','.']
    
    for sep in common_seperators:
        
        df = pd.read_csv(filepath, sep=sep)
        assert not df.empty, "Dataset is empty!"
        
        if len(df.columns) > 1:
            break
    
    assert len(df.columns) > 1, f'Only one column detected from import using these seperators: {common_seperators}. See if csv template is correct.'
    return df




def init_datatable_with_default(settingslist):
    """ Create the initial source, with default values for the IO page datatable """

    df = pd.DataFrame()
    for setting in settingslist:
        
        #cleanup by selecting the first element if a list value is given (i.g. units)
        clean_setting = {}
        for key, itemdict in setting.items():
            clean_setting[key] = {}
            for column, value in itemdict.items():
                if isinstance(value, list):
                    value = value[0]
                clean_setting[key][column] = value
        
        clean_df = pd.DataFrame(clean_setting).transpose()
        df = pd.concat([df, clean_df])
    
    df['template_column_name'] = np.nan
    df = df.reset_index()
    df = df.rename(columns={'index': 'toolkit_name'})
    
    return ColumnDataSource(df)
    
    

