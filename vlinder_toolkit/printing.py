#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 14:56:26 2023

@author: thoverga
"""

import pandas as pd
from datetime import datetime
from .settings import Settings


def print_dataset_info(df, outliersdf, gapsdf):

    if df.empty:
        print("This dataset is empty!")
        # logger.error('The dataset is empty!')
    else: 
        print('\n','--------  General ---------', '\n')
        print(f' .... ')
        
        
        
        print('\n','--------  Observations ---------', '\n')
        starttimestr = datetime.strftime(min(df.index.get_level_values(level='datetime')),
                                         Settings.print_fmt_datetime)
        endtimestr = datetime.strftime(max(df.index.get_level_values(level='datetime')), Settings.print_fmt_datetime)
        
        stations_available = list(df.index.get_level_values(level='name').unique())
        print(f'Observations found for period: {starttimestr} --> {endtimestr}')
        # logger.debug(f'Observations found for period: {starttimestr} --> {endtimestr}')
        print(f'Following stations are in dataset: {stations_available}')
        # logger.debug(f'Following stations are in dataset: {stations_available}')
        
        print('\n', '--------  Outliers ---------', '\n')
        print(f'There are {outliersdf.shape[0]} flagged observations found in total. They occure in these stations: {list(outliersdf.index.get_level_values("name").unique())}')
        
        print('\n', '--------  Gaps ---------', '\n')
        print(f'There are {gapsdf.shape[0]} gaps found in total. They occure in these stations: {list(gapsdf.index.unique())}')
        