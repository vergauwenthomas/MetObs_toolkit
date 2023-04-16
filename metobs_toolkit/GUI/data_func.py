#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 11:00:52 2023

@author: thoverga
"""

import os, sys
import pandas as pd



def pass_return(returnvalue, is_error=False,  err_message='Error detected!'):

    return returnvalue, not(is_error), err_message





def isvalidfile(filepath, filetype=None):
    if not os.path.isfile(filepath):
        return False, "f{filepath} is not a file."

    if not isinstance(filetype, type(None)):
        if not filepath.endswith(filetype):
            return False, f'{filepath} is not of the {filetype} format.'

    return True, 'Ok'






def readfile(filepath):
    common_seperators = [';',',','    ','.']

    df=pd.DataFrame()
    for sep in common_seperators:
        try:
            df = pd.read_csv(filepath, sep=sep, engine='python')
        except:
            pass

        if not df.empty:
            break #df found

    if df.empty:
        return pass_return(df, True, f'Datafile is empyt.')
    return pass_return(df)








def get_columns(filepath):

    # check if filepath is file
    _isfile, _err_m = isvalidfile(filepath, filetype='.csv')
    if not _isfile:
        #pass errormessage
        return pass_return([], True, _err_m)


    # read file into dataframe
    df, _ok, _err_m = readfile(filepath)
    if not _ok:
        #pass errormessage
        return pass_return([], True, _err_m)

    return pass_return(list(df.columns))