#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 08:12:02 2022

@author: thoverga
"""
import os, sys
from pathlib import Path
import pandas as pd

#%%


csv_templates_file = os.path.join(str(Path(__file__).parent), 'csv_templates.xlsx')

# %%
# All templates or combined in a list, so if the template is not specified, the corresponding template can be found by iteration.

csv_templates_list = [] #Note that the order of elements is of importance.

# =============================================================================
# VLINDER CSV templates 
# =============================================================================


#templates have nested dict structure where the keys are the column names in the csv file, and the 
# values contain the mapping information to the toolkit classes and names. 



# =============================================================================
# Read the template excel file and convert it to nested dict + add it to the template list
# =============================================================================

templ_dict = pd.read_excel(csv_templates_file, sheet_name=None)




for templ_name in templ_dict.keys():
    data = templ_dict[templ_name]
    #Drop emty rows 
    data = data.dropna(axis='index', how='all')
    
    #Drop variables that are not present in data
    data = data[data['template column name'].notna()]
    
    #create dictionary from dataframe
    data = data.set_index('template column name')
    
    
    
    #create a dict from the dataframe, remove Nan value row wise
    template = {}
    for idx, row in data.iterrows():
        template[idx] = row[~row.isnull()].to_dict()
    
    csv_templates_list.append(template)




# =============================================================================
# add template of static meta data
# =============================================================================

vlinder_static_meta_data = {
    'ID': {'varname': '_ID',
                'dtype': 'object'},
    'VLINDER': {'varname': 'name',
                'dtype': 'object'},
    'lat': {'varname': 'lat',
                'dtype': 'float64'},
    'lon': {'varname': 'lon',
                'dtype': 'float64'},
    'stad': {'varname': 'location',
                'dtype': 'object'},
    'benaming': {'varname': 'call_name',
                'dtype': 'object'},
    'Network': {'varname': 'network',
                'dtype': 'object'},
    }

csv_templates_list.append(vlinder_static_meta_data)




        
# =============================================================================
# Check if templates column names are unique
# =============================================================================
def check_if_templates_are_unique_defined(templatelist):
    templ_df = pd.DataFrame()
    for templ in templatelist:
        templ_df = templ_df.append(pd.Series(templ.keys()), ignore_index=True)
    
    #if all columnnames are identical to other template (type 1)
    if templ_df.duplicated().any():
        sys.exit('Duplicated csv_template column names found (type 1). Make shure the templates have unique column identifiers!')


    #If all columnnames are in a subset of the columnnames of another template (type 2).
    for _idx, row in templ_df.iterrows():
        others = templ_df.loc[templ_df.index != _idx]
        for _idx_others, row_others in others.iterrows():
            counts = row.dropna().isin(row_others.dropna()).value_counts()
            if True in counts.index:
                n_cols_row = len(row.dropna())
                if counts[True] == n_cols_row:
                    sys.exit('Duplicated csv_template column names found (type 2). Make shure the templates have unique column identifiers!')
                    
            

check_if_templates_are_unique_defined(csv_templates_list)