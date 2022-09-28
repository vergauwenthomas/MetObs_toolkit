#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 08:12:02 2022

@author: thoverga
"""


# %%
# All templates or combined in a list, so if the template is not specified, the corresponding template can be found by iteration.

csv_templates_list = [] #Note that the order of elements is of importance.

# =============================================================================
# VLINDER CSV templates 
# =============================================================================


#templates have nested dict structure where the keys are the column names in the csv file, and the 
# values contain the mapping information to the toolkit classes and names. 



vlinder_brian_csv_template = {
        'Vlinder': {'varname': 'name',
                    'dtype': 'object'},
        
         'Datum': {'varname': '_date',
                   'fmt':'%Y-%m-%d',
                   'dtype': 'object' },
         'Tijd (UTC)': {'varname':'_time',
                        'fmt': '%H:%M:%S',
                        'dtype': 'object'},
         'Temperatuur': {'varname':'temp',
                        'units': r'$^o$C',
                        'dtype': 'float64',
                        'description': 'temperature' },
         'Vochtigheid': {'varname':'humidity',
                         'units': '%',
                         'dtype': 'float64',
                         'description': 'relative humidity'},
         'Luchtdruk': {'varname':'pressure', 
                       'units': 'pa',
                       'dtype': 'float64',
                       'description': 'airpressure'},
         'Neerslagintensiteit': {'varname':'precip',
                                 'units': r'l/m$^2$',
                                 'dtype': 'float64',
                                 'description': 'precipitation intensity'},
         'Neerslagsom': {'varname':'precip_sum', 
                         'units': r'l/m^2',
                         'dtype': 'float64',
                         'description': 'precipitation cumulated from midnight'},
         'Windrichting': {'varname': 'wind_direction',
                          'units': r'° from North (CW)',
                          'dtype': 'float64',
                          'description': 'Wind direction'},
         'Windsnelheid': {'varname':'wind_speed',
                          'units': r'm/s',
                          'dtype': 'float64',
                          'description': 'windspeed'},
         'Rukwind': {'varname': 'wind_gust',
                     'units': r'm/s',
                     'dtype': 'float64',
                     'description': 'windgust'},
         'Luchtdruk_Zeeniveau': {'varname': 'pressure_at_sea_level',
                                 'units': 'pa',
                                 'dtype': 'float64',
                                 'description': 'pressure at sea level'},
         'Globe Temperatuur': {'varname': 'radiation_temp',
                               'units': r'celscius denk ik??',
                               'dtype': 'float64',
                               'description': 'Radiative temperature'},
         
        }

csv_templates_list.append(vlinder_brian_csv_template)

def template_to_package_space(template):
    returndict = {}
    
    for orig_fieldname, item in template.items():
        # print(item)
        new_name = item['varname']
        
        
        returndict[new_name] = item
        returndict[new_name].pop('varname')
        returndict[new_nam]['orig_name'] = orig_fieldname #maybe needed for fruther developm.?_
        
        
    return returndict


test = template_to_package_space(vlinder_brian_csv_template)


    
def compress_dict(nested_dict, valuesname):
    """
    This function unnests a nested dictionary for a specific valuename that is a key in the nested dict. 

    Parameters
    ----------
    
    nested_dict : dict 
        Nested dictionary
    
    valuesname : str 
        Nested dict Key-name of nested dict.

    Returns
    -------
    returndict : DICT
        A dictionarry where the keys are kept that have the valuesname as a nesteddict key, 
        and values are the values of the values of the valuesname. 
        {[key-nested_dict-if-exists]: nested_dict[key-nested_dict-if-exists][valuesname]}

    """
    returndict = {}
    for key, item in nested_dict.items():
        if valuesname in item:
            returndict[key] = item[valuesname]
    return returndict    
    

def get_template_from_df_columns(columns, csv_templates_list=csv_templates_list):
    template_found= False
    for csv_template in csv_templates_list:
        columnnames = compress_dict(csv_template, 'varname')
        template_found = [col_name_csv in columnnames for col_name_csv in columns].all()
        if template_found:
            return csv_template
        
    
    
    