#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 13:50:17 2023

@author: thoverga
"""
import pandas as pd

from vlinder_toolkit.df_helpers import init_multiindexdf
from vlinder_toolkit.settings import Settings

from vlinder_toolkit.landcover_functions import (connect_to_gee,
                                                 gee_extract_timeseries)

from vlinder_toolkit.convertors import convert_to_toolkit_units


# =============================================================================
# Class Model data (collection of external model data)
# =============================================================================

class Modeldata():
    def __init__(self, modelname):
        self.df = init_multiindexdf
        self.modelname = modelname
        
    def __repr__(self):
        return f'ModelData instance: {self.modelname} model data of {list(self.df.columns)}'
    
    
    
    def get_ERA5_data(self, metadf, startdt, enddt, obstype='temp'):
        mapinfo = Settings.gee_dataset_info['ERA5_hourly']
        # Connect to Gee
        connect_to_gee()        
        # Get data using GEE
        df = gee_extract_timeseries(metadf = metadf,
                                    mapinfo=mapinfo,
                                    startdt=startdt,
                                    enddt=enddt,
                                    obstype=obstype,
                                    latcolname='lat',
                                    loncolname='lon')
        # Convert to toolkit units
        df[obstype], _tlk_unit = convert_to_toolkit_units(data=df[obstype],
                                               data_unit = mapinfo['band_of_use'][obstype]['units'])
        
        
        self.df = df
        self.modelname = 'ERA5_hourly'
    
       
        
    def interpolate_modeldata(self, to_multiidx, obstype='temp'):
        returndf = init_multiindexdf()
        
        recordsdf = init_multiindexdf()
        recordsdf.index = to_multiidx
        # iterate over stations check to avoid extrapolation is done per stations
        for sta in recordsdf.index.get_level_values('name').unique():
            sta_recordsdf = recordsdf.xs(sta, level='name', drop_level=False) 
            sta_moddf = self.df.xs(sta, level='name', drop_level=False)
            
            
            # check if modeldata is will not be extrapolated !
            if min(sta_recordsdf.index.get_level_values('datetime')) < min(sta_moddf.index.get_level_values('datetime')):
                print('Extrapolation')
            if max(sta_recordsdf.index.get_level_values('datetime')) > max(sta_moddf.index.get_level_values('datetime')):
                print('Extrapolation')
        
            # combine model and records
            mergedf = sta_recordsdf.merge(sta_moddf, how='outer', 
                                          left_index=True, right_index=True)
            
            # reset index for time interpolation
            mergedf = mergedf.reset_index().set_index('datetime').sort_index()
            
            # interpolate missing modeldata
            mergedf[obstype].interpolate(method='time',
                                        limit_area='inside',
                                        inplace=True)
            # convert back to multiindex
            mergedf = mergedf.reset_index().set_index(['name', 'datetime']).sort_index()
            #filter only records 
            mergedf = mergedf.loc[sta_recordsdf.index]
        
            returndf = pd.concat([returndf, mergedf])
        return returndf
    