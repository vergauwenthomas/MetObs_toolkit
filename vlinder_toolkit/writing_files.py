#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 15:30:55 2023

@author: thoverga
"""

import os
import pandas as pd
from datetime import datetime
from .settings import Settings





def write_dataset_to_csv(df, outliersdf, metadf, 
                         observation_types, location_info,
                         filename, include_outliers):
    """
        Write the dataset to a file where the observations, metadata and (if available)
        the quality labels per observation type are merged together. 
        
        A final qualty controll label for each quality-controlled-observation type
        can be added in the outputfile.
        
        The file will be writen to the Settings.outputfolder.

        Parameters
        ----------
        filename : string, optional
            The name of the output csv file. If none, a standard-filename is generated
            based on the period of data. The default is None.
        add_final_labels : Bool, optional
            If True, a final qualty control label per observation type
            is added as a column. The default is True.

        Returns
        -------
        None

        """
    
    
     
   
    #Get observations and metadata columns in the right order
    # logger.debug('Merging data and metadata')
    
    
    #make column ordering
    df_columns = observation_types.copy() #observations
    df_columns.extend(location_info) #metadata
    qc_columns = [col for col in outliersdf if col.endswith('_label')] #add qc labels
    df_columns.extend(qc_columns)
    df_columns.insert(0, 'datetime') # timestamp as first column
    df_columns.insert(1, 'name') #station name as second column
    
    

    
    #unstack observations and merge with metadf
    df[qc_columns] = 'ok'
    df = pd.concat([df, outliersdf])
    df = df.reset_index()
    
    metadf = metadf.reset_index()
    df = df.merge(metadf, how='left', on='name')
    
    #sort and subset columns
    df = df[df_columns]
    
            
    #find observation type that are not present
    ignore_obstypes = [col for col in observation_types if df[col].isnull().all()]
    
    df = df.drop(columns=ignore_obstypes)
    
    # logger.debug(f'Skip quality labels for obstypes: {ignore_obstypes}.')
    
    df = df.sort_values(['name', 'datetime'])
   
    #make filename
    if isinstance(filename, type(None)):
        startstr = df.index.min().strftime('%Y%m%d') 
        endstr = df.index.max().strftime('%Y%m%d') 
        filename= 'dataset_' + startstr + '_' + endstr
    else:
        if filename.endswith('.csv'):
            filename = filename[:-4] #to avoid two times .csv.csv
        
    filepath = os.path.join(Settings.output_folder, filename + '.csv')
    
    #write to csv in output folder
    # logger.info(f'write dataset to file: {filepath}')
    df.to_csv(path_or_buf=filepath,
                   sep=';',
                   na_rep='NaN',
                   index=True)    