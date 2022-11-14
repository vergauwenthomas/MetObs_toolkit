#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
All needed setting are combined in a settings class.

@author: thoverga
"""
import os
import json
from pathlib import Path

class Settings:
    
    #These are CLass variables!! Acceseble by station and dataset classes! 
    
    #make settingsfiles path
    _settings_files_path = os.path.join(str(Path(__file__).parent), 'settings_files') 
    
    #Output 
    output_folder = None
   
    #input data settings
    input_data_file = None
    input_csv_template = None
   
    input_metadata_file = None
    input_metadata_template = None
    
    #Geo datasets templates and info
    geo_datasets_templates = None
    geo_lcz_file = None
    
    #String mappers for display
    display_name_mapper = None
    
    
    #plot settings
    plot_settings = None
    world_boundary_map = None
    
    #time resolution settings
    target_time_res = None # "H"
    resample_method = None #"bfill"
    
    #Database settings (these will be updated in init)
    db_host = None
    db_database = None
    db_obs_table = None
    db_meta_table = None
    db_user = None 
    db_passw = None        

    #App print settings (these will be updated in init)
    print_fmt_datetime = None
    
    #Load templates
    template_list = None #List with available templates
    vlinder_csv_template = None
    vlinder_db_meta_template = None
    vlinder_db_obs_template = None
    
    
    def __init__(self):
     
        #Update (instance and class variables) what can be updated by setingsfiles
        self.update_db_settings() 
        self.update_time_res_settings()
        self.update_app_settings()
        self.update_templates()
        
        

    # =============================================================================
    #     Update settings from files
    # =============================================================================
    @classmethod
    def update_db_settings(self):
        f = open(os.path.join(Settings._settings_files_path, "server_login.json"))
        login_data = json.load(f)
        f.close()
        
        
        Settings.db_host=login_data['host']
        
        # self.db_host = Settings.db_host
        Settings.db_database = login_data['database']
        Settings.db_obs_table = login_data['obs_table']
        Settings.db_meta_table = login_data['meta_table']
        
        Settings.db_user = os.getenv('VLINDER_DB_USER_NAME')
        Settings.db_passw = os.getenv('VLINDER_DB_USER_PASW')
        
    @classmethod
    def update_time_res_settings(self):
        f = open(os.path.join(Settings._settings_files_path, "dataset_resolution_settings.json"))
        res_settings = json.load(f)
        f.close()
        
        Settings.target_time_res = res_settings['target_time_resolution']
        Settings.resample_method = res_settings['method']
       
        
    @classmethod
    def update_app_settings(self):
        from .settings_files.default_formats_settings import plot_settings, print_settings, vars_display
        #1. Print settings
        Settings.print_fmt_datetime = print_settings['fmt_datetime']
        Settings.print_max_n = int(print_settings["max_print_per_line"])
        #2. Plot settings
        
        Settings.plot_settings = plot_settings
        Settings.world_boundary_map = os.path.join(Settings._settings_files_path,
                                                   "world_boundaries",
                                                   "WB_countries_Admin0_10m.shp")
    
        # 3. display name mappers
        Settings.display_name_mapper = vars_display
    
    @classmethod
    def update_templates(self):
       # from .data_templates.csv_templates import vlinder_brian_csv_template, csv_templates_list, vlinder_static_meta_data
       from .data_templates.import_templates import csv_templates_list
       from .data_templates.db_templates import vlinder_metadata_db_template, vlinder_observations_db_template
       from .data_templates.geo_datasets_templates import geo_datasets
       
       #import csv templates
       Settings.template_list = csv_templates_list
       # Settings.vlinder_csv_template = vlinder_brian_csv_template
      
       #import db templates
       Settings.vlinder_db_meta_template = vlinder_metadata_db_template
       Settings.vlinder_db_obs_template = vlinder_observations_db_template
       
       #import geo datasets templates
       Settings.geo_datasets_templates = geo_datasets
       
       
       #Set standard templates
       
       # Settings.input_csv_template = Settings.vlinder_csv_template #set vlindertemplate as standard
       # Settings.input_metadata_template = vlinder_static_meta_data #set vlindertemplate as standard
    
    
    @classmethod
    def update_settings(self, output_folder=None, input_data_file=None,
                        input_metadata_file=None, geotiff_lcz_file=None):
        #TODO: posibility to add a default data template + metadata template
        
        if not isinstance(output_folder, type(None)):    
            print('Update output_folder: ', self.output_folder, ' --> ', output_folder)
            Settings.output_folder = output_folder
            
        if not isinstance(input_data_file, type(None)):    
            print('Update input_data_file: ', self.input_data_file, ' --> ', input_data_file)
            Settings.input_data_file = input_data_file
        
        if not isinstance(input_metadata_file, type(None)):    
            print('Update input_metadata_file: ', self.input_metadata_file, ' --> ', input_metadata_file)
            Settings.input_metadata_file = input_metadata_file
        
        if not isinstance(geotiff_lcz_file, type(None)):    
            print('Update input_metadata_file: ', self.geo_lcz_file, ' --> ', geotiff_lcz_file)
            Settings.geo_lcz_file = geotiff_lcz_file
        
        
        
      

    # =============================================================================
    #     Check settings
    # =============================================================================
    def check_settings(self):
        class_vars_name = [attr for attr in dir(Settings) if not callable(getattr(Settings, attr)) and not attr.startswith("__")]
        
        #Drop variables starting with _
        class_vars_name = [mem for mem in class_vars_name if not mem.startswith('_')]
        
        for setting_attr in class_vars_name:
            if isinstance(getattr(Settings,setting_attr), type(None)):
                print('Setting: ', setting_attr, ' is not known! Please specify this.')
            
        
        
    
    def show(self):
        
        class_vars_name = [attr for attr in dir(Settings) if not callable(getattr(Settings, attr)) and not attr.startswith("__")]
        
        #Drop variables starting with _
        class_vars_name = [mem for mem in class_vars_name if not mem.startswith('_')]
        print("All settings:")
        
        print('')
        print(' ---------------------------------------')
        print('')
        for var in class_vars_name:
            value =  getattr(Settings, var)
            if isinstance(value, str):
                print('  * ', var, '  :  ', value)
            else:
                value=str(value)[0:Settings.print_max_n]
                print('  * ', var, '  :  ', value + '...')
        print('')
        print(' ---------------------------------------')
        
       
        