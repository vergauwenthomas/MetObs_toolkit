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
    output_data_folder = None
   
    #CSV input settings
    input_file = None
    input_csv_template = None
   
    
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
        self.update_app_settings()
        self.update_templates()
        

    # =============================================================================
    #     Update settings
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
    def update_app_settings(self):
        
        #1. Print settings
        f = open(os.path.join(Settings._settings_files_path, "app_print_settings.json"))
        print_settings = json.load(f)
        f.close()
        Settings.print_fmt_datetime = print_settings['fmt_datetime']
        Settings.print_max_n = int(print_settings["max_print_per_line"])
        #2. Plot settings
    
    @classmethod
    def update_templates(self):
       from .data_templates.csv_templates import vlinder_brian_csv_template, csv_templates_list
       from .data_templates.db_templates import vlinder_metadata_db_template, vlinder_observations_db_template
       
       Settings.template_list = csv_templates_list
       Settings.vlinder_csv_template = vlinder_brian_csv_template
       Settings.vlinder_db_meta_template = vlinder_metadata_db_template
       Settings.vlinder_db_obs_template = vlinder_observations_db_template
       
       
       Settings.input_csv_template = Settings.vlinder_csv_template #set vlindertemplate as standard
    
    
    
    @classmethod
    def update_settings(self, output_data_folder=None, input_file=None):
        if not isinstance(output_data_folder, type(None)):    
            print('Update output_data_folder: ', self.output_data_folder, ' --> ', output_data_folder)
            Settings.output_data_folder = output_data_folder
            
        if not isinstance(input_file, type(None)):    
            print('Update input_file: ', self.input_file, ' --> ', input_file)
            Settings.input_file = input_file
        
        print('settings inputfile: ', Settings.input_file)
        
      

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
        
       
        