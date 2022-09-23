#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
All needed setting are combined in a settings class.

@author: thoverga
"""
# import os

class Settings:
    _user = None
    _output_data_folder = None
    
    _input_file = None
    
    
    def __init__(self, user = _user, output_folder = _output_data_folder, input_file= _input_file):
        
        self.output_data_folder = output_folder #where outputfiles are saved
        self.user = user
        self.input_file = input_file
    
       
        
        
        
        # =============================================================================
        #     Update settings
        # =============================================================================
    def update_settings(self, output_data_folder=None, input_file=None):
        if not isinstance(output_data_folder, type(None)):    
            self.output_data_folder = output_data_folder
            print('Update output_data_folder: ', self.output_data_folder, ' --> ', output_data_folder)
            
        if not isinstance(input_file, type(None)):    
            self.input_file = input_file
            print('Update input_file: ', self.input_file, ' --> ', input_file)
        
        
        
      

    # =============================================================================
    #     Check settings
    # =============================================================================
    def check_settings(self):
        for setting_attr in self.__dict__.keys():
            if isinstance(getattr(self,setting_attr), type(None)):
                print('Setting: ', setting_attr, ' is not known! Please specify this.')
                
            
        
    
    def show(self):
        print('User: ', self.user)
        print('Output data folder: ', self.output_data_folder)
        print('Input observations file: ', self.input_file)
        