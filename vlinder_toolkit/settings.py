#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
All needed setting are combined in a settings class.

@author: thoverga
"""
import os
import json
from pathlib import Path
import logging

#connect to logger
logger = logging.getLogger(__name__)
class Settings:

    #These are CLass variables!! Acceseble by station and dataset classes!

    #make settingsfiles path
    _settings_files_path = os.path.join(str(Path(__file__).parent), 'settings_files')





    def __init__(self):
        logger.info('Initialising settings')

        # define thematics in settings. Corresponds to settings files.
        self.db = {}
        self.time_resolution = {}
        self.app = {}
        self.qc = {}
        self.gap = {}
        self.templates = {}
        self.gee = {}
        self.IO = {'output_folder': None,
                   'input_data_file': None,
                   'input_metadata_file': None}

        #Update (instance and class variables) what can be updated by setingsfiles
        self._update_db_settings()
        self._update_time_res_settings()
        self._update_app_settings()
        self._update_qc_settings()
        self._update_gap_settings()
        self._update_templates()
        self._update_gee_settings()



    # =============================================================================
    #     Update settings from files in initialisation
    # =============================================================================

    def _update_db_settings(self):
        logger.debug('Updating Database settings.')
        f = open(os.path.join(Settings._settings_files_path, "server_login.json"))
        login_data = json.load(f)
        f.close()

        self.db['db_host']=login_data['host']

        # self.db_host = Settings.db_host
        self.db['db_database'] = login_data['database']
        self.db['db_obs_table']= login_data['obs_table']
        self.db['db_meta_table'] = login_data['meta_table']

        self.db['db_user'] = os.getenv('VLINDER_DB_USER_NAME')
        self.db['db_passw'] = os.getenv('VLINDER_DB_USER_PASW')

        #import db templates
        from .data_templates.db_templates import vlinder_metadata_db_template, vlinder_observations_db_template
        self.db['vlinder_db_meta_template'] = vlinder_metadata_db_template
        self.db['vlinder_db_obs_template'] = vlinder_observations_db_template

    def _update_time_res_settings(self):
        logger.debug('Updating time resolution settings.')
        f = open(os.path.join(Settings._settings_files_path, "dataset_resolution_settings.json"))
        res_settings = json.load(f)
        f.close()

        self.time_resolution['target_time_res'] = res_settings['target_time_resolution']
        self.time_resolution['resample_method'] = res_settings['method']
        self.time_resolution['resample_limit'] = res_settings['limit']


    def _update_app_settings(self):
        logger.debug('Updating app settings.')
        from .settings_files.default_formats_settings import plot_settings, print_settings, vars_display
        from .settings_files.default_formats_settings import static_fields, categorical_fields, observation_types, location_info


        #1. Print settings
        self.app['print_fmt_datetime'] = print_settings['fmt_datetime']
        self.app['print_max_n'] = int(print_settings["max_print_per_line"])
        #2. Plot settings
        self.app['plot_settings'] = plot_settings
        self.app['world_boundary_map'] = os.path.join(Settings._settings_files_path,
                                                   "world_boundaries",
                                                   "WB_countries_Admin0_10m.shp")

        # 3. display name mappers
        self.app['display_name_mapper'] = vars_display

        #4 Fields settings
        self.app['static_fields'] = static_fields #fields without timeevolution
        self.app['categorical_fields'] = categorical_fields #wind and lcz
        self.app['observation_types'] = observation_types #order of all possible observations
        self.app['location_info'] = location_info #all possible metadata


    def _update_qc_settings(self):
        logger.debug('Updating QC settings.')
        from .settings_files.qc_settings import check_settings, checks_info
        self.qc['qc_check_settings'] = check_settings
        self.qc['qc_checks_info'] = checks_info


    def _update_gap_settings(self):
        logger.debug('Updating gap settings.')
        from .settings_files.gaps_settings import (gaps_settings, gaps_info,
                                                   gaps_fill_settings,
                                                   gaps_fill_info)

        self.gap['gaps_settings'] = gaps_settings
        self.gap['gaps_info'] = gaps_info
        self.gap['gaps_fill_settings'] = gaps_fill_settings
        self.gap['gaps_fill_info'] = gaps_fill_info



    def _update_templates(self):
       logger.debug('Updating data templates settings.')
       from .data_templates.import_templates import csv_templates_list


       #import csv templates
       self.templates['template_list'] = csv_templates_list

       self.templates['input_metadata_template'] = None #force this template to use
       self.templates['input_csv_template'] = None #force this template to use




    def _update_gee_settings(self):
        logger.debug('Updating gee settings.')
        from .settings_files.gee_settings import gee_datasets

        self.gee['gee_dataset_info'] = gee_datasets




    # @classmethod
    def update_IO(self, output_folder=None, input_data_file=None,
                        input_metadata_file=None):

        logger.info('Updating settings with input: ')

        #TODO: posibility to add a default data template + metadata template

        if not isinstance(output_folder, type(None)):
            print('Update output_folder: ', self.IO['output_folder'], ' --> ', output_folder)
            logger.info(f'Update output_folder:  {self.IO["output_folder"]}  -->  {output_folder}')
            self.IO['output_folder'] = output_folder

        if not isinstance(input_data_file, type(None)):
            print('Update input_data_file: ', self.IO['input_data_file'], ' --> ', input_data_file)
            logger.info(f'Update input_data_file:  {self.IO["input_data_file"]}  -->  {input_data_file}')
            self.IO['input_data_file'] = input_data_file

        if not isinstance(input_metadata_file, type(None)):
            print('Update input_metadata_file: ', self.IO['input_metadata_file'], ' --> ', input_metadata_file)
            logger.info(f'Update meta_data_file:  {self.IO["input_metadata_file"]}  -->  {input_metadata_file}')
            self.IO['input_metadata_file'] = input_metadata_file



    def add_csv_template(self, csv_file):
         """
         Add a template-csv to the templates.
         The Settings class will be updated.

         Parameters
         ----------
         csv_file : String
             csv-template file path.

         Returns
         -------
         None.

         """

         logger.info(f'Adding template from csv file: {csv_file}')

         from .data_templates.import_templates import read_csv_template, check_if_templates_are_unique_defined

         template = read_csv_template(csv_file)
         logger.debug(f'Added teplate: {template}')

         self.templates['template_list'].append(template)

         #Check if all templates are still unique
         check_if_templates_are_unique_defined(self.templates['template_list'])

         # force import to use this template
         self.templates['input_csv_template'] = template


    def copy_template_csv_files(self, target_folder):
         import shutil
         from .data_templates.import_templates import csv_templates_dir, find_csv_filenames

         default_files = find_csv_filenames(csv_templates_dir) #list all default templ files

         #copy all defaults to target folder with dummy index
         _idx = 1
         for templ_file in default_files:
             target_file = os.path.join(target_folder, 'default_templates' + str(_idx) + '.csv')

             shutil.copy2(templ_file, target_file)
             _idx += 1
             logger.info(f'Templatates copied to: {target_file}')

         print("Templates copied to : ", target_file)

    # =============================================================================
    #     Check settings
    # =============================================================================



    def show(self):
        logger.info('Show settings.')
        class_vars_name = [attr for attr in dir(Settings) if not callable(getattr(Settings, attr)) and not attr.startswith("__")]


        attr_list = ['IO', 'db', 'time_resolution', 'app', 'qc', 'gap', 'templates', 'gee']

        #Drop variables starting with _
        class_vars_name = [mem for mem in class_vars_name if not mem.startswith('_')]
        print("All settings:")
        print(' \n ---------------------------------------\n')

        for theme in attr_list:
            print(f' ---------------- {theme} ----------------------\n')
            printdict = getattr(self, theme)
            for key1, item1 in printdict.items():
                print(f'* {key1}: \n')
                if isinstance(item1, type({})):
                    # nested dict level 1
                    for key2, item2 in item1.items():
                        print(f'  - {key2}: \n')
                        print(f'    -{item2} \n')
                else:
                    # not nested
                    print(f'  -{item1} \n')




