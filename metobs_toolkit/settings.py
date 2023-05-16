#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
All needed setting are combined in a settings class.

@author: thoverga
"""
import os
import json
import shutil
from pathlib import Path
import logging
from pytz import all_timezones, common_timezones

# connect to logger
logger = logging.getLogger(__name__)


class Settings:
    # make settingsfiles path
    _settings_files_path = os.path.join(str(Path(__file__).parent), "settings_files")

    def __init__(self):
        logger.info("Initialising settings")

        # define thematics in settings. Corresponds to settings files.
        self.db = {}
        self.time_settings = {}
        self.app = {}
        self.qc = {}
        self.gap = {}
        self.missing_obs = {}
        self.templates = {}
        self.gee = {}
        self.IO = {
            "output_folder": None,
            "input_data_file": None,
            "input_metadata_file": None,
        }

        # Update (instance and class variables) what can be updated by setingsfiles
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
        """
        Update the database settings of self using the default settings templates
        and the 'db_user' and 'db_passw' envrionment variables if available.
        :return: No return
        :rtype: No return
        """
        logger.debug("Updating Database settings.")
        f = open(os.path.join(Settings._settings_files_path, "server_login.json"))
        login_data = json.load(f)
        f.close()

        self.db["db_host"] = login_data["host"]

        # self.db_host = Settings.db_host
        self.db["db_database"] = login_data["database"]
        self.db["db_obs_table"] = login_data["obs_table"]
        self.db["db_meta_table"] = login_data["meta_table"]

        self.db["db_user"] = os.getenv("VLINDER_DB_USER_NAME")
        self.db["db_passw"] = os.getenv("VLINDER_DB_USER_PASW")

        # import db templates
        from .data_templates.db_templates import (
            vlinder_metadata_db_template,
            vlinder_observations_db_template,
        )

        self.db["vlinder_db_meta_template"] = vlinder_metadata_db_template
        self.db["vlinder_db_obs_template"] = vlinder_observations_db_template

    def _update_time_res_settings(self):
        """
        Update settings on time resolutions of self using the default settings templates.

        :return: No return
        :rtype: No return
        """
        logger.debug("Updating time resolution settings.")
        f = open(
            os.path.join(
                Settings._settings_files_path, "dataset_resolution_settings.json"
            )
        )
        res_settings = json.load(f)
        f.close()

        self.time_settings["target_time_res"] = res_settings["target_time_resolution"]
        self.time_settings["resample_method"] = res_settings["method"]
        self.time_settings["resample_limit"] = res_settings["limit"]
        self.time_settings["timezone"] = res_settings["timezone"]

        # Freq estimation
        self.time_settings['freq_estimation_method'] = res_settings["freq_estimation_method"]
        self.time_settings['freq_estimation_simplify'] = bool(res_settings["freq_estimation_simplify"])
        self.time_settings['freq_estimation_simplify_error'] = res_settings["freq_estimation_simplify_error"]


    def _update_app_settings(self):
        """
        Update prefered display, print, plot and staticinfo settings of self using the default settings templates.
        :return: No return
        :rtype: No return
        """
        logger.debug("Updating app settings.")
        from .settings_files.default_formats_settings import (
            plot_settings,
            print_settings,
            vars_display,
            default_name,
        )
        from .settings_files.default_formats_settings import (
            static_fields,
            categorical_fields,
            location_info,
        )

        # 1. Print settings
        self.app["print_fmt_datetime"] = print_settings["fmt_datetime"]
        self.app["print_max_n"] = int(print_settings["max_print_per_line"])
        # 2. Plot settings
        self.app["plot_settings"] = plot_settings
        self.app["world_boundary_map"] = os.path.join(
            Settings._settings_files_path,
            "world_boundaries",
            "WB_countries_Admin0_10m.shp",
        )

        # 3. display name mappers
        self.app["display_name_mapper"] = vars_display

        # 4 Fields settings
        # fields without timeevolution
        self.app["static_fields"] = static_fields
        self.app["categorical_fields"] = categorical_fields  # wind and lcz
        self.app["location_info"] = location_info  # all possible metadata

        # 5. default name (when station name is not present in dataset)
        self.app["default_name"] = default_name

    def _update_qc_settings(self):
        """
        Update quality control settings of self using the default settings templates.
        :return: No return
        :rtype: No return
        """
        logger.debug("Updating QC settings.")
        from .settings_files.qc_settings import check_settings, checks_info

        self.qc["qc_check_settings"] = check_settings
        self.qc["qc_checks_info"] = checks_info

    def _update_gap_settings(self):
        """
        Update gap defenition and fill settings of self using the default settings templates.
        :return: No return
        :rtype: No return
        """
        logger.debug("Updating gap settings.")
        from .settings_files.gaps_and_missing_settings import (
            gaps_settings,
            gaps_info,
            gaps_fill_settings,
            gaps_fill_info,
            missing_obs_fill_settings,
            missing_obs_fill_info
        )

        self.gap["gaps_settings"] = gaps_settings
        self.gap["gaps_info"] = gaps_info
        self.gap["gaps_fill_settings"] = gaps_fill_settings
        self.gap["gaps_fill_info"] = gaps_fill_info

        self.missing_obs['missing_obs_fill_settings'] = missing_obs_fill_settings
        self.missing_obs['missing_obs_fill_info'] = missing_obs_fill_info

    def _update_templates(self):
        """
        Import the default mapper-template, and assign it to be used on the
        observations and metadata.
        :return: No return
        :rtype: No return

        """
        logger.debug("Updating data templates settings.")
        from .data_templates.import_templates import default_template_file

        # Set default templates
        self.templates["data_template_file"] = default_template_file
        self.templates["metadata_template_file"] = default_template_file

    def _update_gee_settings(self):
        """
        Update the google earth enginge settings using the default settings templates.
        :return: No return
        :rtype: No return
        """
        logger.debug("Updating gee settings.")
        from .settings_files.gee_settings import gee_datasets

        self.gee["gee_dataset_info"] = gee_datasets

    def update_timezone(self, timezonestr):
        """
        Change the timezone of the input data.

        :param timezonestr: Timezone string of the input observations.
        :type timezonestr: String
        :return: None
        :rtype: None
        """
        if not timezonestr in all_timezones:
            print(
                f"timezone: {timezonestr}, is not a valid timezone. Select one of the following:"
            )
            print(f"{common_timezones}")
            return
        else:
            print(
                f'Update timezone: {self.time_settings["timezone"]} --> {timezonestr}'
            )
            self.time_settings["timezone"] = timezonestr

    def update_IO(
        self,
        output_folder=None,
        input_data_file=None,
        input_metadata_file=None,
        data_template_file=None,
        metadata_template_file=None,
    ):
        """
        Update some settings that are relevent before data is imported. The self
        object will be updated.

        :param output_folder: a directory to store the output to, defaults to None.
        :type output_folder: String, optional
        :param input_data_file: Path to the input data file, defaults to None.
        :type input_data_file: String, optional
        :param input_metadata_file: Path to the input metadata file, defaults to None
        :type input_metadata_file: String, optional
        :param data_template_file: Path to the mapper-template csv file to be used on the observations. If not given, the default template is used.
        :type data_template_file: String, optional
        :param metadata_template_file: Path to the mapper-template csv file to be used on the metadata. If not given, the default template is used.
        :type metadata_template_file: String, optional
        :return: No return
        :rtype: No return

        """

        logger.info("Updating settings with input: ")

        if not isinstance(output_folder, type(None)):
            print(
                "Update output_folder: ",
                self.IO["output_folder"],
                " --> ",
                output_folder,
            )
            logger.info(
                f'Update output_folder:  {self.IO["output_folder"]}  -->  {output_folder}'
            )
            self.IO["output_folder"] = output_folder

        if not isinstance(input_data_file, type(None)):
            print(
                "Update input_data_file: ",
                self.IO["input_data_file"],
                " --> ",
                input_data_file,
            )
            logger.info(
                f'Update input_data_file:  {self.IO["input_data_file"]}  -->  {input_data_file}'
            )
            self.IO["input_data_file"] = input_data_file

        if not isinstance(input_metadata_file, type(None)):
            print(
                "Update input_metadata_file: ",
                self.IO["input_metadata_file"],
                " --> ",
                input_metadata_file,
            )
            logger.info(
                f'Update meta_data_file:  {self.IO["input_metadata_file"]}  -->  {input_metadata_file}'
            )
            self.IO["input_metadata_file"] = input_metadata_file

        if not isinstance(data_template_file, type(None)):
            print(
                "Update data template file: ",
                self.templates["data_template_file"],
                " --> ",
                data_template_file,
            )
            logger.info(
                f'Update data template file:  {self.templates["data_template_file"]}  -->  {data_template_file}'
            )
            self.templates["data_template_file"] = data_template_file

        if not isinstance(metadata_template_file, type(None)):
            print(
                "Update metadata template file: ",
                self.templates["metadata_template_file"],
                " --> ",
                metadata_template_file,
            )
            logger.info(
                f'Update metadata template file:  {self.templates["metadata_template_file"]}  -->  {metadata_template_file}'
            )
            self.templates["metadata_template_file"] = metadata_template_file

    def copy_template_csv_files(self, target_folder):
        """
        A function to copy the default template file to an other location. This
        can be of use when creating a template file to start from the default.

        :param target_folder: Directory to copy the default template to (default_template.csv).
        :type target_folder: String
        :return: No return
        :rtype: No return
        """

        from .data_templates.import_templates import default_template_file

        # test if target_folder is a folder
        assert os.path.isdir(target_folder), f"{target_folder} is not a folder"

        target_file = os.path.join(target_folder, "default_template.csv")

        shutil.copy2(default_template_file, target_file)

        print("Templates copied to : ", target_file)

    # =============================================================================
    #     Check settings
    # =============================================================================

    def show(self):
        """
        A function that prints out all the settings, structured per thematic.
        :return: No return
        :rtype: No return

        """
        logger.info("Show settings.")
        class_vars_name = [
            attr
            for attr in dir(Settings)
            if not callable(getattr(Settings, attr)) and not attr.startswith("__")
        ]

        attr_list = [
            "IO",
            "db",
            "time_settings",
            "app",
            "qc",
            "gap",
            "templates",
            "gee",
        ]

        # Drop variables starting with _
        class_vars_name = [mem for mem in class_vars_name if not mem.startswith("_")]
        print("All settings:")
        print(" \n ---------------------------------------\n")

        for theme in attr_list:
            print(f" ---------------- {theme} ----------------------\n")
            printdict = getattr(self, theme)
            for key1, item1 in printdict.items():
                print(f"* {key1}: \n")
                if isinstance(item1, type({})):
                    # nested dict level 1
                    for key2, item2 in item1.items():
                        print(f"  - {key2}: \n")
                        print(f"    -{item2} \n")
                else:
                    # not nested
                    print(f"  -{item1} \n")
