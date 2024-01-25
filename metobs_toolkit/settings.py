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
    """Class defenition to store all settings."""

    # make settingsfiles path
    _settings_files_path = os.path.join(str(Path(__file__).parent), "settings_files")

    def __init__(self):
        """Initiate the settings."""
        logger.info("Initialising settings")

        # define thematics in settings. Corresponds to settings files.
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
        self._update_time_res_settings()
        self._update_app_settings()
        self._update_qc_settings()
        self._update_gap_settings()
        self._update_templates()
        self._update_gee_settings()

    # =============================================================================
    #     Update settings from files in initialisation
    # =============================================================================

    def _update_time_res_settings(self):
        """
        Update settings on time resolutions of self using the default settings templates.

        Returns
        -------
        None.
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
        self.time_settings["freq_estimation_method"] = res_settings[
            "freq_estimation_method"
        ]
        self.time_settings["freq_estimation_simplify"] = bool(
            res_settings["freq_estimation_simplify"]
        )
        self.time_settings["freq_estimation_simplify_error"] = res_settings[
            "freq_estimation_simplify_error"
        ]

    def _update_app_settings(self):
        """
        Update prefered display, print, plot and staticinfo settings of self using the default settings templates.

        Returns
        -------
        None.
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

        Returns
        -------
        None.
        """
        logger.debug("Updating QC settings.")
        from .settings_files.qc_settings import (
            check_settings,
            checks_info,
            titan_check_settings,
            titan_specific_labeler,
        )

        self.qc["qc_check_settings"] = check_settings
        self.qc["qc_checks_info"] = checks_info
        self.qc["titan_check_settings"] = titan_check_settings
        self.qc["titan_specific_labeler"] = titan_specific_labeler

    def _update_gap_settings(self):
        """
        Update gap defenition and fill settings of self using the default settings templates.

        Returns
        -------
        None.
        """
        logger.debug("Updating gap settings.")
        from .settings_files.gaps_and_missing_settings import (
            gaps_settings,
            gaps_info,
            gaps_fill_settings,
            gaps_fill_info,
            missing_obs_fill_settings,
            missing_obs_fill_info,
        )

        self.gap["gaps_settings"] = gaps_settings
        self.gap["gaps_info"] = gaps_info
        self.gap["gaps_fill_settings"] = gaps_fill_settings
        self.gap["gaps_fill_info"] = gaps_fill_info

        self.missing_obs["missing_obs_fill_settings"] = missing_obs_fill_settings
        self.missing_obs["missing_obs_fill_info"] = missing_obs_fill_info

    def _update_templates(self):
        """
        Import the default mapper-template, and used it on the observations and metadata.

        Returns
        -------
        None.

        """
        logger.debug("Updating data templates settings.")
        from .data_templates.import_templates import default_template_file

        # Set default templates
        self.templates["template_file"] = default_template_file

    def _update_gee_settings(self):
        """
        Update the google earth enginge settings using the default settings templates.

        Returns
        -------
        None.
        """
        logger.debug("Updating gee settings.")
        from .settings_files.gee_settings import gee_datasets

        self.gee["gee_dataset_info"] = gee_datasets

    def update_timezone(self, timezonestr):
        """
        Change the timezone of the input data.

        Parameters
        ------------
        timezonestr : str
            Timezone string of the input observations.

        Returns
        -------
        None.
        """
        if timezonestr not in all_timezones:
            print(
                f"timezone: {timezonestr}, is not a valid timezone. Select one of the following:"
            )
            print(f"{common_timezones}")
            return
        else:
            logger.info(
                f'Update timezone: {self.time_settings["timezone"]} --> {timezonestr}'
            )
            self.time_settings["timezone"] = timezonestr

    def update_IO(
        self,
        output_folder=None,
        input_data_file=None,
        input_metadata_file=None,
        template_file=None,
    ):
        """
        Update some settings that are relevent before data is imported.

        When a argument is None, no update of that settings is performed.
        The self object will be updated.

        Parameters
        ----------
        output_folder : str, optional
            A directory to store the output to, defaults to None.
        input_data_file : str, optional
            Path to the input data file, defaults to None.
        input_metadata_file : str, optional
            Path to the input metadata file, defaults to None
        template_file : str, optional
            Path to the mapper-template csv file to be used on the observations
            and metadata. If not given, the default template is used. The
            default is None.

        Returns
        -------
        None.

        """
        logger.info("Updating settings with input: ")

        if not isinstance(output_folder, type(None)):
            logger.info(
                f'Update output_folder:  {self.IO["output_folder"]}  -->  {output_folder}'
            )
            self.IO["output_folder"] = output_folder

        if not isinstance(input_data_file, type(None)):
            logger.info(
                f'Update input_data_file:  {self.IO["input_data_file"]}  -->  {input_data_file}'
            )
            self.IO["input_data_file"] = input_data_file

        if not isinstance(input_metadata_file, type(None)):
            logger.info(
                f'Update meta_data_file:  {self.IO["input_metadata_file"]}  -->  {input_metadata_file}'
            )
            self.IO["input_metadata_file"] = input_metadata_file

        if not isinstance(template_file, type(None)):
            logger.info(
                f'Update template file:  {self.templates["template_file"]}  -->  {template_file}'
            )
            self.templates["template_file"] = template_file

    def copy_template_csv_files(self, target_folder):
        """Copy the default template.

        A function to copy the default template file to an other location. This
        can be of use when creating a template file to start from the default.

        Parameters
        ----------
        target_folder : str
            Directory to copy the default template to (default_template.csv).

        Returns
        -------
        None.

        """
        from .data_templates.import_templates import default_template_file

        # test if target_folder is a folder
        assert os.path.isdir(target_folder), f"{target_folder} is not a folder"

        target_file = os.path.join(target_folder, "default_template.csv")

        shutil.copy2(default_template_file, target_file)

        logger.info("Templates copied to : ", target_file)

    # =============================================================================
    #     Check settings
    # =============================================================================

    def show(self):
        """Print out an overview of the settings.

        Returns
        -------
        None.

        """
        logger.info("Show settings.")
        class_vars_name = [
            attr
            for attr in dir(Settings)
            if not callable(getattr(Settings, attr)) and not attr.startswith("__")
        ]

        attr_list = [
            "IO",
            "time_settings",
            "app",
            "qc",
            "gap",
            "missing_obs",
            "templates",
            "gee",
        ]

        # Drop variables starting with _
        class_vars_name = [mem for mem in class_vars_name if not mem.startswith("_")]
        print("All settings:")
        print(" \n ---------------------------------------\n")

        for theme in attr_list:
            print(f" ---------------- {theme} (settings) ----------------------\n")
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
