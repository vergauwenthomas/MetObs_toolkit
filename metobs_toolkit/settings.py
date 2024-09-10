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
    """Class definition to store all settings."""

    # make settingsfiles path
    _settings_files_path = os.path.join(str(Path(__file__).parent), "settings_files")

    def __init__(self):
        """Initiate the settings."""
        logger.info("Initialising settings")

        # define thematics in settings. Corresponds to settings files.
        self.app = {}
        self.qc = {}
        self.templatefile = None  # filepath
        self.IO = {
            "output_folder": None,
            "input_data_file": None,
            "input_metadata_file": None,
        }

        # Update (instance and class variables) what can be updated by setingsfiles

        self._update_app_settings()
        self._update_qc_settings()

    # =============================================================================
    #     Update settings from files in initialisation
    # =============================================================================

    def _update_app_settings(self):
        """
        Update preferred display, print, plot and static info settings.

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
            label_def,
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
            # checks_info,
            titan_check_settings,
            titan_specific_labeler,
        )

        self.qc["qc_check_settings"] = check_settings
        # self.qc["qc_checks_info"] = checks_info
        self.qc["titan_check_settings"] = titan_check_settings
        self.qc["titan_specific_labeler"] = titan_specific_labeler

    def update_IO(
        self,
        input_data_file=None,
        input_metadata_file=None,
        template_file=None,
    ):
        """Update the paths to the input files.

        This method will set the path to your data file, (metadata file) and
        template file.

         * input_data_file:  The path to your raw observations (CSV)
         * input_metadata_file: (Optional) The path to your metadata file (CSV)
         * template_file: The path to the template file (JSON). (Use the
           `metobs_toolkit.build_template_prompt()` method to create this file.)

        (This should be applied before importing the observations.)

        Parameters
        ----------
        input_data_file : str, optional
            Path to the input data file, defaults to None.
        input_metadata_file : str, optional
            Path to the input metadata file, defaults to None
        template_file : str, optional
            Path to the mapper-template JSON file to be used on the observations
            and metadata. If not given, the default template is used. The
            default is None.

         Returns
         -------
         None.

        """
        logger.info("Updating settings with input: ")

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
                f"Update template file:  {self.templatefile}  -->  {template_file}"
            )
            self.templatefile = str(template_file)

    # =============================================================================
    #     Check settings
    # =============================================================================

    def show(self):
        """Alias of get_info()."""
        self.get_info()

    def get_info(self):
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

        attr_list = ["IO", "app", "qc", "templatefile"]

        # Drop variables starting with _
        class_vars_name = [mem for mem in class_vars_name if not mem.startswith("_")]
        print("All settings:")
        print(" \n ---------------------------------------\n")

        for theme in attr_list:
            print(f" ---------------- {theme} (settings) ----------------------\n")
            printdict = getattr(self, theme)

            # if attribute is a dict
            if isinstance(printdict, dict):
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
            elif isinstance(printdict, str):
                print(printdict)
            else:
                continue


# =============================================================================
# Docstring test
# =============================================================================
if __name__ == "__main__":
    from metobs_toolkit.doctest_fmt import setup_and_run_doctest

    setup_and_run_doctest()
