#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Extension of the Dataset class (methods for updating settings).
@author: thoverga
"""
import logging
import os
import pandas as pd

logger = logging.getLogger(__name__)


class DatasetSettingsCore:
    """Extension on the metobs_toolkit.Dataset class with updaters."""

    def update_file_paths(
        self,
        input_data_file=None,
        input_metadata_file=None,
        template_file=None,
    ):
        """Update the paths to the input files.

        This method will set the path to your data file, metadata file and
        template file if provided.

         * input_data_file:  The path to your raw observations (CSV)
         * input_metadata_file: The path to your metadata file (CSV)
         * template_file: The path to the template file (JSON). (Use the
           `metobs_toolkit.build_template_prompt()` method to create this file.)

        (This should be applied before importing the observations.)


        Parameters
        ----------

        input_data_file : string, optional
            Path to the input data file with observations (CSV). If None, the
            path is not updated. The default is None.
        input_metadata_file : string, optional
            Path to the input metadata file (CSV). If None, the
            path is not updated. The default is None.
        template_file : string, optional
            Path to the template file (JSON) to be used on the observations
            and metadata. The default is None.

        Returns
        -------
        None.

        See Also
        -----------
        Dataset.update_output_dir: Update the default output directory.

        Note
        -----
        This method is redundant if you specify the paths in the
        `metobs_toolkit.Dataset.import_data_from_file()` or
        `metobs_toolkit.Dataset.import_only_metadata_from_file()` methods.

        Note
        -------
        To create the template file (JSON), use the `metobs_toolkit.build_template_prompt()` method.

        Warning
        ---------
        In previous versions ( <= v0.2.1) the templatefile was a CSV file. Thus
        you have to create the template again to be compatible with this
        version of the toolkit.

        Examples
        --------
        Start by creating a template file. This file contains all the info for
        the toolkit on how to interpret your raw datafile.

        To create the template run the `metobs_toolkit.build_template_prompt()`
        function.

        >>> import metobs_toolkit
        >>> metobs_toolkit.build_template_prompt() # doctest: +SKIP

        You will be prompted some questions on your data file (and metadatafile if you have one).
        In the end, the toolkit will write the templatefile (JSON format) to
        a location of your choice.

        Now you can create a Dataset (empty), update the paths on it, and import the data.

        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> dataset = update_file_paths(input_data_file=" ... ",  # doctest: +SKIP
        ...                             template_file=" ... ",  # doctest: +SKIP
        ...                             input_metadata_file=" ... ")  # doctest: +SKIP
        >>> dataset.import_data_from_file()  # doctest: +SKIP

        To be more concise, you can avoid the `Dataset.update_file_paths()` method
        and pass the paths directly to the `Dataset.update_file_paths()` method:

        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> dataset.import_data_from_file(input_data_file=" ... ",  # doctest: +SKIP
        ...                               input_metadata_file=" ... ",  # doctest: +SKIP
        ...                               template_file=" ... ")  # doctest: +SKIP

        """
        self.settings.update_IO(
            input_data_file=input_data_file,
            input_metadata_file=input_metadata_file,
            template_file=template_file,
        )

    def update_output_dir(self, dirpath):
        """Update the path to the output directory.

        The output directory is used if an output is created (figure, datafiles, ...),
        without specifying the output folder in the specific method.


        Parameters
        ----------
        dirpath : str
            Path to existing directory that will be used as default output directory.

        Returns
        -------
        None.

        See Also
        -----------
        Dataset.update_file_paths: Update the input file paths

        """

        if not os.path.isdir(dirpath):
            raise MetobsDatasetSettingsUpdaterError(
                f"{dirpath} is not a path to an existing directory."
            )
        self.settings.IO["output_folder"] = str(dirpath)

    def update_default_name(self, default_name):
        """Update the default name (the name of the station).

        This name will be used when only one station is detected. Applying
        this method will overwrite the 'default_name' in the templatefile.

        (All observations are assumed to come from one station.)

        Parameters
        ----------
        default_name : string
            Default name to use when no names are present in the data.

        Returns
        -------
        None.

        """
        self.settings.app["default_name"] = str(default_name)

    def update_qc_settings(
        self,
        obstype="temp",
        dupl_timestamp_keep=None,
        persis_time_win_to_check=None,
        persis_min_num_obs=None,
        rep_max_valid_repetitions=None,
        gross_value_min_value=None,
        gross_value_max_value=None,
        win_var_max_increase_per_sec=None,
        win_var_max_decrease_per_sec=None,
        win_var_time_win_to_check=None,
        win_var_min_num_obs=None,
        step_max_increase_per_sec=None,
        step_max_decrease_per_sec=None,
        buddy_radius=None,
        buddy_min_sample_size=None,
        buddy_max_elev_diff=None,
        buddy_min_std=None,
        buddy_threshold=None,
        buddy_elev_gradient=None,
    ):
        """Update the QC settings for the specified observation type.

        If an argument value is None, the default settings will not be updated.

        Parameters
        ----------
        obstype : str, optional
            The observation type to update the quality control settings for.
            The default is 'temp'.
        dupl_timestamp_keep : bool, optional
            Setting that determines to keep, or remove duplicated timestamps. The default is None.
        persis_time_win_to_check : Timedelta or str, optional
            Time window for persistence check. The default is None.
        persis_min_num_obs : int (> 0), optional
            Minimal window members for persistence check. The default is None.
        rep_max_valid_repetitions : int (> 0), optional
            Maximal valid repetitions for repetitions check. The default is None.
        gross_value_min_value : numeric, optional
            Minimal value for gross value check. The default is None.
        gross_value_max_value : numeric, optional
            Maximal value for gross value check. The default is None.
        win_var_max_increase_per_sec : numeric (> 0), optional
            Maximal increase per second for window variation check. The default is None.
        win_var_max_decrease_per_sec : numeric (> 0), optional
            Maximal decrease per second for window variation check. The default is None.
        win_var_time_win_to_check : Timedelta or str, optional
            Time window for window variation check. The default is None.
        win_var_min_num_obs : int (> 0), optional
            Minimal window members for window variation check. The default is None.
        step_max_increase_per_sec : numeric, optional
            Maximal increase per second for step check. The default is None.
        step_max_decrease_per_sec : numeric (< 0), optional
            Maximal decrease per second for step check. The default is None.
        buddy_radius : numeric (> 0), optional
            The radius to define neighbors in meters. The default is None.
        buddy_min_sample_size : int (> 2), optional
            The minimum sample size to calculate statistics on. The default is
            None.
        buddy_max_elev_diff : numeric (> 0), optional
            The maximum altitude difference allowed for buddies. The default is
            None.
        buddy_min_std : numeric (> 0), optional
            The minimum standard deviation for sample statistics. This should
            represent the accuracy of the observations. The default is None.
        buddy_threshold : numeric (> 0), optional
            The threshold (std units) for flagging observations as buddy
            outliers. The default is None.
        buddy_elev_gradient : numeric, optional
            Describe how the obstype changes with altitude (in meters). The
            default is -0.0065. The default is None.

        Returns
        -------
        None.

        See Also
        -----------
        Dataset.update_titan_qc_settings: Update the QC settings for TITAN checks
        Dataset.apply_quality_control: Apply the default QC pipeline.
        Dataset.apply_buddy_check: Apply spatial buddy check.

        Examples
        --------

        We start by creating a Dataset, and importing data.

        >>> import metobs_toolkit
        >>>
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> dataset.import_data_from_file(
        ...                         input_data_file=metobs_toolkit.demo_datafile,
        ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                         template_file=metobs_toolkit.demo_template,
        ...                         )
        >>> print(dataset)
        Dataset instance containing:
             *28 stations
             *['humidity', 'temp', 'wind_direction', 'wind_speed'] observation types present
             *483828 observation records (not Nan's)
             *0 records labeled as outliers
             *8 gaps
             *records range: 2022-09-01 00:00:00+00:00 --> 2022-09-15 23:55:00+00:00 (total duration:  14 days 23:55:00)
             *time zone of the records: UTC
             *Coordinates are available for all stations.
             *Known GEE datasets for: ['lcz', 'altitude', 'worldcover', 'ERA5-land']



        For this example, we reduce the data by coarsening the time resolution
        to hourly. It is important to resample the time resolution in advance of
        applying quality control since some checks depend on the frequency of the records!

        >>> dataset.coarsen_time_resolution(freq='1h')

        There are default settings for quality control (for temperature). These
        are stored in the `Dataset.settings` attribute. We can inspect them directly,
        or by using the `Datatest.show_settings()` method.

        >>> dataset.settings.qc['qc_check_settings']
        {'duplicated_timestamp': {'keep': False}, 'persistence': {'temp': {'time_window_to_check': '1h', 'min_num_obs': 5}}, 'repetitions': {'temp': {'max_valid_repetitions': 5}}, 'gross_value': {'temp': {'min_value': -15.0, 'max_value': 39.0}}, 'window_variation': {'temp': {'max_increase_per_second': 0.0022222222222222222, 'max_decrease_per_second': 0.002777777777777778, 'time_window_to_check': '1h', 'min_window_members': 3}}, 'step': {'temp': {'max_increase_per_second': 0.0022222222222222222, 'max_decrease_per_second': -0.002777777777777778}}, 'buddy_check': {'temp': {'radius': 15000, 'num_min': 2, 'threshold': 1.5, 'max_elev_diff': 200, 'elev_gradient': -0.0065, 'min_std': 1.0}}}
        >>> dataset.show_settings()
        All settings: ...


        We can change the (default) settings for QC using the `Dataset.update_qc_settings()`
        method. These settings are observationtype dependent!

        >>> dataset.update_qc_settings(
        ...                obstype='temp',
        ...                gross_value_max_value=26.0,
        ...                step_max_increase_per_sec=6.5/3600,
        ...                rep_max_valid_repetitions=4) #depends highly on records frequency!

        >>> dataset.update_qc_settings(
        ...                obstype='humidity',
        ...                gross_value_min_value = 0.,
        ...                gross_value_max_value=100.0)


        Now that the QC settings are adjusted, we can apply standard QC checks.

        >>> dataset.apply_quality_control(obstype='temp')
        >>> dataset.apply_quality_control(obstype='humidity')


        """
        assert (
            obstype in self.obstypes.keys()
        ), f"{obstype} is not a known observation type"

        def _updater(dictionary, obstype, argname, value):
            """Update nested dictionaries."""
            if obstype not in dictionary.keys():
                dictionary[obstype] = {}
                printstr = f"{obstype} : unexisting --> {value}"
            elif argname not in dictionary[obstype]:
                printstr = f"{obstype} : unexisting --> {value}"
            else:
                printstr = f"{obstype} : {dictionary[obstype][argname]} --> {value}"

            dictionary[obstype][argname] = value
            return dictionary, printstr

        # Gross value check
        if gross_value_max_value is not None:
            self.settings.qc["qc_check_settings"]["gross_value"], updatestr = _updater(
                self.settings.qc["qc_check_settings"]["gross_value"],
                obstype=obstype,
                argname="max_value",
                value=float(gross_value_max_value),
            )
            logger.info(f"Maximal value for gross value check updated: {updatestr}")

        if gross_value_min_value is not None:
            self.settings.qc["qc_check_settings"]["gross_value"], updatestr = _updater(
                self.settings.qc["qc_check_settings"]["gross_value"],
                obstype=obstype,
                argname="min_value",
                value=float(gross_value_min_value),
            )
            logger.info(f"Minimal value for gross value check updated: {updatestr}")

        # Duplicate check
        if dupl_timestamp_keep is not None:
            logger.info(
                f'Setting to keep (True) are remove (False) duplicate timestamps updated: \
        {self.settings.qc["qc_check_settings"]["duplicated_timestamp"]["keep"]} -->  {bool(dupl_timestamp_keep)}'
            )
            self.settings.qc["qc_check_settings"]["duplicated_timestamp"]["keep"] = (
                bool(dupl_timestamp_keep)
            )

        # Persistence check
        if persis_time_win_to_check is not None:
            if is_timedelta(str(persis_time_win_to_check)):
                (
                    self.settings.qc["qc_check_settings"]["persistence"],
                    updatestr,
                ) = _updater(
                    self.settings.qc["qc_check_settings"]["persistence"],
                    obstype=obstype,
                    argname="time_window_to_check",
                    value=str(persis_time_win_to_check),
                )

                logger.info(
                    f"Time window size for persistence check updated: {updatestr}"
                )

            else:
                logger.warning(
                    f" {str(persis_time_win_to_check)} is not a valid timedelta string. No update on this setting."
                )

        if persis_min_num_obs is not None:
            self.settings.qc["qc_check_settings"]["persistence"], updatestr = _updater(
                self.settings.qc["qc_check_settings"]["persistence"],
                obstype=obstype,
                argname="min_num_obs",
                value=abs(int(persis_min_num_obs)),
            )

            logger.info(
                f"Minimal window members for persistence check updated: {updatestr}"
            )

        # Repetitions check
        if rep_max_valid_repetitions is not None:
            self.settings.qc["qc_check_settings"]["repetitions"], updatestr = _updater(
                self.settings.qc["qc_check_settings"]["repetitions"],
                obstype=obstype,
                argname="max_valid_repetitions",
                value=abs(int(rep_max_valid_repetitions)),
            )
            logger.info(
                f"Maximal valid repetitions for repetitions check updated: {updatestr}"
            )

        # Window variation check
        if win_var_max_increase_per_sec is not None:
            (
                self.settings.qc["qc_check_settings"]["window_variation"],
                updatestr,
            ) = _updater(
                self.settings.qc["qc_check_settings"]["window_variation"],
                obstype=obstype,
                argname="max_increase_per_second",
                value=abs(float(win_var_max_increase_per_sec)),
            )

            logger.info(
                f"Maximal increase per second for window variation check updated: {updatestr}"
            )

        if win_var_max_decrease_per_sec is not None:
            (
                self.settings.qc["qc_check_settings"]["window_variation"],
                updatestr,
            ) = _updater(
                self.settings.qc["qc_check_settings"]["window_variation"],
                obstype=obstype,
                argname="max_decrease_per_second",
                value=abs(float(win_var_max_decrease_per_sec)),
            )
            logger.info(
                f"Maximal decrease per second for window variation check updated: {updatestr}"
            )

        if win_var_time_win_to_check is not None:
            if is_timedelta(str(win_var_time_win_to_check)):
                (
                    self.settings.qc["qc_check_settings"]["window_variation"],
                    updatestr,
                ) = _updater(
                    self.settings.qc["qc_check_settings"]["window_variation"],
                    obstype=obstype,
                    argname="time_window_to_check",
                    value=str(win_var_time_win_to_check),
                )
                logger.info(
                    f"Time window for window variation check updated: {updatestr}"
                )
            else:
                logger.warning(
                    f" {str(persis_time_win_to_check)} is not a valid timedelta string. No update on this setting."
                )

        if win_var_min_num_obs is not None:
            (
                self.settings.qc["qc_check_settings"]["window_variation"],
                updatestr,
            ) = _updater(
                self.settings.qc["qc_check_settings"]["window_variation"],
                obstype=obstype,
                argname="min_window_members",
                value=abs(int(win_var_min_num_obs)),
            )
            logger.info(
                f"Minimal window members for window variation check updated: {updatestr}"
            )

        # Step check
        if step_max_increase_per_sec is not None:
            self.settings.qc["qc_check_settings"]["step"], updatestr = _updater(
                self.settings.qc["qc_check_settings"]["step"],
                obstype=obstype,
                argname="max_increase_per_second",
                value=abs(float(step_max_increase_per_sec)),
            )

            logger.info(
                f"Maximal increase per second for step check updated: {updatestr}"
            )

        if step_max_decrease_per_sec is not None:
            self.settings.qc["qc_check_settings"]["step"], updatestr = _updater(
                self.settings.qc["qc_check_settings"]["step"],
                obstype=obstype,
                argname="max_decrease_per_second",
                value=-1.0 * abs(float(step_max_decrease_per_sec)),
            )

            logger.info(
                f"Maximal decrease per second for step check updated: {updatestr}"
            )

        # Buddy check
        buddy_elev_gradient = None
        if buddy_radius is not None:
            self.settings.qc["qc_check_settings"]["buddy_check"], updatestr = _updater(
                self.settings.qc["qc_check_settings"]["buddy_check"],
                obstype=obstype,
                argname="radius",
                value=abs(float(buddy_radius)),
            )
            logger.info(f"Buddy radius for buddy check updated: {updatestr}")

        if buddy_min_sample_size is not None:
            value = abs(int(buddy_min_sample_size))
            if value >= 2:
                (
                    self.settings.qc["qc_check_settings"]["buddy_check"],
                    updatestr,
                ) = _updater(
                    self.settings.qc["qc_check_settings"]["buddy_check"],
                    obstype=obstype,
                    argname="num_min",
                    value=value,
                )
                logger.info(
                    f"Minimum number of buddies for buddy check updated: {updatestr}"
                )
            else:
                logger.warning(
                    f"Minimum number of buddies must be >= 2, but {value} is given. Not updated."
                )

        if buddy_max_elev_diff is not None:
            self.settings.qc["qc_check_settings"]["buddy_check"], updatestr = _updater(
                self.settings.qc["qc_check_settings"]["buddy_check"],
                obstype=obstype,
                argname="max_elev_diff",
                value=abs(float(buddy_max_elev_diff)),
            )
            logger.info(
                f"Max elevation differences for buddy check updated: {updatestr}"
            )

        if buddy_min_std is not None:
            self.settings.qc["qc_check_settings"]["buddy_check"], updatestr = _updater(
                self.settings.qc["qc_check_settings"]["buddy_check"],
                obstype=obstype,
                argname="min_std",
                value=abs(float(buddy_min_std)),
            )
            logger.info(f"Minimum std in sample for buddy check updated: {updatestr}")

        if buddy_threshold is not None:
            self.settings.qc["qc_check_settings"]["buddy_check"], updatestr = _updater(
                self.settings.qc["qc_check_settings"]["buddy_check"],
                obstype=obstype,
                argname="threshold",
                value=abs(float(buddy_threshold)),
            )
            logger.info(
                f"Outlier threshold (in sigma) for buddy check updated: {updatestr}"
            )

        if buddy_elev_gradient is not None:
            self.settings.qc["qc_check_settings"]["buddy_check"], updatestr = _updater(
                self.settings.qc["qc_check_settings"]["buddy_check"],
                obstype=obstype,
                argname="elev_gradient",
                value=float(buddy_max_elev_diff),
            )
            logger.info(f"Elevation gradient for buddy check updated: {updatestr}")

    # def update_titan_qc_settings(
    #     self,
    #     obstype="temp",
    #     # buddy settings
    #     buddy_radius=None,
    #     buddy_num_min=None,
    #     buddy_threshold=None,
    #     buddy_max_elev_diff=None,
    #     buddy_elev_gradient=None,
    #     buddy_min_std=None,
    #     buddy_num_iterations=None,
    #     buddy_debug=None,
    #     # sct settings
    #     sct_num_min_outer=None,
    #     sct_num_max_outer=None,
    #     sct_inner_radius=None,
    #     sct_outer_radius=None,
    #     sct_num_iterations=None,
    #     sct_num_min_prof=None,
    #     sct_min_elev_diff=None,
    #     sct_min_horizontal_scale=None,
    #     sct_max_horizontal_scale=None,
    #     sct_kth_closest_obs_horizontal_scale=None,
    #     sct_vertical_scale=None,
    #     sct_mina_deviation=None,  # vec Minimum admissible value
    #     sct_maxa_deviation=None,  # vec Maximum admissible value
    #     sct_minv_deviation=None,  # vec Minimum valid value
    #     sct_maxv_deviation=None,  # vec Maximum valid value
    #     sct_eps2=None,  # Ratio of observation error variance to background variance
    #     sct_tpos=None,  # vec Positive deviation allowed
    #     sct_tneg=None,  # vec Negative deviation allowed
    #     sct_basic=None,
    #     sct_debug=None,
    # ):
    # """Update the TITAN QC settings for the specified observation type.

    # If a argument value is None, the default settings will not be updated.

    # For a detailed explanation of the settings, we refer to the
    # [TITAN documetation](https://github.com/metno/titanlib/wiki)

    # Parameters
    # ----------
    # The observation type to update the quality control settings for.
    # The default is 'temp'.
    # buddy_radius : int (> 0), optional
    #     Search radius in m. The default is None.
    # buddy_num_min : int (> 0), optional
    #     The minimum number of buddies a station can have. The default is
    #     None.
    # buddy_threshold : num (> 0), optional
    #     The variance threshold for flagging a station. The default is None.
    # buddy_max_elev_diff : num, optional
    #     The maximum difference in elevation for a buddy (if negative will not check for heigh difference). The default is None.
    # buddy_elev_gradient : num, optional
    #     Linear elevation temperature gradient with height. The default is None.
    # buddy_min_std : num (> 0), optional
    #     If the standard deviation of values in a neighborhood are less than min_std, min_std will be used instead. The default is None.
    # buddy_num_iterations : int (> 0), optional
    #     The number of iterations to perform. The default is None.
    # buddy_debug : bool, optional
    #     If True, print out debug information. The default is None.
    # sct_num_min_outer : int (> 0), optional
    #     Minimal points in outer circle. The default is None.
    # sct_num_max_outer : int (> 0), optional
    #     Maximal points in outer circle. The default is None.
    # sct_inner_radius : num (> 0), optional
    #     Radius of inner circle. The default is None.
    # sct_outer_radius : num (> 0), optional
    #     Radius of outer circle. The default is None.
    # sct_num_iterations : int (> 0), optional
    #     Number of iterations. The default is None.
    # sct_num_min_prof : int (> 0), optional
    #     Minimum number of observations to compute vertical profile. The default is None.
    # sct_min_elev_diff : num (> 0), optional
    #     Minimum elevation difference to compute vertical profile. The default is None.
    # sct_min_horizontal_scale : num (> 0), optional
    #     Minimum horizontal decorrelation length. The default is None.
    # sct_max_horizontal_scale : num (> 0), optional
    #     Maximum horizontal decorrelation length. The default is None.
    # sct_kth_closest_obs_horizontal_scale : int (> 0), optional
    #     Number of closest observations to consider. The default is None.
    # sct_vertical_scale : num (> 0), optional
    #     Vertical decorrelation length. The default is None.
    # sct_mina_deviation : num (> 0), optional
    #     Minimum admissible value deviation. The default is None.
    # sct_maxa_deviation : num (> 0), optional
    #     Maximum admissible value deviation. The default is None.
    # sct_minv_deviation : num (> 0), optional
    #     Minimum valid value deviation. The default is None.
    # sct_maxv_deviation : num (> 0), optional
    #     Maximum valid value deviation. The default is None.
    # sct_eps2 : num (> 0), optional
    #     Ratio of observation error variance to background variance. The default is None.
    # sct_tpos : num (> 0), optional
    #     Positive deviation allowed. The default is None.
    # sct_tneg : num (> 0), optional
    #     Positive deviation allowed. The default is None.
    # sct_basic : bool, optional
    #     Basic mode. The default is None.
    # sct_debug : bool, optional
    #     If True, print out debug information. The default is None.

    # Returns
    # -------
    # None.

    # """
    # assert (
    #     obstype in self.obstypes.keys()
    # ), f"{obstype} is not a known observation type"

    # # check buddy settings for updates
    # buddy_attrs = {
    #     "buddy_radius": {"new_value": buddy_radius, "dtype": "numeric"},
    #     "buddy_num_min": {"new_value": buddy_num_min, "dtype": "int"},
    #     "buddy_threshold": {"new_value": buddy_threshold, "dtype": "numeric"},
    #     "buddy_max_elev_diff": {
    #         "new_value": buddy_max_elev_diff,
    #         "dtype": "numeric",
    #     },
    #     "buddy_elev_gradient": {
    #         "new_value": buddy_elev_gradient,
    #         "dtype": "numeric",
    #     },
    #     "buddy_min_std": {"new_value": buddy_min_std, "dtype": "numeric"},
    #     "buddy_num_iterations": {"new_value": buddy_num_iterations, "dtype": "int"},
    #     "buddy_debug": {"new_value": buddy_debug, "dtype": "bool"},
    # }

    # sct_attrs = {
    #     "sct_num_min_outer": {"new_value": sct_num_min_outer, "dtype": "int"},
    #     "sct_num_max_outer": {"new_value": sct_num_max_outer, "dtype": "int"},
    #     "sct_inner_radius": {"new_value": sct_inner_radius, "dtype": "numeric"},
    #     "sct_outer_radius": {"new_value": sct_outer_radius, "dtype": "numeric"},
    #     "sct_num_iterations": {"new_value": sct_num_iterations, "dtype": "int"},
    #     "sct_num_min_prof": {"new_value": sct_num_min_prof, "dtype": "int"},
    #     "sct_min_elev_diff": {"new_value": sct_min_elev_diff, "dtype": "numeric"},
    #     "sct_min_horizontal_scale": {
    #         "new_value": sct_min_horizontal_scale,
    #         "dtype": "numeric",
    #     },
    #     "sct_max_horizontal_scale": {
    #         "new_value": sct_max_horizontal_scale,
    #         "dtype": "numeric",
    #     },
    #     "sct_kth_closest_obs_horizontal_scale": {
    #         "new_value": sct_kth_closest_obs_horizontal_scale,
    #         "dtype": "int",
    #     },
    #     "sct_vertical_scale": {"new_value": sct_vertical_scale, "dtype": "numeric"},
    #     "sct_mina_deviation": {"new_value": sct_mina_deviation, "dtype": "numeric"},
    #     "sct_minv_deviation": {"new_value": sct_minv_deviation, "dtype": "numeric"},
    #     "sct_maxv_deviation": {"new_value": sct_maxv_deviation, "dtype": "numeric"},
    #     "sct_eps2": {"new_value": sct_eps2, "dtype": "numeric"},
    #     "sct_tpos": {"new_value": sct_tpos, "dtype": "numeric"},
    #     "sct_tneg": {"new_value": sct_tneg, "dtype": "numeric"},
    #     "sct_basic": {"new_value": sct_basic, "dtype": "bool"},
    #     "sct_debug": {"new_value": sct_debug, "dtype": "bool"},
    # }
    def update_titan_qc_settings(
        self,
        obstype="temp",
        # buddy settings
        buddy_radius=None,
        buddy_num_min=None,
        buddy_threshold=None,
        buddy_max_elev_diff=None,
        buddy_elev_gradient=None,
        buddy_min_std=None,
        buddy_num_iterations=None,
        buddy_debug=None,
    ):
        """Update the TITAN QC settings for the specified observation type.

        If an argument value is None, the default settings will not be updated.

        For a detailed explanation of the settings, we refer to the
        [TITAN documetation](https://github.com/metno/titanlib/wiki)

        Parameters
        ----------
        obstype : str, optional
        The observation type to update the quality control settings for.
        The default is 'temp'.
        buddy_radius : int (> 0), optional
            Search radius in m. The default is None.
        buddy_num_min : int (> 0), optional
            The minimum number of buddies a station can have. The default is
            None.
        buddy_threshold : num (> 0), optional
            The variance threshold for flagging a station. The default is None.
        buddy_max_elev_diff : num, optional
            The maximum difference in elevation for a buddy (if negative will
            not check for height difference). The default is None.
        buddy_elev_gradient : num, optional
            Linear elevation temperature gradient with height. The default is None.
        buddy_min_std : num (> 0), optional
            If the standard deviation of values in a neighborhood is less than min_std, min_std will be used instead. The default is None.
        buddy_num_iterations : int (> 0), optional
            The number of iterations to perform. The default is None.
        buddy_debug : bool, optional
            If True, print out debug information. The default is None.

        Returns
        -------
        None.

        See Also
        -----------
        Dataset.update_qc_settings: Update the settings for the default QC pipeline.
        Dataset.apply_titan_buddy_check: Apply TITAN buddy check

        Examples
        --------

        We start by creating a Dataset, and importing data.

        >>> import metobs_toolkit
        >>>
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> dataset.import_data_from_file(
        ...                         input_data_file=metobs_toolkit.demo_datafile,
        ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                         template_file=metobs_toolkit.demo_template,
        ...                         )
        >>> print(dataset)
        Dataset instance containing:
             *28 stations
             *['humidity', 'temp', 'wind_direction', 'wind_speed'] observation types present
             *483828 observation records (not Nan's)
             *0 records labeled as outliers
             *8 gaps
             *records range: 2022-09-01 00:00:00+00:00 --> 2022-09-15 23:55:00+00:00 (total duration:  14 days 23:55:00)
             *time zone of the records: UTC
             *Coordinates are available for all stations.
             *Known GEE datasets for: ['lcz', 'altitude', 'worldcover', 'ERA5-land']

        For this example, we reduce the data by coarsening the time resolution
        to hourly.

        >>> dataset.coarsen_time_resolution(freq='1h')

        There are default TITAN buddy check settings (for temperature). These
        are stored in the `Dataset.settings` attribute. We can inspect them directly,
        or by using the `Datatest.show_settings()` method.

        >>> dataset.settings.qc['titan_check_settings']['titan_buddy_check']
        {'temp': {'radius': 50000, 'num_min': 2, 'threshold': 1.5, 'max_elev_diff': 200, 'elev_gradient': -0.0065, 'min_std': 1.0, 'num_iterations': 1}}

        We can change the (default) settings using the `Dataset.update_titan_qc_settings()`
        method. These settings are observationtype dependent!

        >>> # The following settings are illustrative, do not copy blindly
        >>> dataset.update_titan_qc_settings(
        ...                obstype='temp',
        ...                buddy_radius=20000,
        ...                buddy_num_min=4,
        ...                buddy_threshold=2.5,
        ...                buddy_min_std=1.)
        buddy radius for the TITAN buddy check updated:  50000--> 20000.0
        buddy num min for the TITAN buddy check updated:  2--> 4
        buddy threshold for the TITAN buddy check updated:  1.5--> 2.5
        buddy min std for the TITAN buddy check updated:  1.0--> 1.0

        """
        assert (
            obstype in self.obstypes.keys()
        ), f"{obstype} is not a known observation type"

        # check buddy settings for updates
        buddy_attrs = {
            "buddy_radius": {"new_value": buddy_radius, "dtype": "numeric"},
            "buddy_num_min": {"new_value": buddy_num_min, "dtype": "int"},
            "buddy_threshold": {"new_value": buddy_threshold, "dtype": "numeric"},
            "buddy_max_elev_diff": {
                "new_value": buddy_max_elev_diff,
                "dtype": "numeric",
            },
            "buddy_elev_gradient": {
                "new_value": buddy_elev_gradient,
                "dtype": "numeric",
            },
            "buddy_min_std": {"new_value": buddy_min_std, "dtype": "numeric"},
            "buddy_num_iterations": {"new_value": buddy_num_iterations, "dtype": "int"},
            "buddy_debug": {"new_value": buddy_debug, "dtype": "bool"},
        }

        def _iterate_attributes(obstype, attr_dict, attr_prefix, checkname):

            if obstype not in self.settings.qc["titan_check_settings"][checkname]:
                self.settings.qc["titan_check_settings"][checkname][obstype] = {}

            for key, val in attr_dict.items():
                if not val["new_value"] is None:
                    settings_key = key.split(attr_prefix)[1]  # remove 'buddy_'
                    if val["dtype"] == "numeric":
                        new_val = float(val["new_value"])
                    elif val["dtype"] == "int":
                        new_val = int(val["new_value"])
                    elif val["dtype"] == "bool":
                        new_val = bool(val["new_value"])
                    else:  # val['dtype'] == 'str':
                        new_val = str(val["new_value"])

                    try:
                        old_value = self.settings.qc["titan_check_settings"][checkname][
                            obstype
                        ][settings_key]
                        print(
                            f'{key.replace("_", " ")} for the TITAN buddy check updated:  {old_value}--> {new_val}'
                        )
                    except KeyError:
                        print(
                            f'{key.replace("_", " ")} for the TITAN buddy check added:  --> {new_val}'
                        )

                    self.settings.qc["titan_check_settings"][checkname][obstype][
                        settings_key
                    ] = new_val

        _iterate_attributes(obstype, buddy_attrs, "buddy_", "titan_buddy_check")
        # _iterate_attributes(obstype, sct_attrs, "sct_", "titan_sct_resistant_check")


# =============================================================================
# dtype check functions
# =============================================================================


def is_timedelta(timedeltastr):
    """Test if string can be timedelta representation.

    Parameters
    ----------
    timedeltastr : str
        Representation of timedelta.

    Returns
    -------
    bool


    """
    try:
        pd.to_timedelta(timedeltastr)
        return True
    except ValueError:
        return False


# =============================================================================
# Exceptions
# =============================================================================


class MetobsDatasetSettingsUpdaterError(Exception):
    """Exception raised for errors when updating settings."""

    pass


# =============================================================================
# Docstring test
# =============================================================================
if __name__ == "__main__":
    from metobs_toolkit.doctest_fmt import setup_and_run_doctest

    setup_and_run_doctest()
