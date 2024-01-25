#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Extension of the Dataset class (methods for updating settings).
@author: thoverga
"""
import logging
import pandas as pd


import metobs_toolkit.dataset as dataset

logger = logging.getLogger(__name__)


class Dataset(dataset.Dataset):
    """Extension on the metobs_toolkit.Dataset class with updaters."""

    def update_settings(
        self,
        output_folder=None,
        input_data_file=None,
        input_metadata_file=None,
        template_file=None,
    ):
        """Update the most common input-output (IO) settings.

        (This should be applied before importing the observations.)

        When an update value is None, the specific setting will not be updated.

        Parameters
        ----------
        output_folder : string, optional
            A directory to store the output to. The default is None.
        input_data_file : string, optional
            Path to the input data file with observations. The default is None.
        input_metadata_file : string, optional
            Path to the input metadata file. The default is None.
        template_file : string, optional
            Path to the mapper-template csv file to be used on the observations
            and metadata. The default is None.

        Returns
        -------
        None.

        """
        self.settings.update_IO(
            output_folder=output_folder,
            input_data_file=input_data_file,
            input_metadata_file=input_metadata_file,
            template_file=template_file,
        )

    def update_timezone(self, timezonestr):
        """Change the timezone of the input data.

        By default UTC is assumed.
        A valid timezonestring is an element of the pytz.all_timezones.

        Parameters
        ----------
        timezonestr : string
            Timezone string of the input observations. Element of pytz.all_timezones.

        Returns
        -------
        None.

        """
        self.settings.update_timezone(timezonestr)

    def update_default_name(self, default_name):
        """Update the default name (the name of the station).

        This name will be used when no names are found in the observational dataset.

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

    def update_gap_and_missing_fill_settings(
        self,
        gap_interpolation_method=None,
        gap_interpolation_max_consec_fill=None,
        gap_debias_prefered_leading_period_hours=None,
        gap_debias_prefered_trailing_period_hours=None,
        gap_debias_minimum_leading_period_hours=None,
        gap_debias_minimum_trailing_period_hours=None,
        automatic_max_interpolation_duration_str=None,
        missing_obs_interpolation_method=None,
    ):
        """Update fill settings for gaps and missing observations.

        If None, the current setting is not updated.

        Parameters
        ----------
        gap_interpolation_method : str, optional
            The interpolation method to pass to numpy.interpolate. The default is None.
        gap_interpolation_max_consec_fill : int, optional
            Maximum number of lacking observations to interpolate. This is
            passed to the limit argument of Numpy.interpolate. The default is
            None.
        gap_debias_prefered_leading_period_hours : int, optional
            The preferd size of the leading period for calculating hourly
            biasses wrt the model. The default is None.
        gap_debias_prefered_trailing_period_hours : int, optional
            The preferd size of the trailing period for calculating hourly
            biasses wrt the model. The default is None.
        gap_debias_minimum_leading_period_hours : int, optional
            The minimum size of the leading period for calculating hourly
            biasses wrt the model. The default is None.
        gap_debias_minimum_trailing_period_hours : int, optional
            The minimum size of the trailing period for calculating hourly
            biasses wrt the model. The default is None.
        automatic_max_interpolation_duration_str : Timedelta or str, optional
            Maximum duration to apply interpolation for gapfill when using the
            automatic gapfill method. Gaps with longer durations will be filled
            using debiased modeldata. The default is None.
        missing_obs_interpolation_method : str, optional
            The interpolation method to pass to numpy.interpolate. The default is None.

        Returns
        -------
        None.

        """
        # Gap linear interpolation
        if gap_interpolation_method is not None:
            logger.info(
                f' The gap interpolation method is updated: \
        {self.settings.gap["gaps_fill_settings"]["linear"]["method"]} --> {str(gap_interpolation_method)}'
            )
            self.settings.gap["gaps_fill_settings"]["linear"]["method"] = str(
                gap_interpolation_method
            )

        if gap_interpolation_max_consec_fill is not None:
            logger.info(
                f' The gap max number of consecutive interpolations is updated: \
        {self.settings.gap["gaps_fill_settings"]["linear"]["max_consec_fill"]} --> {abs(int(gap_interpolation_max_consec_fill))}'
            )
            self.settings.gap["gaps_fill_settings"]["linear"]["max_consec_fill"] = abs(
                int(gap_interpolation_max_consec_fill)
            )

        # Gap debias fill
        if gap_debias_prefered_leading_period_hours is not None:
            logger.info(
                f' The size of the prefered leading period for debias gapfill is updated: \
        {self.settings.gap["gaps_fill_settings"]["model_debias"]["debias_period"]["prefered_leading_sample_duration_hours"]} --> {abs(int(gap_debias_prefered_leading_period_hours))}'
            )
            self.settings.gap["gaps_fill_settings"]["model_debias"]["debias_period"][
                "prefered_leading_sample_duration_hours"
            ] = abs(int(gap_debias_prefered_leading_period_hours))

        if gap_debias_prefered_trailing_period_hours is not None:
            logger.info(
                f' The size of the prefered trailing period for debias gapfill is updated: \
        {self.settings.gap["gaps_fill_settings"]["model_debias"]["debias_period"]["prefered_trailing_sample_duration_hours"]} --> {abs(int(gap_debias_prefered_trailing_period_hours))}'
            )
            self.settings.gap["gaps_fill_settings"]["model_debias"]["debias_period"][
                "prefered_trailing_sample_duration_hours"
            ] = abs(int(gap_debias_prefered_trailing_period_hours))

        if gap_debias_minimum_leading_period_hours is not None:
            logger.info(
                f' The minimum size of the leading period for debias gapfill is updated: \
        {self.settings.gap["gaps_fill_settings"]["model_debias"]["debias_period"]["minimum_leading_sample_duration_hours"]} --> {abs(int(gap_debias_minimum_leading_period_hours))}'
            )
            self.settings.gap["gaps_fill_settings"]["model_debias"]["debias_period"][
                "minimum_leading_sample_duration_hours"
            ] = abs(int(gap_debias_minimum_leading_period_hours))

        if gap_debias_minimum_trailing_period_hours is not None:
            logger.info(
                f' The minimum size of the trailing period for debias gapfill is updated: \
        {self.settings.gap["gaps_fill_settings"]["model_debias"]["debias_period"]["minimum_trailing_sample_duration_hours"]} --> {abs(int(gap_debias_minimum_trailing_period_hours))}'
            )
            self.settings.gap["gaps_fill_settings"]["model_debias"]["debias_period"][
                "minimum_trailing_sample_duration_hours"
            ] = abs(int(gap_debias_minimum_trailing_period_hours))

        # Gapfill automatic
        if automatic_max_interpolation_duration_str is not None:
            if is_timedelta(str(automatic_max_interpolation_duration_str)):
                logger.info(
                    f' The maximum interpolation duration for automatic gapfill is updated: \
            {self.settings.gap["gaps_fill_settings"]["automatic"]["max_interpolation_duration_str"]} --> {str(automatic_max_interpolation_duration_str)}'
                )
                self.settings.gap["gaps_fill_settings"]["automatic"][
                    "max_interpolation_duration_str"
                ] = str(automatic_max_interpolation_duration_str)
            else:
                logger.warning(
                    f" {str(automatic_max_interpolation_duration_str)} is not a valid timedelta string. No update on this setting."
                )

        # Missing obs interpolation
        if missing_obs_interpolation_method is not None:
            logger.info(
                f' The missing observations interpolation method is updated: \
        {self.settings.missing_obs["missing_obs_fill_settings"]["linear"]["method"]} --> {str(missing_obs_interpolation_method)}'
            )
            self.settings.missing_obs["missing_obs_fill_settings"]["linear"][
                "method"
            ] = str(missing_obs_interpolation_method)

    def update_qc_settings(
        self,
        obstype="temp",
        gapsize_in_records=None,
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

        If a argument value is None, the default settings will not be updated.

        Parameters
        ----------
        obstype : str, optional
            The observation type to update the quality control settings for.
            The default is 'temp'.
        gapsize_in_records : int (> 0), optional
            A gap is defined as a sequence of missing observations with a length
            greater or equal to this number, on the input frequencies. The default is None.
        dupl_timestamp_keep : bool, optional
            Setting that determines to keep, or remove duplicated timestamps. The default is None.
        persis_time_win_to_check :automatic_max_interpolation_duration_str
            Time window for persistance check. The default is None.
        persis_min_num_obs : int (> 0), optional
            Minimal window members for persistance check. The default is None.
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
            The radius to define neighbours in meters. The default is None.
        buddy_min_sample_size : int (> 2), optional
            The minimum sample size to calculate statistics on. The default is
            None.
        buddy_max_elev_diff : numeric (> 0), optional
            The maximum altitude difference allowed for buddies. The default is
            None.
        buddy_min_std : numeric (> 0), optional
            The minimum standard deviation for sample statistics. This should
            represent the accuracty of the observations. The default is None.
        buddy_threshold : numeric (> 0), optional
            The threshold (std units) for flaggging observations as buddy
            outliers. The default is None.
        buddy_elev_gradient : numeric, optional
            Describes how the obstype changes with altitude (in meters). The
            default is -0.0065. The default is None.

        Returns
        -------
        None.

        Note
        -------
        The gap defenition is independend of the observation type, and is thus set for
        all the observation types.

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

        # Gap defenition
        if gapsize_in_records is not None:
            logger.info(
                f' The defenition of a gap (=gapsize) is updated: \
        {self.settings.gap["gaps_settings"]["gaps_finder"]["gapsize_n"]} --> {abs(int(gapsize_in_records))}'
            )
            self.settings.gap["gaps_settings"]["gaps_finder"]["gapsize_n"] = abs(
                int(gapsize_in_records)
            )

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
            self.settings.qc["qc_check_settings"]["duplicated_timestamp"][
                "keep"
            ] = bool(dupl_timestamp_keep)

        # Persistance check
        if persis_time_win_to_check is not None:
            if is_timedelta(str(persis_time_win_to_check)):
                (
                    self.settings.qc["qc_check_settings"]["persistance"],
                    updatestr,
                ) = _updater(
                    self.settings.qc["qc_check_settings"]["persistance"],
                    obstype=obstype,
                    argname="time_window_to_check",
                    value=str(persis_time_win_to_check),
                )

                logger.info(
                    f"Time window size for persistance check updated: {updatestr}"
                )

            else:
                logger.warning(
                    f" {str(persis_time_win_to_check)} is not a valid timedelta string. No update on this setting."
                )

        if persis_min_num_obs is not None:
            self.settings.qc["qc_check_settings"]["persistance"], updatestr = _updater(
                self.settings.qc["qc_check_settings"]["persistance"],
                obstype=obstype,
                argname="min_num_obs",
                value=abs(int(persis_min_num_obs)),
            )

            logger.info(
                f"Minimal window members for persistance check updated: {updatestr}"
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
        # sct settings
        sct_num_min_outer=None,
        sct_num_max_outer=None,
        sct_inner_radius=None,
        sct_outer_radius=None,
        sct_num_iterations=None,
        sct_num_min_prof=None,
        sct_min_elev_diff=None,
        sct_min_horizontal_scale=None,
        sct_max_horizontal_scale=None,
        sct_kth_closest_obs_horizontal_scale=None,
        sct_vertical_scale=None,
        sct_mina_deviation=None,  # vec Minimum admissible value
        sct_maxa_deviation=None,  # vec Maximum admissible value
        sct_minv_deviation=None,  # vec Minimum valid value
        sct_maxv_deviation=None,  # vec Maximum valid value
        sct_eps2=None,  # Ratio of observation error variance to background variance
        sct_tpos=None,  # vec Positive deviation allowed
        sct_tneg=None,  # vec Negative deviation allowed
        sct_basic=None,
        sct_debug=None,
    ):
        """Update the TITAN QC settings for the specified observation type.

        If a argument value is None, the default settings will not be updated.

        For a detailed explanation of the settings, we refer to the
        [TITAN documetation](https://github.com/metno/titanlib/wiki)

        Parameters
        ----------
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
            The maximum difference in elevation for a buddy (if negative will not check for heigh difference). The default is None.
        buddy_elev_gradient : num, optional
            Linear elevation temperature gradient with height. The default is None.
        buddy_min_std : num (> 0), optional
            If the standard deviation of values in a neighborhood are less than min_std, min_std will be used instead. The default is None.
        buddy_num_iterations : int (> 0), optional
            The number of iterations to perform. The default is None.
        buddy_debug : bool, optional
            If True, print out debug information. The default is None.
        sct_num_min_outer : int (> 0), optional
            Minimal points in outer circle. The default is None.
        sct_num_max_outer : int (> 0), optional
            Maximal points in outer circle. The default is None.
        sct_inner_radius : num (> 0), optional
            Radius of inner circle. The default is None.
        sct_outer_radius : num (> 0), optional
            Radius of outer circle. The default is None.
        sct_num_iterations : int (> 0), optional
            Number of iterations. The default is None.
        sct_num_min_prof : int (> 0), optional
            Minimum number of observations to compute vertical profile. The default is None.
        sct_min_elev_diff : num (> 0), optional
            Minimum elevation difference to compute vertical profile. The default is None.
        sct_min_horizontal_scale : num (> 0), optional
            Minimum horizontal decorrelation length. The default is None.
        sct_max_horizontal_scale : num (> 0), optional
            Maximum horizontal decorrelation length. The default is None.
        sct_kth_closest_obs_horizontal_scale : int (> 0), optional
            Number of closest observations to consider. The default is None.
        sct_vertical_scale : num (> 0), optional
            Vertical decorrelation length. The default is None.
        sct_mina_deviation : num (> 0), optional
            Minimum admissible value deviation. The default is None.
        sct_maxa_deviation : num (> 0), optional
            Maximum admissible value deviation. The default is None.
        sct_minv_deviation : num (> 0), optional
            Minimum valid value deviation. The default is None.
        sct_maxv_deviation : num (> 0), optional
            Maximum valid value deviation. The default is None.
        sct_eps2 : num (> 0), optional
            Ratio of observation error variance to background variance. The default is None.
        sct_tpos : num (> 0), optional
            Positive deviation allowed. The default is None.
        sct_tneg : num (> 0), optional
            Positive deviation allowed. The default is None.
        sct_basic : bool, optional
            Basic mode. The default is None.
        sct_debug : bool, optional
            If True, print out debug information. The default is None.

        Returns
        -------
        None.

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

        sct_attrs = {
            "sct_num_min_outer": {"new_value": sct_num_min_outer, "dtype": "int"},
            "sct_num_max_outer": {"new_value": sct_num_max_outer, "dtype": "int"},
            "sct_inner_radius": {"new_value": sct_inner_radius, "dtype": "numeric"},
            "sct_outer_radius": {"new_value": sct_outer_radius, "dtype": "numeric"},
            "sct_num_iterations": {"new_value": sct_num_iterations, "dtype": "int"},
            "sct_num_min_prof": {"new_value": sct_num_min_prof, "dtype": "int"},
            "sct_min_elev_diff": {"new_value": sct_min_elev_diff, "dtype": "numeric"},
            "sct_min_horizontal_scale": {
                "new_value": sct_min_horizontal_scale,
                "dtype": "numeric",
            },
            "sct_max_horizontal_scale": {
                "new_value": sct_max_horizontal_scale,
                "dtype": "numeric",
            },
            "sct_kth_closest_obs_horizontal_scale": {
                "new_value": sct_kth_closest_obs_horizontal_scale,
                "dtype": "int",
            },
            "sct_vertical_scale": {"new_value": sct_vertical_scale, "dtype": "numeric"},
            "sct_mina_deviation": {"new_value": sct_mina_deviation, "dtype": "numeric"},
            "sct_minv_deviation": {"new_value": sct_minv_deviation, "dtype": "numeric"},
            "sct_maxv_deviation": {"new_value": sct_maxv_deviation, "dtype": "numeric"},
            "sct_eps2": {"new_value": sct_eps2, "dtype": "numeric"},
            "sct_tpos": {"new_value": sct_tpos, "dtype": "numeric"},
            "sct_tneg": {"new_value": sct_tneg, "dtype": "numeric"},
            "sct_basic": {"new_value": sct_basic, "dtype": "bool"},
            "sct_debug": {"new_value": sct_debug, "dtype": "bool"},
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
        _iterate_attributes(obstype, sct_attrs, "sct_", "titan_sct_resistant_check")


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
