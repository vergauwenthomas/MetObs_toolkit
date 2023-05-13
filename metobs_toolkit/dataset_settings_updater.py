
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Seperate interface for the methods that update dataset settings
@author: thoverga
"""
import logging
import pandas as pd


import metobs_toolkit.dataset as dataset
from metobs_toolkit import observation_types
logger = logging.getLogger(__name__)

class Dataset(dataset.Dataset):

    def update_settings(
        self,
        output_folder=None,
        input_data_file=None,
        input_metadata_file=None,
        data_template_file=None,
        metadata_template_file=None,
    ):
        """
        Update the most common input-output (IO) settings.
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
        data_template_file : string, optional
            Path to the mapper-template csv file to be used on the observations.. The default is None.
        metadata_template_file : string, optional
            Path to the mapper-template csv file to be used on the metadata.. The default is None.

        Returns
        -------
        None.

        """

        self.settings.update_IO(
            output_folder=output_folder,
            input_data_file=input_data_file,
            input_metadata_file=input_metadata_file,
            data_template_file=data_template_file,
            metadata_template_file=metadata_template_file,
        )


    def update_timezone(self, timezonestr):
        """
        Change the timezone of the input data. By default the Brussels timezone is assumed.
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
        """
        Update the default name (the name of the station). This name will be
        used when no names are found in the observational dataset.

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



    def update_gap_and_missing_fill_settings(self, gap_interpolation_method=None, gap_interpolation_max_consec_fill=None,
                                             gap_debias_prefered_leading_period_hours=None, gap_debias_prefered_trailing_period_hours=None,
                                             gap_debias_minimum_leading_period_hours=None, gap_debias_minimum_trailing_period_hours=None,
                                             automatic_max_interpolation_duration_str=None, missing_obs_interpolation_method=None):
        """
        Update the settings on the filling methods for missing observations and gaps.

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
        if not gap_interpolation_method is None:
            logger.info(f' The gap interpolation method is updated: \
        {self.settings.gap["gaps_fill_settings"]["linear"]["method"]} --> {str(gap_interpolation_method)}')
            self.settings.gap["gaps_fill_settings"]["linear"]["method"] = str(gap_interpolation_method)

        if not gap_interpolation_max_consec_fill is None:
            logger.info(f' The gap max number of consecutive interpolations is updated: \
        {self.settings.gap["gaps_fill_settings"]["linear"]["max_consec_fill"]} --> {abs(int(gap_interpolation_max_consec_fill))}')
            self.settings.gap["gaps_fill_settings"]["linear"]["max_consec_fill"] = abs(int(gap_interpolation_max_consec_fill))


        # Gap debias fill
        if not gap_debias_prefered_leading_period_hours is None:
            logger.info(f' The size of the prefered leading period for debias gapfill is updated: \
        {self.settings.gap["gaps_fill_settings"]["model_debias"]["debias_period"]["prefered_leading_sample_duration_hours"]} --> {abs(int(gap_debias_prefered_leading_period_hours))}')
            self.settings.gap["gaps_fill_settings"]["model_debias"]["debias_period"]["prefered_leading_sample_duration_hours"] = abs(int(gap_debias_prefered_leading_period_hours))

        if not gap_debias_prefered_trailing_period_hours is None:
            logger.info(f' The size of the prefered trailing period for debias gapfill is updated: \
        {self.settings.gap["gaps_fill_settings"]["model_debias"]["debias_period"]["prefered_trailing_sample_duration_hours"]} --> {abs(int(gap_debias_prefered_trailing_period_hours))}')
            self.settings.gap["gaps_fill_settings"]["model_debias"]["debias_period"]["prefered_trailing_sample_duration_hours"] = abs(int(gap_debias_prefered_trailing_period_hours))

        if not gap_debias_minimum_leading_period_hours is None:
            logger.info(f' The minimum size of the leading period for debias gapfill is updated: \
        {self.settings.gap["gaps_fill_settings"]["model_debias"]["debias_period"]["minimum_leading_sample_duration_hours"]} --> {abs(int(gap_debias_minimum_leading_period_hours))}')
            self.settings.gap["gaps_fill_settings"]["model_debias"]["debias_period"]["minimum_leading_sample_duration_hours"] = abs(int(gap_debias_minimum_leading_period_hours))

        if not gap_debias_minimum_trailing_period_hours is None:
            logger.info(f' The minimum size of the trailing period for debias gapfill is updated: \
        {self.settings.gap["gaps_fill_settings"]["model_debias"]["debias_period"]["minimum_trailing_sample_duration_hours"]} --> {abs(int(gap_debias_minimum_trailing_period_hours))}')
            self.settings.gap["gaps_fill_settings"]["model_debias"]["debias_period"]["minimum_trailing_sample_duration_hours"] = abs(int(gap_debias_minimum_trailing_period_hours))

        # Gapfill automatic
        if not automatic_max_interpolation_duration_str is None:
            if is_timedelta(str(automatic_max_interpolation_duration_str)):
                logger.info(f' The maximum interpolation duration for automatic gapfill is updated: \
            {self.settings.gap["gaps_fill_settings"]["automatic"]["max_interpolation_duration_str"]} --> {str(automatic_max_interpolation_duration_str)}')
                self.settings.gap["gaps_fill_settings"]["automatic"]["max_interpolation_duration_str"] = str(automatic_max_interpolation_duration_str)
            else:
                logger.warning(f' {str(automatic_max_interpolation_duration_str)} is not a valid timedelta string. No update on this setting.')

        # Missing obs interpolation
        if not missing_obs_interpolation_method is None:
            logger.info(f' The missing observations interpolation method is updated: \
        {self.settings.missing_obs["missing_obs_fill_settings"]["linear"]["method"]} --> {str(missing_obs_interpolation_method)}')
            self.settings.missing_obs['missing_obs_fill_settings']['linear']['method'] = str(missing_obs_interpolation_method)



    def update_qc_settings(self, obstype='temp', gapsize_in_records=None, dupl_timestamp_keep=None, persis_time_win_to_check=None, persis_min_num_obs=None,
                            rep_max_valid_repetitions=None, gross_value_min_value=None, gross_value_max_value=None,
                            win_var_max_increase_per_sec=None, win_var_max_decrease_per_sec=None, win_var_time_win_to_check=None,
                            win_var_min_num_obs=None, step_max_increase_per_sec=None, step_max_decrease_per_sec=None):

        """
        Update the QC settings for the specified observation type.

        If a argument value is None, the default settings will not be updated.

        Parameters
        ----------
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

        Returns
        -------
        None.

        Note
        -------
        The gap defenition is independend of the observation type, and is thus set for
        all the observation types.

        """


        assert obstype in observation_types, f'{obstype} is not a known observation type'

        # Gap defenition
        if not gapsize_in_records is None:
            logger.info(f' The defenition of a gap (=gapsize) is updated: \
        {self.settings.gap["gaps_settings"]["gaps_finder"]["gapsize_n"]} --> {abs(int(gapsize_in_records))}')
            self.settings.gap["gaps_settings"]["gaps_finder"]["gapsize_n"] = abs(int(gapsize_in_records))

        # Gross value check
        if not gross_value_max_value is None:
            logger.info(f'Maximal value for gross value check updated: \
        {self.settings.qc["qc_check_settings"]["gross_value"][obstype]["max_value"]} --> {float(gross_value_max_value)}')
            self.settings.qc['qc_check_settings']["gross_value"][obstype]['max_value'] = float(gross_value_max_value)

        if not gross_value_min_value is None:
            logger.info(f'Manimal value for gross value check updated: \
        {self.settings.qc["qc_check_settings"]["gross_value"][obstype]["min_value"]} --> {float(gross_value_min_value)}')
            self.settings.qc['qc_check_settings']["gross_value"][obstype]['min_value'] = float(gross_value_min_value)

        # Duplicate check
        if not dupl_timestamp_keep is None:
            logger.info(f'Setting to keep (True) are remove (False) duplicate timestamps updated: \
        {self.settings.qc["qc_check_settings"]["duplicated_timestamp"]["keep"]} -->  {bool(dupl_timestamp_keep)}')
            self.settings.qc['qc_check_settings']["duplicated_timestamp"]['keep'] = bool(dupl_timestamp_keep)

        # Persistance check
        if not persis_time_win_to_check is None:
            if is_timedelta(str(persis_time_win_to_check)):
                logger.info(f'Time window size for persistance check updated:\
                {self.settings.qc["qc_check_settings"]["persistance"][obstype]["time_window_to_check"]}--> {str(persis_time_win_to_check)}')
                self.settings.qc['qc_check_settings']["persistance"][obstype]['time_window_to_check'] = str(persis_time_win_to_check)
            else:
                logger.warning(f' {str(persis_time_win_to_check)} is not a valid timedelta string. No update on this setting.')

        if not persis_min_num_obs is None:
            logger.info(f'Minimal window members for persistance check updated:\
            {self.settings.qc["qc_check_settings"]["persistance"][obstype]["min_num_obs"]} --> {abs(int(persis_min_num_obs))}')
            self.settings.qc['qc_check_settings']["persistance"][obstype]['min_num_obs'] = abs(int(persis_min_num_obs))

        # Repetitions check
        if not rep_max_valid_repetitions is None:
            logger.info(f'Maximal valid repetitions for repetitions check updated: \
                        {self.settings.qc["qc_check_settings"]["repetitions"][obstype]["max_valid_repetitions"]} --> {abs(int(rep_max_valid_repetitions))}')
            self.settings.qc['qc_check_settings']["repetitions"][obstype]['max_valid_repetitions'] = abs(int(rep_max_valid_repetitions))

        # Window variation check
        if not win_var_max_increase_per_sec is None:
            logger.info(f'Maximal increase per second for window variation check updated:\
                        {self.settings.qc["qc_check_settings"]["window_variation"][obstype]["max_increase_per_second"]} --> {abs(float(win_var_max_increase_per_sec))}')
            self.settings.qc['qc_check_settings']["window_variation"][obstype]['max_increase_per_second'] = abs(float(win_var_max_increase_per_sec))

        if not win_var_max_decrease_per_sec is None:
            logger.info(f'Maximal decrease per second for window variation check updated:\
                        {self.settings.qc["qc_check_settings"]["window_variation"][obstype]["max_decrease_per_second"] } --> {abs(float(win_var_max_decrease_per_sec))}')
            self.settings.qc['qc_check_settings']["window_variation"][obstype]['max_decrease_per_second'] = abs(float(win_var_max_decrease_per_sec))

        if not win_var_time_win_to_check is None:
            if is_timedelta(str(win_var_time_win_to_check)):
                logger.info(f'Time window for window variation check updated:\
                            {self.settings.qc["qc_check_settings"]["window_variation"][obstype]["time_window_to_check"]} --> {str(win_var_time_win_to_check)}')
                self.settings.qc['qc_check_settings']["window_variation"][obstype]['time_window_to_check'] = str(win_var_time_win_to_check)
            else:
                logger.warning(f' {str(persis_time_win_to_check)} is not a valid timedelta string. No update on this setting.')

        if not win_var_min_num_obs is None:
            logger.info(f'Minimal window members for window variation check updated:\
                         {self.settings.qc["qc_check_settings"]["window_variation"][obstype]["min_window_members"]}--> {abs(int(win_var_min_num_obs))}')
            self.settings.qc['qc_check_settings']["window_variation"][obstype]['min_window_members'] = abs(int(win_var_min_num_obs))

        # Step check
        if not step_max_increase_per_sec is None:
            logger.info(f'Maximal increase per second for step check updated:\
                        {self.settings.qc["qc_check_settings"]["step"][obstype]["max_increase_per_second"]}--> {abs(float(step_max_increase_per_sec))}')
            self.settings.qc['qc_check_settings']["step"][obstype]['max_increase_per_second'] = abs(float(step_max_increase_per_sec))

        if not step_max_decrease_per_sec is None:
            logger.info(f'Maximal decrease per second for step check updated:\
                       {self.settings.qc["qc_check_settings"]["step"][obstype]["max_decrease_per_second"]} --> {-1.0 * abs(float(step_max_decrease_per_sec))}')
            self.settings.qc['qc_check_settings']["step"][obstype]['max_decrease_per_second'] = -1.0 * abs(float(step_max_decrease_per_sec))


# =============================================================================
# dtype check functions
# =============================================================================
def is_timedelta(timedeltastr):
    try:
        pd.to_timedelta(timedeltastr)
        return True
    except:
        return False