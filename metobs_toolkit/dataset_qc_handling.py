#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 13:44:49 2024

@author: thoverga
"""

import logging
from metobs_toolkit import Dataset


logger = logging.getLogger(__name__)

from metobs_toolkit.qc_checks import (
    gross_value_check,
    persistance_check,
    repetitions_check,
    # duplicate_timestamp_check,
    step_check,
    window_variation_check,
    # invalid_input_check,
    toolkit_buddy_check,
    titan_buddy_check,
    titan_sct_resistant_check,
)

from metobs_toolkit.plotting_functions import qc_stats_pie
from metobs_toolkit.qc_statistics import get_freq_statistics
from metobs_toolkit.df_helpers import (
    multiindexdf_datetime_subsetting,
    fmt_datetime_argument,
    init_multiindex,
    init_multiindexdf,
    init_triple_multiindexdf,
    metadf_to_gdf,
    conv_applied_qc_to_df,
    get_freqency_series,
    value_labeled_doubleidxdf_to_triple_idxdf,
    xs_save,
    concat_save,
)


class Dataset(Dataset):
    """Extension on the metobs_toolkit.Dataset class with QC related methods"""

    def get_qc_stats(self, obstype="temp", stationname=None, make_plot=True):
        """Get quality control statistics.

        Compute frequency statistics on the qc labels for an observationtype.
        The output is a dataframe containing the frequency statistics presented
        as percentages.

        These frequencies can also be presented as a collection of piecharts
        per check.

        With stationnames you can subset the data to one ore multiple stations.

        Parameters
        -----------
        obstype : str, optional
            Observation type to analyse the QC labels on. The default is
            'temp'.
        stationname : str, optional
            Stationname to subset the quality labels on. If None, all
            stations are used. The default is None.
        make_plot : Bool, optional
            If True, a plot with piecharts is generated. The default is True.

        Returns
        ---------
        dataset_qc_stats : pandas.DataFrame
            A table containing the label frequencies per check presented
            as percentages.

        Examples
        --------
        .. code-block:: python

            >>> import metobs_toolkit
            >>>
            >>> # Import data into a Dataset
            >>> dataset = metobs_toolkit.Dataset()
            >>> dataset.update_settings(
            ...                         input_data_file=metobs_toolkit.demo_datafile,
            ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
            ...                         template_file=metobs_toolkit.demo_template,
            ...                         )
            >>> dataset.import_data_from_file()
            >>> dataset.coarsen_time_resolution(freq='1h')
            >>>
            >>> # Apply quality control on the temperature observations
            >>> dataset.apply_quality_control(obstype='temp') #Using the default QC settings
            >>> dataset
            Dataset instance containing:
                 *28 stations
                 *['temp', 'humidity', 'wind_speed', 'wind_direction'] observation types
                 *10080 observation records
                 *1676 records labeled as outliers
                 *0 gaps
                 *3 missing observations
                 *records range: 2022-09-01 00:00:00+00:00 --> 2022-09-15 23:00:00+00:00 (total duration:  14 days 23:00:00)
                 *time zone of the records: UTC
                 *Coordinates are available for all stations.
            >>>
            >>> #Get quality control statistics
            >>> stats = dataset.get_qc_stats(make_plot=False)
            >>> stats
            ({'ok': 83.37301587301587, 'QC outliers': 16.6269841269...

        """
        # cobmine all and get final label
        comb_df = self.combine_all_to_obsspace()

        # subset to relevant columnt
        comb_df = xs_save(comb_df, obstype, level="obstype")[["label"]]

        # subset to stationnames
        if stationname is not None:
            assert stationname in comb_df.index.get_level_values(
                "name"
            ), f" stationnames: {stationname} is not a list."

            comb_df = comb_df.loc[stationname]

        # compute freq statistics
        final_freq, outl_freq, specific_freq = get_freq_statistics(
            comb_df=comb_df,
            obstype=obstype,
            checks_info=self.settings.qc["qc_checks_info"],
            gaps_info=self.settings.gap["gaps_info"],
            applied_qc_order=self._applied_qc,
        )

        if any([stat is None for stat in [final_freq, outl_freq, specific_freq]]):
            return None

        # make title
        orig_obstype = self.obstypes[obstype].get_orig_name()

        if stationname is None:
            title = f"Label frequency statistics on all stations for {orig_obstype}."
        else:
            title = f"Label frequency statistics for {stationname} for {orig_obstype}."

        if make_plot:
            # make pie plots
            qc_stats_pie(
                final_stats=final_freq,
                outlier_stats=outl_freq,
                specific_stats=specific_freq,
                plot_settings=self.settings.app["plot_settings"],
                qc_check_info=self.settings.qc["qc_checks_info"],
                title=title,
            )

        return (final_freq, outl_freq, specific_freq)

    def apply_quality_control(
        self,
        obstype="temp",
        gross_value=True,
        persistance=True,
        repetitions=True,
        step=True,
        window_variation=True,
    ):
        """Apply quality control methods to the dataset.

        The default settings are used, and can be changed in the
        settings_files/qc_settings.py

        The checks are performed in a sequence: gross_vallue -->
        persistance --> ..., Outliers by a previous check are ignored in the
        following checks!

        The dataset is updated inline.

        Parameters
        ----------
        obstype : String, optional
            Name of the observationtype you want to apply the checks on. The
            default is 'temp'.
        gross_value : Bool, optional
            If True the gross_value check is applied if False not. The default
            is True.
        persistance : Bool, optional
            If True the persistance check is applied if False not. The default
            is True.. The default is True.
        repetition : Bool, optional
            If True the repetations check is applied if False not. The default
            is True.
        step : Bool, optional
            If True the step check is applied if False not. The default is True.
        window_variation : Bool, optional
            If True the window_variation check is applied if False not. The
            default is True.

        Returns
        ---------
        None.

        Notes
        -----
        A schematic description of the quality control checks.

        Gross value check
        ==================
        This check looks for outliers based on unrealistic values

        1. Find observations that exceed a minimum and maximum value threshold.
        2. These observations are labeled as outliers.

        Persistence check
        =================
        Test observations to change over a specific period.

        1. Find the stations that have a maximum assumed observation frequency
           that does not exceed the minimum number of records for moving window
           size. The window size is defined by a duration.
        2. Subset to those stations.
        3. For each station, a moving window scan is applied that validates if
           there is variation in the observations (NaN's are excluded). The
           validation is only applied when a sufficient amount of records are
           found in the window specified by a threshold.
        4. After the scan, all records found in the windows without variation
           are labeled as outliers.

        Repetitions check
        =================
        Test if observation changes after a number of records.

        1. For each station, make a group of consecutive records for which
           the values do not change.
        2. Filter those groups that have more records than the maximum valid
           repetitions.
        3. All the records in these groups are labeled as outliers

        Note
        -----
          The repetitions check is similar to the persistence check, but not identical.
          The persistence check uses thresholds that are meteorologically based (i.g. the moving window is defined by a duration),
          in contrast to the repetitions check whose thresholds are instrumentally based (i.g. the "window" is defined by a number of records.)

        Step check
        ============
        Test if observations do not produce unphysical spikes in time series.

        1. Iterate over all the stations.
        2. Get the observations of the stations (i.g. drop the previously labeled outliers represented by NaN's).
        3. Find the observations for which:

           * The increase between two consecutive records is larger than the
             threshold. This threshold is defined by a maximum increase per second
             multiplied by the timedelta (in seconds) between the consecutive
             records.
           * Similar filter for a decrease.
        4. The found observations are labeled as outliers.

        Note
        -----
          In general, for temperatures,  the decrease threshold is set less stringent than the increase
          threshold. This is because a temperature drop is meteorologycally more
          common than a sudden increase which is often the result of a radiation error.

        Window Variation check
        =======================
        Test if the variation is found in a moving window.

        1. Find the stations that have a maximum assumed observation frequency
           that does not exceed the minimum number of records for moving window
           size. The window size is defined by a duration.
        2. Compute the maximum increase and decrease thresholds for a window.
           This is done by multiplying the maximum increase per second by the
           window size in seconds.
        3. For each station, a moving window scan is applied that validates if
           the maximum increase/decrease thresholds are exceeded. This is done
           by comparison of the minimum and maximum values inside the window. The
           validation is only applied when a sufficient amount of records are
           found in the window specified by a threshold.
        4. After the scan, *all* records found in the window that exceed one
           of these thresholds are labeled as outliers.


        Examples
        --------
        .. code-block:: python

            >>> import metobs_toolkit
            >>>
            >>> # Import data into a Dataset
            >>> dataset = metobs_toolkit.Dataset()
            >>> dataset.update_settings(
            ...                         input_data_file=metobs_toolkit.demo_datafile,
            ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
            ...                         template_file=metobs_toolkit.demo_template,
            ...                         )
            >>> dataset.import_data_from_file()
            >>> dataset.coarsen_time_resolution(freq='1h')
            >>>
            >>> #Update some temperature QC settings
            >>> dataset.update_qc_settings(obstype='temp',
            ...                            gross_value_max_value=42.,
            ...                            persis_time_win_to_check='4h',
            ...                            buddy_min_std = 1.5)

            >>> # Apply quality control on the temperature observations
            >>> dataset.apply_quality_control(obstype='temp')
            >>> dataset
            Dataset instance containing:
                 *28 stations
                 *['temp', 'humidity', 'wind_speed', 'wind_direction'] observation types
                 *10080 observation records
                 *1676 records labeled as outliers
                 *0 gaps
                 *3 missing observations
                 *records range: 2022-09-01 00:00:00+00:00 --> 2022-09-15 23:00:00+00:00 (total duration:  14 days 23:00:00)
                 *time zone of the records: UTC
                 *Coordinates are available for all stations.

        """
        if repetitions:
            apliable = _can_qc_be_applied(self, obstype, "repetitions")
            if apliable:
                logger.info("Applying repetitions check.")

                obsdf, outl_df = repetitions_check(
                    obsdf=self.df,
                    obstype=obstype,
                    checks_info=self.settings.qc["qc_checks_info"],
                    checks_settings=self.settings.qc["qc_check_settings"],
                )

                # update the dataset and outliers
                self.df = obsdf
                if not outl_df.empty:
                    self.outliersdf = concat_save([self.outliersdf, outl_df])

                # add this check to the applied checks
                self._applied_qc = concat_save(
                    [
                        self._applied_qc,
                        conv_applied_qc_to_df(
                            obstypes=obstype, ordered_checknames="repetitions"
                        ),
                    ],
                    ignore_index=True,
                )

        if gross_value:
            apliable = _can_qc_be_applied(self, obstype, "gross_value")

            if apliable:
                logger.info("Applying gross value check.")

                obsdf, outl_df = gross_value_check(
                    obsdf=self.df,
                    obstype=obstype,
                    checks_info=self.settings.qc["qc_checks_info"],
                    checks_settings=self.settings.qc["qc_check_settings"],
                )

                # update the dataset and outliers
                self.df = obsdf
                if not outl_df.empty:
                    self.outliersdf = concat_save([self.outliersdf, outl_df])

                # add this check to the applied checks
                self._applied_qc = concat_save(
                    [
                        self._applied_qc,
                        conv_applied_qc_to_df(
                            obstypes=obstype, ordered_checknames="gross_value"
                        ),
                    ],
                    ignore_index=True,
                )

        if persistance:
            apliable = _can_qc_be_applied(self, obstype, "persistance")

            if apliable:
                logger.info("Applying persistance check.")
                obsdf, outl_df = persistance_check(
                    station_frequencies=self.metadf["dataset_resolution"],
                    obsdf=self.df,
                    obstype=obstype,
                    checks_info=self.settings.qc["qc_checks_info"],
                    checks_settings=self.settings.qc["qc_check_settings"],
                )

                # update the dataset and outliers
                self.df = obsdf
                if not outl_df.empty:
                    self.outliersdf = concat_save([self.outliersdf, outl_df])

                # add this check to the applied checks
                self._applied_qc = concat_save(
                    [
                        self._applied_qc,
                        conv_applied_qc_to_df(
                            obstypes=obstype, ordered_checknames="persistance"
                        ),
                    ],
                    ignore_index=True,
                )

        if step:
            apliable = _can_qc_be_applied(self, obstype, "step")

            if apliable:
                logger.info("Applying step-check.")
                obsdf, outl_df = step_check(
                    obsdf=self.df,
                    obstype=obstype,
                    checks_info=self.settings.qc["qc_checks_info"],
                    checks_settings=self.settings.qc["qc_check_settings"],
                )

                # update the dataset and outliers
                self.df = obsdf
                if not outl_df.empty:
                    self.outliersdf = concat_save([self.outliersdf, outl_df])

                # add this check to the applied checks
                self._applied_qc = concat_save(
                    [
                        self._applied_qc,
                        conv_applied_qc_to_df(
                            obstypes=obstype, ordered_checknames="step"
                        ),
                    ],
                    ignore_index=True,
                )

        if window_variation:
            apliable = _can_qc_be_applied(self, obstype, "window_variation")
            if apliable:
                logger.info("Applying window variation-check.")
                obsdf, outl_df = window_variation_check(
                    station_frequencies=self.metadf["dataset_resolution"],
                    obsdf=self.df,
                    obstype=obstype,
                    checks_info=self.settings.qc["qc_checks_info"],
                    checks_settings=self.settings.qc["qc_check_settings"],
                )

                # update the dataset and outliers
                self.df = obsdf
                if not outl_df.empty:
                    self.outliersdf = concat_save([self.outliersdf, outl_df])

                # add this check to the applied checks
                self._applied_qc = concat_save(
                    [
                        self._applied_qc,
                        conv_applied_qc_to_df(
                            obstypes=obstype,
                            ordered_checknames="window_variation",
                        ),
                    ],
                    ignore_index=True,
                )

        self._qc_checked_obstypes.append(obstype)
        self._qc_checked_obstypes = list(set(self._qc_checked_obstypes))
        self.outliersdf = self.outliersdf.sort_index()

    def apply_buddy_check(
        self,
        obstype="temp",
        use_constant_altitude=False,
        haversine_approx=True,
        metric_epsg="31370",
    ):
        """Apply the buddy check on the observations.

        The buddy check compares an observation against its neighbours (i.e.
        buddies). The check looks for buddies in a neighbourhood specified by
        a certain radius. The buddy check flags observations if the
        (absolute value of the) difference between the observations and the
        average of the neighbours normalized by the standard deviation in the
        circle is greater than a predefined threshold.

        This check is based on the buddy check from titanlib. Documentation on
        the titanlib buddy check can be found
        `here <https://github.com/metno/titanlib/wiki/Buddy-check>`_.


        The observation and outliers attributes will be updated accordingly.

        Parameters
        ----------
        obstype : String, optional
            Name of the observationtype you want to apply the checks on. The
            default is 'temp'.
        use_constant_altitude : bool, optional
            Use a constant altitude for all stations. The default is False.
        haversine_approx : bool, optional
            Use the haversine approximation (earth is a sphere) to calculate
            distances between stations. The default is True.
        metric_epsg : str, optional
            EPSG code for the metric CRS to calculate distances in. Only used when
            haversine approximation is set to False. Thus becoming a better
            distance approximation but not global applicable The default is '31370'
            (which is suitable for Belgium).

        Returns
        -------
        None.

        Notes
        -----
        A schematic step-by-step description of the buddy check:

        1. A distance matrix is constructed for all inter distances between the stations. This is done using the haversine approximation, or by first converting the Coordinate Reference System (CRS) to a metric one, specified by an EPSG code.
        2. A set of all (spatial) buddies per station is created by filtering out all stations that are too far.
        3. The buddies are further filtered based on altitude differences with respect to the reference station.
        4. For each station:

           * Observations of buddies are extracted from all observations.
           * These observations are corrected for altitude differences by assuming a constant lapse rate.
           * For each reference record, the mean, standard deviation (std), and sample size of the corrected buddies’ observations are computed.
           * If the std is lower than the minimum std, it is replaced by the minimum std.
           * Chi values are calculated for all reference records.
           * If the Chi value is larger than the std_threshold, the record is accepted, otherwise it is marked as an outlier.

        Examples
        --------
        .. code-block:: python

            >>> import metobs_toolkit
            >>>
            >>> # Import data into a Dataset
            >>> dataset = metobs_toolkit.Dataset()
            >>> dataset.update_settings(
            ...                         input_data_file=metobs_toolkit.demo_datafile,
            ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
            ...                         template_file=metobs_toolkit.demo_template,
            ...                         )
            >>> dataset.import_data_from_file()
            >>> dataset.coarsen_time_resolution(freq='1h')
            >>>
            >>> #Update some temperature QC settings
            >>> dataset.update_qc_settings(obstype='temp',
            ...                            buddy_min_std=1.5,
            ...                            buddy_threshold=3.2)

            >>> # Apply buddy check on the temperature observations
            >>> dataset.apply_buddy_check(obstype='temp',
            ...                           use_constant_altitude=True)
            >>> dataset
            Dataset instance containing:
                 *28 stations
                 *['temp', 'humidity', 'wind_speed', 'wind_direction'] observation types
                 *10080 observation records
                 *69 records labeled as outliers
                 *0 gaps
                 *3 missing observations
                 *records range: 2022-09-01 00:00:00+00:00 --> 2022-09-15 23:00:00+00:00 (total duration:  14 days 23:00:00)
                 *time zone of the records: UTC
                 *Coordinates are available for all stations.
        """

        logger.info("Applying the toolkit buddy check")

        checkname = "buddy_check"

        # 1. coordinates are available?
        if self.metadf["lat"].isnull().any():
            logger.warning(
                f"Not all coordinates are available, the {checkname} cannot be executed!"
            )
            return
        if self.metadf["lon"].isnull().any():
            logger.warning(
                f"Not all coordinates are available, the {checkname} cannot be executed!"
            )
            return

        # set constant altitude if needed:

        # if altitude is already available, save it to restore it after this check
        restore_altitude = False
        if use_constant_altitude:
            if "altitulde" in self.metadf.columns:
                self.metadf["altitude_backup"] = self.metadf["altitude"]
                restore_altitude = True

            self.metadf["altitude"] = 2.0  # absolut value does not matter

        # 2. altitude available?
        if (not use_constant_altitude) & ("altitude" not in self.metadf.columns):
            logger.warning(
                f"The altitude is not known for all stations. The {checkname} cannot be executed!"
            )
            logger.info(
                '(To resolve this error you can: \n *Use the Dataset.get_altitude() method \n *Set use_constant_altitude to True \n update the "altitude" column in the metadf attribute of your Dataset.'
            )
            return
        if (not use_constant_altitude) & (self.metadf["altitude"].isnull().any()):
            logger.warning(
                f"The altitude is not known for all stations. The {checkname} cannot be executed!"
            )
            logger.info(
                '(To resolve this error you can: \n *Use the Dataset.get_altitude() method \n *Set use_constant_altitude to True \n *Update the "altitude" column in the metadf attribute of your Dataset.)'
            )
            return

        apliable = _can_qc_be_applied(self, obstype, checkname)
        if apliable:
            buddy_set = self.settings.qc["qc_check_settings"][checkname][obstype]
            outl_flag = self.settings.qc["qc_checks_info"][checkname]["outlier_flag"]
            obsdf, outliersdf = toolkit_buddy_check(
                obsdf=self.df,
                metadf=self.metadf,
                obstype=obstype,
                buddy_radius=buddy_set["radius"],
                min_sample_size=buddy_set["num_min"],
                max_alt_diff=buddy_set["max_elev_diff"],
                min_std=buddy_set["min_std"],
                std_threshold=buddy_set["threshold"],
                metric_epsg=metric_epsg,
                lapserate=buddy_set["elev_gradient"],
                outl_flag=outl_flag,
                haversine_approx=haversine_approx,
            )

            # update the dataset and outliers
            self.df = obsdf
            if not outliersdf.empty:
                self.outliersdf = concat_save([self.outliersdf, outliersdf])

            # add this check to the applied checks
            self._applied_qc = concat_save(
                [
                    self._applied_qc,
                    conv_applied_qc_to_df(
                        obstypes=obstype, ordered_checknames=checkname
                    ),
                ],
                ignore_index=True,
            )

        else:
            logger.warning(
                f"The {checkname} can NOT be applied on {obstype} because it was already applied on this observation type!"
            )

        # Revert artificial data that has been added if needed
        if restore_altitude:  # altitude was overwritten, thus revert it
            self.metadf["altitude"] = self.metadf["altitude_backup"]
            self.metadf = self.metadf.drop(columns=["altitude_backup"])

        elif use_constant_altitude:
            # when no alitude was available apriori, remove the fake constant altitude column
            self.metadf = self.metadf.drop(columns=["altitude"])

    def apply_titan_buddy_check(self, obstype="temp", use_constant_altitude=False):
        """Apply the TITAN buddy check on the observations.

        The buddy check compares an observation against its neighbours (i.e. buddies). The check looks for
        buddies in a neighbourhood specified by a certain radius. The buddy check flags observations if the
        (absolute value of the) difference between the observations and the average of the neighbours
        normalized by the standard deviation in the circle is greater than a predefined threshold.

        See the `titanlib documentation on the buddy check <https://github.com/metno/titanlib/wiki/Buddy-check>`_
        for futher details.

        The observation and outliers attributes will be updated accordingly.

        Parameters
        ----------
        obstype : String, optional
            Name of the observationtype you want to apply the checks on. The
            default is 'temp'.
        use_constant_altitude : bool, optional
            Use a constant altitude for all stations. The default is False.

        Returns
        -------
        None.

        Note
        -------
        To update the check settings, use the update_titan_qc_settings method
        of the Dataset class.

        Warning
        --------
        To use this method, you must install titanlib. Windows users must have
        a c++ compiler installed. See the titanlib documentation: https://github.com/metno/titanlib/wiki/Installation.

        Examples
        --------
        .. code-block:: python

            >>> import metobs_toolkit
            >>>
            >>> # Import data into a Dataset
            >>> dataset = metobs_toolkit.Dataset()
            >>> dataset.update_settings(
            ...                         input_data_file=metobs_toolkit.demo_datafile,
            ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
            ...                         template_file=metobs_toolkit.demo_template,
            ...                         )
            >>> dataset.import_data_from_file()
            >>> dataset.coarsen_time_resolution(freq='1h')
            >>>
            >>> #Update some temperature QC settings
            >>> dataset.update_titan_qc_settings(obstype='temp',
            ...                                  buddy_min_std=1.5,
            ...                                  buddy_threshold=3.2,
            ...                                  buddy_num_min=5)
            buddy num min for the TITAN buddy check updated:  2--> 5
            buddy threshold for the TITAN buddy check updated:  1.5--> 3.2
            buddy min std for the TITAN buddy check updated:  1.0--> 1.5
            >>> # Apply buddy check on the temperature observations
            >>> dataset.apply_titan_buddy_check(obstype='temp',
            ...                                 use_constant_altitude=True)
            >>> dataset
            Dataset instance containing:
                 *28 stations
                 *['temp', 'humidity', 'wind_speed', 'wind_direction'] observation types
                 *10080 observation records
                 *35 records labeled as outliers
                 *0 gaps
                 *3 missing observations
                 *records range: 2022-09-01 00:00:00+00:00 --> 2022-09-15 23:00:00+00:00 (total duration:  14 days 23:00:00)
                 *time zone of the records: UTC
                 *Coordinates are available for all stations.

        """
        logger.info("Applying the titan buddy check")

        try:
            import titanlib

            # Add version restrictions??
        except ModuleNotFoundError:
            logger.warning(
                "Titanlib is not installed, install it manually if you want to use this functionallity."
            )
            return

        checkname = "titan_buddy_check"

        # 1. coordinates are available?
        if self.metadf["lat"].isnull().any():
            logger.warning(
                f"Not all coordinates are available, the {checkname} cannot be executed!"
            )
            return
        if self.metadf["lon"].isnull().any():
            logger.warning(
                f"Not all coordinates are available, the {checkname} cannot be executed!"
            )
            return

        # set constant altitude if needed:

        # if altitude is already available, save it to restore it after this check
        restore_altitude = False
        if use_constant_altitude:
            if "altitulde" in self.metadf.columns:
                self.metadf["altitude_backup"] = self.metadf["altitude"]
                restore_altitude = True

            self.metadf["altitude"] = 2.0  # absolut value does not matter

        # 2. altitude available?
        if (not use_constant_altitude) & ("altitude" not in self.metadf.columns):
            logger.warning(
                f"The altitude is not known for all stations. The {checkname} cannot be executed!"
            )
            logger.info(
                '(To resolve this error you can: \n *Use the Dataset.get_altitude() method \n *Set use_constant_altitude to True \n update the "altitude" column in the metadf attribute of your Dataset.'
            )
            return
        if (not use_constant_altitude) & (self.metadf["altitude"].isnull().any()):
            logger.warning(
                f"The altitude is not known for all stations. The {checkname} cannot be executed!"
            )
            logger.info(
                '(To resolve this error you can: \n *Use the Dataset.get_altitude() method \n *Set use_constant_altitude to True \n *Update the "altitude" column in the metadf attribute of your Dataset.)'
            )
            return

        apliable = _can_qc_be_applied(self, obstype, checkname)
        if apliable:
            obsdf, outliersdf = titan_buddy_check(
                obsdf=self.df,
                metadf=self.metadf,
                obstype=obstype,
                checks_info=self.settings.qc["qc_checks_info"],
                checks_settings=self.settings.qc["titan_check_settings"][checkname][
                    obstype
                ],
                titan_specific_labeler=self.settings.qc["titan_specific_labeler"][
                    checkname
                ],
            )

            # update the dataset and outliers
            self.df = obsdf
            if not outliersdf.empty:
                self.outliersdf = concat_save([self.outliersdf, outliersdf])

            # add this check to the applied checks
            self._applied_qc = concat_save(
                [
                    self._applied_qc,
                    conv_applied_qc_to_df(
                        obstypes=obstype, ordered_checknames=checkname
                    ),
                ],
                ignore_index=True,
            )

        else:
            logger.warning(
                f"The {checkname} can NOT be applied on {obstype} because it was already applied on this observation type!"
            )

        # Revert artificial data that has been added if needed
        if restore_altitude:  # altitude was overwritten, thus revert it
            self.metadf["altitude"] = self.metadf["altitude_backup"]
            self.metadf = self.metadf.drop(columns=["altitude_backup"])

        elif use_constant_altitude:
            # when no alitude was available apriori, remove the fake constant altitude column
            self.metadf = self.metadf.drop(columns=["altitude"])

    def apply_titan_sct_resistant_check(self, obstype="temp"):
        """Apply the TITAN spatial consistency test (resistant).

        The SCT resistant check is a spatial consistency check which compares each observations to what is expected given the other observations in the
        nearby area. If the deviation is large, the observation is removed. The SCT uses optimal interpolation
        (OI) to compute an expected value for each observation. The background for the OI is computed from
        a general vertical profile of observations in the area.

        See the `titanlib documentation on the sct check <https://github.com/metno/titanlib/wiki/Spatial-consistency-test-resistant>`_
        for futher details.

        The observation and outliers attributes will be updated accordingly.


        Parameters
        ----------
        obstype : String, optional
            Name of the observationtype you want to apply the checks on. The
            default is 'temp'.

        Returns
        -------
        None.

        Note
        -------
        To update the check settings, use the update_titan_qc_settings method
        of the Dataset class.

        Warning
        --------
        To use this method, you must install titanlib. Windows users must have
        a c++ compiler installed. See the titanlib documentation: https://github.com/metno/titanlib/wiki/Installation.

        Warning
        -------
        This method is a python wrapper on titanlib c++ scripts, and it is prone
        to segmentation faults. The perfomance of this check is thus not
        guaranteed!

        Examples
        --------
        .. code-block:: python

             import metobs_toolkit

             # Import data into a Dataset
             dataset = metobs_toolkit.Dataset()
             dataset.update_settings(
                                     input_data_file=metobs_toolkit.demo_datafile,
                                     input_metadata_file=metobs_toolkit.demo_metadatafile,
                                     template_file=metobs_toolkit.demo_template,
                                     )
             dataset.import_data_from_file()
             dataset.coarsen_time_resolution(freq='1h')

             #Get altitude of all stations
             dataset.get_altitude()

             #Update some temperature QC settings
             dataset.update_titan_qc_settings(obstype='temp',
                                              sct_outer_radius=25000)


             # Apply buddy check on the temperature observations
             dataset.apply_titan_sct_resistant_check(obstype='temp')

        """
        logger.info("Applying the titan SCT check")

        try:
            import titanlib

            # Add version restrictions??
        except ModuleNotFoundError:
            logger.warning(
                "Titanlib is not installed, install it manually if you want to use this functionallity."
            )
            return

        checkname = "titan_sct_resistant_check"
        # check if required metadata is available:

        # 1. coordinates are available?
        if self.metadf["lat"].isnull().any():
            logger.warning(
                f"Not all coordinates are available, the {checkname} cannot be executed!"
            )
            return
        if self.metadf["lon"].isnull().any():
            logger.warning(
                f"Not all coordinates are available, the {checkname} cannot be executed!"
            )
            return

        # 2. altitude available?
        if "altitude" not in self.metadf.columns:
            logger.warning(
                f"The altitude is not known for all stations. The {checkname} cannot be executed!"
            )
            logger.info(
                '(To resolve this error you can: \n *Use the Dataset.get_altitude() method \n *Set use_constant_altitude to True \n update the "altitude" column in the metadf attribute of your Dataset.'
            )
            return
        if self.metadf["altitude"].isnull().any():
            logger.warning(
                f"The altitude is not known for all stations. The {checkname} cannot be executed!"
            )
            logger.info(
                '(To resolve this error you can: \n *Use the Dataset.get_altitude() method \n *Set use_constant_altitude to True \n *Update the "altitude" column in the metadf attribute of your Dataset.)'
            )
            return

        apliable = _can_qc_be_applied(self, obstype, checkname)
        if apliable:
            obsdf, outliersdf = titan_sct_resistant_check(
                obsdf=self.df,
                metadf=self.metadf,
                obstype=obstype,
                checks_info=self.settings.qc["qc_checks_info"],
                checks_settings=self.settings.qc["titan_check_settings"][checkname][
                    obstype
                ],
                titan_specific_labeler=self.settings.qc["titan_specific_labeler"][
                    checkname
                ],
            )

            # update the dataset and outliers
            self.df = obsdf
            if not outliersdf.empty:
                self.outliersdf = concat_save([self.outliersdf, outliersdf])

            # add this check to the applied checks
            self._applied_qc = concat_save(
                [
                    self._applied_qc,
                    conv_applied_qc_to_df(
                        obstypes=obstype, ordered_checknames=checkname
                    ),
                ],
                ignore_index=True,
            )

        else:
            logger.warning(
                f"The {checkname} can NOT be applied on {obstype} because it was already applied on this observation type!"
            )


# =============================================================================
# Helpers
# =============================================================================


def _can_qc_be_applied(dataset, obstype, checkname):
    """Test if a qc check can be applied."""
    # test if check is already applied on the obstype
    applied_df = dataset._applied_qc
    can_be_applied = (
        not applied_df[
            (applied_df["obstype"] == obstype) & (applied_df["checkname"] == checkname)
        ].shape[0]
        > 0
    )

    if not can_be_applied:
        logger.warning(
            f"The {checkname} check can NOT be applied on {obstype} because it was already applied on this observation type!"
        )
        return False
    # test of all settings are present for the check on the obstype
    if checkname not in [
        "duplicated_timestamp",
        "titan_buddy_check",
        "titan_sct_resistant_check",
    ]:
        # these checks are obstype depending,
        required_keys = list(
            dataset.settings.qc["qc_check_settings"][checkname]["temp"].keys()
        )  # use temp to find all required settings
        if obstype not in dataset.settings.qc["qc_check_settings"][checkname].keys():
            logger.warning(
                f"The {checkname} check can NOT be applied on {obstype} because none of the required check settings are found. The following are missing: {required_keys}"
            )
            return False

        if not all(
            [
                req_key
                in dataset.settings.qc["qc_check_settings"][checkname][obstype].keys()
                for req_key in required_keys
            ]
        ):
            # not all required settings are available
            missing_settings = [
                req_key
                for req_key in required_keys
                if req_key
                not in dataset.settings.qc["qc_check_settings"][checkname][
                    obstype
                ].keys()
            ]
            logger.warning(
                f"The {checkname} check can NOT be applied on {obstype} because not all required check settings ar found. The following are missing: {missing_settings}"
            )
            return False

    return True
