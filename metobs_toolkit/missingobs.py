#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 16:12:41 2023

@author: thoverga
"""


import pandas as pd
import numpy as np
from datetime import datetime, timedelta
import logging

from metobs_toolkit.df_helpers import _find_closes_occuring_date

logger = logging.getLogger(__name__)

# =============================================================================
# Missing observation

# a missing observation is a missing timestamp
# =============================================================================


class Missingob_collection:
    """ Class object handling a set of missing observations. """
    def __init__(self, missing_obs_series):

        missing_obs_series.name = 'datetime'
        missing_obs_series.index.name = 'name'

        missing_obs_series = missing_obs_series.drop_duplicates()
        missing_obs_series = missing_obs_series.sort_index()

        missing_idx = missing_obs_series.reset_index()
        missing_idx = missing_idx.set_index(['name', 'datetime'])

        self.series = missing_obs_series
        self.idx = missing_idx.index

        # gap fill (only for conventional saving)
        self.fill_df = pd.DataFrame()
        self.fill_technique = None

    def __add__(self, other):
        comb_series = pd.concat([self.series, other.series])

        comb_series = comb_series.drop_duplicates()
        comb_series = comb_series.sort_index()

        self.series = comb_series
        comb_idx = comb_series.reset_index()
        comb_idx = comb_idx.set_index(['name', 'datetime'])
        self.idx = comb_idx.index
        return self


    def __str__(self):
        if self.series.empty:
            return f'Empty missing observations.'
        if not self.fill_df.empty:
            return f'Missing observations with filled ({self.fill_technique}) \
                values: \n {self.fill_df} \n Original missing observations on import: \n {self.idx}'

        return f'Missing observations: \n {self.series}'
    def __repr__(self):
        return self.__str__()


    def get_station_missingobs(self, name):
        """
        Get the missing observations of a specific station.

        Parameters
        ----------
        name : str
            The name of the station to extract the missing observation from.

        Returns
        -------
        Metobs_toolkit.Missingob_collection
            A subset of the missing observations from a specific station.

        """
        if name in self.series.index:
            return Missingob_collection(self.series.loc[[name]])
        else:
            # return empty collection
            series = pd.Series(data=[], name="datetime", dtype=object)
            series.index.name = "name"
            return Missingob_collection(series)

    def remove_missing_from_obs(self, obsdf):
        """
        Drop the missing observation records from an observational dataframe, if
        they are present.

        Parameters
        ----------
        obsdf : pandas.DataFrame
            Multiindex observational dataframe.

        Returns
        -------
        obsdf : pandas.DataFrame
            Multiindex observational dataframe without records linked to missing
            observations.

        """
        # Normally there are no missing records in the obsdf
        missing_multiidx = pd.MultiIndex.from_arrays(
            arrays=[self.series.index.to_list(), self.series.to_list()],
            names=["name", "datetime"],
        )

        obsdf = obsdf.drop(index=missing_multiidx, errors="ignore")

        return obsdf

    def interpolate_missing(self, obsdf, resolutionseries, obstype='temp', method='time'):
        """
        Fill the missing observations using an interpolation method.

        The "fill_df" and "fill_technique" attributes will be updated.

        Parameters
        ----------
        obsdf : Metobs_toolkit.Dataset.df
            The observations that can be used for the interpolation.
        resolutionseries : pd.Series
            The dataset resolution series for all stations..
        obstype : element of Metobs_toolkit.observational_types, optional
            Select which observation type you wish to interpolate. The default is 'temp'.
        method : valid input for pandas.DataFrame.interpolate method arg, optional
            Which interpolation method to use. The default is 'time'.

        Returns
        -------
        None.

        """
        # create fill column for the obstype
        self.fill_df[obstype] = np.nan
        self.fill_technique = 'interpolate'
        # locate the missing observation in observation space
        missing_obsspace = self.get_missing_indx_in_obs_space(obsdf, resolutionseries)

        # Set index for df fill attribute
        self.fill_df = pd.DataFrame(index=missing_obsspace)


        for staname, missingdt in missing_obsspace:
            staobs = obsdf.xs(staname, level='name')[obstype]
            # exclude nan values because they are no good leading/trailing
            staobs = staobs[~staobs.isnull()]
            print(f'staname: {staname}, missingdt: {missingdt}')
            # find leading and trailing datetimes
            leading_seconds =_find_closes_occuring_date(refdt = missingdt,
                                                        series_of_dt = staobs.index,
                                                        where='before')
            if np.isnan(leading_seconds):
                logger.warn(f'missing obs: {staname}, at {missingdt} does not have a leading timestamp.')
                continue

            leading_dt = missingdt - timedelta(seconds=leading_seconds)

            trailing_seconds =_find_closes_occuring_date(
                                                    refdt = missingdt,
                                                    series_of_dt = staobs.index,
                                                    where='after')
            if np.isnan(trailing_seconds):
                logger.warn(f'missing obs: {staname}, at {missingdt} does not have a trailing timestamp.')
                continue
            trailing_dt = missingdt + timedelta(seconds=trailing_seconds)

            # extract the values and combine them in a dataframe
            leading_val = staobs.loc[leading_dt]
            trailing_val = staobs.loc[trailing_dt]

            stadf = pd.DataFrame(
                index=[leading_dt, missingdt, trailing_dt],
                data={obstype: [leading_val, np.nan, trailing_val]}
            )


            # interpolate the missing obs
            stadf['interp'] = stadf[obstype].interpolate(
                                                method=method,
                )

            self.fill_df.loc[(staname, missingdt), obstype] = stadf.loc[missingdt, 'interp']




    def get_missing_indx_in_obs_space(self, obsdf, resolutionseries):
        """
        Function to found which missing timestamps are expected in the observation space.
        Because of time coarsening not all missing timestamps are expected in observation space.

        This function handles each station seperatly because stations can have differnent resolution/timerange.



        Parameters
        ----------
        obsdf : pandas.DataFrame()
            Dataset.df.
        resolutionseries : pandas.Series() or Timedelta
            Dataset.metadf['dataset_resolution'].

        Returns
        -------
        missing_obsspace : pandas.MultiIndex
            The multiindex (name - datetime) is returned with the missing timestamps that are expexted in the observation space.

        """

        missing_obsspace = pd.MultiIndex(
            levels=[["name"], ["datetime"]], codes=[[], []], names=["name", "datetime"]
        )

        # per stationtion because stations can have different resolutions/timerange
        for sta in self.series.index.unique():
            # Get missing observations in IO space
            sta_missing = self.series.loc[sta]
            if not isinstance(sta_missing, type(pd.Series(dtype=object))):
                sta_missing = pd.Series(data=[sta_missing], index=[sta], dtype=object)

            # Get start, end and frequency of the observation in obs space
            startdt = obsdf.xs(sta, level="name").index.min()
            enddt = obsdf.xs(sta, level="name").index.max()
            obs_freq = resolutionseries.loc[sta]

            # Make datetimerange
            obsrange = pd.date_range(
                start=startdt, end=enddt, freq=obs_freq, inclusive="both"
            )

            # # Look which missing timestamps appears obsspace
            sta_missing = sta_missing[sta_missing.isin(obsrange)]

            # Convert to multiindex
            if sta_missing.empty:
                continue

            sta_missing_idx = pd.MultiIndex.from_arrays(
                arrays=[[sta] * len(sta_missing), sta_missing.to_numpy()],
                names=["name", "datetime"],
            )

            missing_obsspace = missing_obsspace.append(sta_missing_idx)

        return missing_obsspace


