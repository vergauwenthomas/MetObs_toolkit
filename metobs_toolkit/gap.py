#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Specific classes are created because gaps and missing obs can have out-of-sync
datetimes wrt dataset.df.

Created on Fri Mar  3 09:15:56 2023

@author: thoverga
"""

import pandas as pd
import numpy as np
import logging
from datetime import datetime, timedelta
import math


from metobs_toolkit.gap_filling import (
    interpolate_gap,
    create_leading_trailing_debias_periods,
    make_era_bias_correction,
)

from metobs_toolkit.df_helpers import (
    format_outliersdf_to_doubleidx,
    get_likely_frequency,
    _find_closes_occuring_date
)

from metobs_toolkit.df_helpers import init_multiindex, init_multiindexdf

from metobs_toolkit.missingobs import Missingob_collection

logger = logging.getLogger(__name__)


# =============================================================================
# Gap class

# a gap is a sequence of repeting missing obs
# =============================================================================


class Gap:
    """Gap class holds all gap information and methods for gaps."""

    def __init__(self, name, startdt, enddt):
        """
        Initiate Gap object with a name, startdt and enddt.

        Parameters
        ----------
        name : String
            Station name where the gap occures.
        startdt : datetime.datetime
            Start datetime of the gap (included).
        enddt : datetime.datetime
            End datetime of the gap (included)

        Returns
        -------
        None

        """
        # init attributes
        self.name = name
        self.startgap = startdt  # in IO space
        self.endgap = enddt  # in IO space
        self.duration = enddt - startdt

        # computed attributes
        self.leading_timestamp = None  # last ob_dt before gap in datset space
        self.leading_val = {} #keys are obstypes
        self.trailing_timestamp = None  # first ob_dt after gap in dataset space
        self.trailing_val = {} #keys are obstypes

        self.exp_gap_idx = None

        # gap fill (only for conventional saving)
        self.gapfill_df = pd.DataFrame() #index: datetime, columns: obstypes, values: fill_values
        self.gapfill_technique = None #will become a string
        self.gapfill_info = None #only for the user
        self.gapfill_errormessage = {} #keys are obstypes


    def __str__(self):
        return f"Gap instance of {self.name} for {self.startgap} --> {self.endgap}, duration: {self.duration}"
    def __repr__(self):
        return self.__str__()

    def get_info(self):
        print(f'Gap for {self.name} with: \n')
        print(f'\n ---- Gap info ----- \n')
        print(f'  * Start gap: {self.startgap} \n')
        print(f'  * End gap: {self.endgap} \n')
        print(f'  * Duration gap: {self.duration} \n')
        print(f'\n ---- Gap fill info ----- \n')

        obstypes = self.gapfill_df.columns.to_list()
        if self.gapfill_df.empty:
            print ('(No gapfill applied)')
        elif self.gapfill_technique == 'interpolation':
            for obstype in obstypes:
                print(f'  * On observation type: {obstype}')
                print(f'  * Technique: {self.gapfill_technique} \n')
                print(f'  * Leading timestamp: {self.leading_timestamp} with  {obstype} = {self.leading_val[obstype]}\n')
                print(f'  * Trailing timestamp: {self.trailing_timestamp} with  {obstype} = {self.trailing_val[obstype]}\n')
                print(f'  * Filled values: {self.gapfill_df[obstype]} \n')
                if obstype in self.gapfill_errormessage:
                    print(f'  * Gapfill message: {self.gapfill_errormessage[obstype]} \n')


        elif self.gapfill_technique == "debias gapfill":
            for obstype in obstypes:
                print(f'  * On observation type: {obstype}')
                print(f'  * Technique: {self.gapfill_technique} \n')
                # print(f'  * Leading timestamp: {self.leading_timestamp} with  {obstype} = {self.leading_val[obstype]}\n')
                # print(f'  * Trailing timestamp: {self.trailing_timestamp} with  {obstype} = {self.trailing_val[obstype]}\n')
                print(f'  * Filled values: {self.gapfill_df[obstype]} \n')
                if obstype in self.gapfill_errormessage:
                    print(f'  * Gapfill message: {self.gapfill_errormessage[obstype]} \n')
        else:
            print('technique not implemented in yet in show')



    def to_df(self):
        """
        Convert a Gap object to a dataframe (with one row). The station name is
        the index and two colums ('start_gap', 'end_gap') are constructed.

        Returns
        -------
        pandas.DataFrame()
            Gap in dataframe format.

        """
        return pd.DataFrame(
            index=[self.name],
            data={"start_gap": self.startgap,
                  "end_gap": self.endgap,
                  "duration": self.duration}
        )

    def update_leading_trailing_obs(self, obsdf, outliersdf, obs_only=False):
        """
        Add the leading (last obs before gap) and trailing (first obs after gap)
        as extra columns to the self.df.

        One can specify to look for leading and trailing in the obsdf or in both
        the obsdf and outliersdf.

        The gap leading and trailing timestamps and value attributes are updated.

        If no leading/trailing timestamp is found, it is set to the gaps startdt/enddt.

        Parameters
        ----------
        obsdf : pandas.DataFrame
            Dataset.df
        outliersdf : pandas.DataFrame
            Dataset.outliersdf
        obs_only: bool, optional
            If True, only the obsdf will be used to search for leading and trailing.

        Returns
        -------
        None.

        """

        sta_obs = obsdf.xs(self.name, level="name").index
        if obs_only:
            sta_comb = sta_obs
        else:
            outliersdf = format_outliersdf_to_doubleidx(outliersdf)

            # combine timestamps of observations and outliers
            sta_outl = outliersdf.xs(self.name, level="name").index
            sta_comb = sta_obs.append(sta_outl)

        # find minimium timediff before
        before_diff = _find_closes_occuring_date(
            refdt=self.startgap, series_of_dt=sta_comb, where="before"
        )

        # if no timestamps are before gap, assume gap at the start of the observations
        if math.isnan(before_diff):
            before_diff = 0.0

        # find minimum timediff after gap
        after_diff = _find_closes_occuring_date(
            refdt=self.endgap, series_of_dt=sta_comb, where="after"
        )
        # if no timestamps are after gap, assume gap at the end of the observations
        if math.isnan(after_diff):
            after_diff = 0.0

        # get before and after timestamps
        self.leading_timestamp = self.startgap - timedelta(seconds=before_diff)
        self.trailing_timestamp = self.endgap + timedelta(seconds=after_diff)

        # get the values
        try:
            self.leading_val = obsdf.loc[(self.name, self.leading_timestamp)].to_dict()
        except KeyError:
            print('LEADING VAL NOT IN OBSDF --> THIS IS NOT WHAT YOU WHANT I THINK ; FIX THIS')
            self.leading_val = {}
        try:
            self.trailing_val = obsdf.loc[(self.name, self.trailing_timestamp)].to_dict()
        except KeyError:
            print('LEADING VAL NOT IN OBSDF --> THIS IS NOT WHAT YOU WHANT I THINK ; FIX THIS')
            self.trailing_val = {}



    def update_gaps_indx_in_obs_space(self, obsdf, outliersdf, dataset_res):
        """

        Explode the gap, to the dataset resolution and format to a multiindex
        with name -- datetime.

        In addition the last observation before the gap (leading), and first
        observation (after) the gap are computed and stored in the df attribute.
        (the outliers are used to look for leading and trailing observations.)


        Parameters
        ----------
        obsdf : Dataset.df
            The Dataset.df attribute. (Needed to extract trailing/leading
                                       observations.)
        outliersdf : Dataset.outliersdf
            The Dataset.outliersdf attribute.(Needed to extract trailing/leading
                                              observations.))
        resolutionseries : Datetime.timedelta
            Resolution of the station observations in the dataset.

        Returns
        -------
        None

        """

        outliersdf = format_outliersdf_to_doubleidx(outliersdf)
        self.update_leading_trailing_obs(obsdf, outliersdf)

        gaprange = pd.date_range(
            start=self.leading_timestamp,
            end=self.trailing_timestamp,
            freq=dataset_res,
            inclusive="neither",
        )

        self.exp_gap_idx = pd.MultiIndex.from_arrays(
            arrays=[[self.name] * len(gaprange), gaprange], names=["name", "datetime"]
        )

    # =============================================================================
    #         Gapfill
    # =============================================================================

    def apply_interpolate_gap(
        self,
        obsdf,
        outliersdf,
        dataset_res,
        obstype="temp",
        method="time",
        max_consec_fill=100,
    ):
        """
        Fill a Gap using a linear interpolation gapfill method for an obstype.

        The filled datetimes (in dataset resolution) are returned in the form
        af a multiindex pandas Series (name -- datetime) as index.

        Parameters
        ----------
        obsdf : Dataset.df
            The Dataset.df attribute. (Needed to extract trailing/leading
                                       observations.)
        outliersdf : Dataset.outliersdf
            The Dataset.outliersdf attribute.(Needed to extract trailing/leading
                                              observations.))
        resolutionseries : Datetime.timedelta
            Resolution of the station observations in the dataset.
        obstype : String, optional
            The observational type to apply gapfilling on. The default is 'temp'.
        method : String, optional
            Method to pass to the Numpy.interpolate function. The default is 'time'.
        max_consec_fill : Integer, optional
            Value to pass to the limit argument of Numpy.interpolate. The default is 100.

        Returns
        -------
        Pandas.Series
            Multiindex Series with filled gap values in dataset space.

        """
        print(f' interpolate on {self}')
        outliersdf = format_outliersdf_to_doubleidx(outliersdf)

        gapfill_series= interpolate_gap(
            gap=self,
            obsdf=obsdf,
            outliersdf=outliersdf,
            dataset_res=dataset_res,
            obstype=obstype,
            method=method,
            max_consec_fill=max_consec_fill,
        )

        # update self
        self.gapfill_technique = 'interpolation'
        self.gapfill_df[obstype] = gapfill_series





# =============================================================================
# Find gaps and missing values
# =============================================================================
def get_station_gaps(gapslist, name):
       """
       Extract a Gap_collection specific to one station. If no gaps are found
       for the station, an empty Gap_collection is returned.

       Parameters
       ----------
       name : String
           Name of the station to extract a Gaps_collection from.

       Returns
       -------
       Gap_collection
           A Gap collection specific of the specified station.

       """
       return [gap for gap in gapslist if gap.name == name]


def get_gaps_indx_in_obs_space(gapslist, obsdf, outliersdf, resolutionseries):
    """

    Explode the gaps, to the dataset resolution and format to a multiindex
    with name -- datetime.

    In addition the last observation before the gap (leading), and first
    observation (after) the gap are computed and stored in the df attribute.
    (the outliers are used to look for leading and trailing observations.)


    Parameters
    ----------
    obsdf : TYPE
        DESCRIPTION.
    outliersdf : TYPE
        DESCRIPTION.
    resolutionseries : TYPE
        DESCRIPTION.

    Returns
    -------
    expanded_gabsidx_obsspace : TYPE
        DESCRIPTION.

    """
    outliersdf = format_outliersdf_to_doubleidx(outliersdf)

    expanded_gabsidx_obsspace = init_multiindex()


    for gap in gapslist:
        gap.update_gaps_indx_in_obs_space(
            obsdf, outliersdf, resolutionseries.loc[gap.name]
        )
        expanded_gabsidx_obsspace = expanded_gabsidx_obsspace.append(
            gap.exp_gap_idx
        )


    return expanded_gabsidx_obsspace


def gaps_to_df(gapslist):
    """
    Combine all gaps into a dataframe as an overview.

    Parameters
    ----------
    gapslist : list
        List of gaps.

    Returns
    -------
    pandas.DataFrame
        A DataFrame with stationnames as index, and the start, end and duretion
        of the gaps as columns.

    """


    gapdflist = []
    for gap in gapslist:
        gapdflist.append(gap.to_df())

    return pd.concat(gapdflist)

def remove_gaps_from_obs(gaplist, obsdf):
    """
    Remove station - datetime records that are in the gaps from the obsdf.

    (Usefull when filling timestamps to a df, and if you whant to remove the
      gaps.)

    Parameters
    ----------
    obsdf : pandas.DataFrame()
        A MultiIndex dataframe with name -- datetime as index.

    Returns
    -------
    obsdf : pandas.DataFrame()
        The same dataframe with records inside gaps removed.

    """

    # Create index for gaps records in the obsdf
    expanded_gabsidx = init_multiindex()
    for gap in gaplist:
        sta_records = obsdf.xs(gap.name, level="name").index  # filter by name

        gaps_dt = sta_records[
            (sta_records >= gap.startgap)
            & (  # filter if the observations are within a gap
                sta_records <= gap.endgap
            )
        ]

        gaps_multiidx = pd.MultiIndex.from_arrays(
            arrays=[[gap.name] * len(gaps_dt), gaps_dt], names=["name", "datetime"]
        )

        expanded_gabsidx = expanded_gabsidx.append(gaps_multiidx)

    # remove gaps idx from the obsdf
    obsdf = obsdf.drop(index=expanded_gabsidx)
    return obsdf

# =============================================================================
# Helpers
# =============================================================================


# def _find_closes_occuring_date(refdt, series_of_dt, where="before"):
#     if where == "before":
#         diff = refdt - (series_of_dt[series_of_dt < refdt])
#     elif where == "after":
#         diff = (series_of_dt[series_of_dt > refdt]) - refdt

#     if diff.empty:
#         # no occurences before of after

#         return np.nan
#     else:
#         return min(diff).total_seconds()



def apply_debias_era5_gapfill(
        gapslist, dataset, eraModelData, debias_settings, obstype="temp",
    ):

        gapfill_settings = dataset.settings.gap['gaps_fill_info']
        expanded_gabsidx_obsspace = init_multiindex()

        filled_gaps_series = pd.Series(
            data=[], index=expanded_gabsidx_obsspace, dtype=object
        )

        # Convert modeldata to the same timzone as the data
        targettz = dataset.df.index.get_level_values('datetime').tz.zone
        eraModelData._conv_to_timezone(targettz)


        for gap in gapslist:
            print(f' Era5 gapfill for {gap}')
            gap.gapfill_technique = "debias gapfill"

            # avoid passing full dataset around
            station = dataset.get_station(gap.name)

            # Update gap attributes
            gap.update_gaps_indx_in_obs_space(
                obsdf=station.df,
                outliersdf=station.outliersdf,
                dataset_res=station.metadf["dataset_resolution"].squeeze(),
            )

            # get leading and trailing period
            leading_obs, trailing_obs = create_leading_trailing_debias_periods(
                station=station,
                gap=gap,
                debias_period_settings=debias_settings["debias_period"],
                obstype=obstype,
            )

            # check if leading/trailing is valid
            if leading_obs.empty | trailing_obs.empty:
                print(
                    "No suitable leading or trailing period found. Gapfill not possible"
                )
                gap.gapfill_errormessage[obstype] = 'gapfill not possible: no leading/trailing period'

                default_return = pd.Series(
                    index=gap.exp_gap_idx, name=obstype, dtype="object"
                )

                default_return.name = obstype
                gapfill_df =default_return.to_frame()
                gapfill_df[obstype + "_" + gapfill_settings["label_columnname"]] = gapfill_settings["label"]["model_debias"]

                # update the gaps attributes
                gap.gapfill_df = gapfill_df

                continue

            # extract model values at leading and trailing period
            leading_model = eraModelData.interpolate_modeldata(leading_obs.index)
            trailing_model = eraModelData.interpolate_modeldata(trailing_obs.index)

            # TODO check if there is modeldata for the leading and trailing + obs period
            if (leading_model[obstype].isnull().any()) | (
                trailing_model[obstype].isnull().any()
            ):
                print(
                    "No modeldata for the full leading/trailing period found. Gapfill not possible"
                )
                gap.gapfill_errormessage[obstype] = 'gapfill not possible: not enough modeldata'

                default_return = pd.Series(
                    index=gap.exp_gap_idx, name=obstype, dtype="object"
                )
                default_return.name = obstype
                gapfill_df =default_return.to_frame()
                gapfill_df[obstype + "_" + gapfill_settings["label_columnname"]] = gapfill_settings["label"]["model_debias"]

                # update the gaps attributes
                gap.gapfill_df = gapfill_df
                continue

            # Get model data for gap timestamps
            gap_model = eraModelData.interpolate_modeldata(gap.exp_gap_idx)

            # apply bias correction
            filled_gap_series, fill_info, err_message = make_era_bias_correction(
                leading_model=leading_model,
                trailing_model=trailing_model,
                gap_model=gap_model,
                leading_obs=leading_obs,
                trailing_obs=trailing_obs,
                obstype=obstype,
            )


            filled_gap_series.name = obstype
            gapfill_df =filled_gap_series.to_frame()
            gapfill_df[obstype + "_" + gapfill_settings["label_columnname"]] = gapfill_settings["label"]["model_debias"]

            # update the gaps attributes
            gap.gapfill_df = gapfill_df
            gap.gapfill_technique = gapfill_settings["label"]["model_debias"]
            gap.gapfill_info = fill_info
            if bool(err_message):
                gap.gapfill_errormessage = err_message





def apply_interpolate_gaps(gapslist, obsdf, outliersdf, dataset_res, gapfill_settings,
                           obstype="temp", method="time", max_consec_fill=100,
                           ):

    """ No return, only update the gaps instances attributes"""

    for gap in gapslist:
        gapfill_series = interpolate_gap(
                        gap=gap,
                        obsdf=obsdf.xs(gap.name, level='name', drop_level=False),
                        outliersdf=outliersdf.xs(gap.name, level='name', drop_level=False),
                        dataset_res=dataset_res.loc[gap.name],
                        obstype=obstype,
                        method=method,
                        max_consec_fill=max_consec_fill,
                        )

        gapfill_series.name = obstype
        gapfill_df = gapfill_series.to_frame()
        gapfill_df[obstype + "_" + gapfill_settings["label_columnname"]] = gapfill_settings["label"]["linear"]

        # update the gaps attributes
        gap.gapfill_df = gapfill_df
        gap.gapfill_technique = gapfill_settings["label"]["linear"]


def make_gapfill_df(gapslist):
    concatlist = []
    for gap in gapslist:
        subgapfill = gap.gapfill_df.reset_index()
        subgapfill['name'] = gap.name
        subgapfill = subgapfill.set_index(['name', 'datetime'])

        concatlist.append(subgapfill)

    return pd.concat(concatlist).sort_index()


def missing_timestamp_and_gap_check(df, gapsize_n):
    # TODO update docstring
    """

    V3
    Looking for missing timestaps by assuming an observation frequency. The assumed frequency is the highest occuring frequency PER STATION.
    If missing observations are detected, they can be catogirized as a missing timestamp or as gap.

    A gap is define as a sequence of missing values with more than N repetitive missing values. N is define in the QC settings.



    Parameters
    ----------
    df : pandas.DataFrame
        The observations dataframe of the dataset object (Dataset.df)

    Returns
    -------

    df : pandas.DataFrame()
        The observations dataframe.
    outlier_df : pandas.DataFrame()
        The dataframe containing the missing timestamps (not gaps) with the outlier label.
    gap_df : pandas.Dataframe()
        The dataframe containing the start and end date of a specific gap.

    """

    gap_list = []
    # gap_df = pd.DataFrame()
    # gap_indices = []
    missing_timestamp_series = pd.Series(dtype=object)
    station_freqs = {}

    # missing timestamp per station (because some stations can have other frequencies!)

    stationnames = df.index.get_level_values(level="name").unique()
    for station in stationnames:
        # find missing timestamps
        timestamps = df.xs(station, level="name").index
        likely_freq = get_likely_frequency(timestamps, method="highest", simplify=False)


        assert likely_freq.seconds > 0, f"The frequency is not positive!"

        station_freqs[station] = likely_freq

        missing_datetimeseries = (pd.date_range(start=timestamps.min(),
                                               end=timestamps.max(),
                                               freq=likely_freq)
                                  .difference(timestamps)
                                  .to_series()
                                  .diff())


        if missing_datetimeseries.empty:
            continue

        # Check for gaps
        gap_defenition = ((missing_datetimeseries != likely_freq)).cumsum()
        consec_missing_groups = missing_datetimeseries.groupby(gap_defenition)
        group_sizes = consec_missing_groups.size()

        gap_groups = group_sizes[group_sizes > gapsize_n]

        # iterate over the gabs and fill the gapsdf
        for gap_idx in gap_groups.index:
            datetime_of_gap_records = consec_missing_groups.get_group(gap_idx).index
            gap = Gap(name=station,
                      startdt=datetime_of_gap_records.min(),
                      enddt=datetime_of_gap_records.max())
            gap_list.append(gap)



        # combine the missing timestams values
        missing_timestamp_groups = group_sizes[group_sizes <= gapsize_n]
        for missing_idx in missing_timestamp_groups.index:
            datetime_of_missing_records = consec_missing_groups.get_group(
                missing_idx
            ).index.to_list()

            missing_timestamp_series = pd.concat(
                [
                    missing_timestamp_series,
                    pd.Series(
                        index=[station] * len(datetime_of_missing_records),
                        data=datetime_of_missing_records,
                    ),
                ]
            )

    missing_obs_collection = Missingob_collection(missing_timestamp_series)
    df = df.sort_index()

    return missing_obs_collection, gap_list


