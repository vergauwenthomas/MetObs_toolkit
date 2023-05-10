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
)



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

        # computed attributes
        self.leading_timestamp = None  # last ob_dt before gap in datset space
        self.trailing_timestamp = None  # first ob_dt after gap in dataset space

        self.exp_gap_idx = None

        # gap fill (only for conventional saving)
        self.gapfill_values = None
        self.gapfill_technique = None

    def __str__(self):
        return f"Gap instance of {self.name} for {self.startgap} --> {self.endgap}"
    def __repr__(self):
        return self.__str__()


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
            index=[self.name], data={"start_gap": self.startgap, "end_gap": self.endgap}
        )

    def update_leading_trailing_obs(self, obsdf, outliersdf):
        """
        Add the leading (last obs before gap) and trailing (first obs after gap)
        as extra columns to the self.df.

        The obsdf and outliersdf are both used to scan for the leading and trailing obs.

        Parameters
        ----------
        obsdf : pandas.DataFrame
            Dataset.df
        outliersdf : pandas.DataFrame
            Dataset.outliersdf

        Returns
        -------
        None.

        """
        outliersdf = format_outliersdf_to_doubleidx(outliersdf)

        # combine timestamps of observations and outliers
        sta_obs = obsdf.xs(self.name, level="name").index
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

        outliersdf = format_outliersdf_to_doubleidx(outliersdf)

        gapfill_series = interpolate_gap(
            gap=self,
            obsdf=obsdf,
            outliersdf=outliersdf,
            dataset_res=dataset_res,
            obstype=obstype,
            method=method,
            max_consec_fill=max_consec_fill,
        )
        gapdf = gapfill_series.to_frame().reset_index()
        gapdf["name"] = self.name
        gapdf.index = pd.MultiIndex.from_arrays(
            arrays=[gapdf["name"].values, gapdf["datetime"].values],
            names=["name", "datetime"],
        )
        return gapdf[obstype]

    def get_leading_trailing_debias_periods(self, station, obstype, debias_periods):
        # get debias periods

        leading_period, trailing_period = create_leading_trailing_debias_periods(
            station=station,
            gap=self,
            debias_period_settings=debias_periods,
            obstype=obstype,
        )

        return leading_period, trailing_period


class Gap_collection:
    def __init__(self, gapsdf):
        self.list = [
            Gap(sta, row["start_gap"], row["end_gap"]) for sta, row in gapsdf.iterrows()
        ]

    def __str__(self):
        if not bool(self.list):
            return f'Empty gap collection'
        longstring = ''
        for gap in self.list:
            longstring += str(gap) + '\n'
        return f"Gap collection for: \n {longstring}"
    def __repr__(self):
        return self.__str__()

    def to_df(self):
        gaps_names = []
        gaps_startdt = []
        gaps_enddt = []
        for gap in self.list:
            gaps_names.append(gap.name)
            gaps_startdt.append(gap.startgap)
            gaps_enddt.append(gap.endgap)

        df = pd.DataFrame(
            index=pd.Index(gaps_names),
            data={"start_gap": gaps_startdt, "end_gap": gaps_enddt},
        )
        df.index.name = "name"

        return df

    def get_station_gaps(self, name):
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
        gapdf = self.to_df()

        if name in gapdf.index:
            return Gap_collection(gapdf.loc[name])
        else:
            return Gap_collection(pd.DataFrame())

    def remove_gaps_from_obs(self, obsdf):
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
        expanded_gabsidx = pd.MultiIndex(
            levels=[["name"], ["datetime"]], codes=[[], []], names=["name", "datetime"]
        )

        for gap in self.list:
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

    def get_gaps_indx_in_obs_space(self, obsdf, outliersdf, resolutionseries):
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

        expanded_gabsidx_obsspace = pd.MultiIndex(
            levels=[["name"], ["datetime"]], codes=[[], []], names=["name", "datetime"]
        )

        for gap in self.list:
            gap.update_gaps_indx_in_obs_space(
                obsdf, outliersdf, resolutionseries.loc[gap.name]
            )
            expanded_gabsidx_obsspace = expanded_gabsidx_obsspace.append(
                gap.exp_gap_idx
            )

        return expanded_gabsidx_obsspace

    def apply_interpolate_gaps(
        self,
        obsdf,
        outliersdf,
        dataset_res,
        obstype="temp",
        method="time",
        max_consec_fill=100,
    ):
        outliersdf = format_outliersdf_to_doubleidx(outliersdf)

        expanded_gabsidx_obsspace = pd.MultiIndex(
            levels=[["name"], ["datetime"]], codes=[[], []], names=["name", "datetime"]
        )
        filled_gaps_series = pd.Series(
            data=[], index=expanded_gabsidx_obsspace, dtype=object
        )

        for gap in self.list:
            gapfill_series = interpolate_gap(
                gap=gap,
                obsdf=obsdf,
                outliersdf=outliersdf,
                dataset_res=dataset_res.loc[gap.name],
                obstype=obstype,
                method=method,
                max_consec_fill=max_consec_fill,
            )

            gapdf = gapfill_series.to_frame().reset_index()
            gapdf["name"] = gap.name
            gapdf = gapdf.set_index(['name', 'datetime'])
            # gapdf.index = pd.MultiIndex.from_arrays(
            #     arrays=[gapdf["name"].values, gapdf["datetime"].values],
            #     names=["name", "datetime"],
            # )

            # Update gap
            gap.gapfill_technique = "interpolation"
            gap.gapfill_values = gapdf[obstype]

            filled_gaps_series = pd.concat([filled_gaps_series, gapdf[obstype]])
        return filled_gaps_series

    def apply_debias_era5_gapfill(
        self, dataset, eraModelData, debias_settings, obstype="temp"
    ):
        expanded_gabsidx_obsspace = pd.MultiIndex(
            levels=[["name"], ["datetime"]], codes=[[], []], names=["name", "datetime"]
        )
        filled_gaps_series = pd.Series(
            data=[], index=expanded_gabsidx_obsspace, dtype=object
        )

        for gap in self.list:
            # avoid passing full dataset around
            station = dataset.get_station(gap.name)

            # Update gap attributes
            gap.update_gaps_indx_in_obs_space(
                obsdf=station.df,
                outliersdf=station.outliersdf,
                dataset_res=station.metadf["dataset_resolution"].squeeze(),
            )

            # get leading and trailing period
            leading_obs, trailing_obs = gap.get_leading_trailing_debias_periods(
                obstype=obstype,
                station=dataset.get_station(gap.name),
                debias_periods=debias_settings["debias_period"],
            )
            # check if leading/trailing is valid
            if leading_obs.empty | trailing_obs.empty:
                print(
                    "No suitable leading or trailing period found. Gapfill not possible"
                )
                gap.gapfill_technique = (
                    "debias era5 gapfill (not possible: no leading/trailing period)"
                )
                default_return = pd.Series(
                    index=gap.exp_gap_idx, name=obstype, dtype="object"
                )
                gap.gapfill_values = default_return
                filled_gaps_series = pd.concat([filled_gaps_series, default_return])
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
                gap.gapfill_technique = (
                    "debias era5 gapfill (not possible: not enough modeldata)"
                )
                default_return = pd.Series(
                    index=gap.exp_gap_idx, name=obstype, dtype="object"
                )
                gap.gapfill_values = default_return
                filled_gaps_series = pd.concat([filled_gaps_series, default_return])

            # Get model data for gap timestamps
            gap_model = eraModelData.interpolate_modeldata(gap.exp_gap_idx)

            # apply bias correction
            filled_gap_series = make_era_bias_correction(
                leading_model=leading_model,
                trailing_model=trailing_model,
                gap_model=gap_model,
                leading_obs=leading_obs,
                trailing_obs=trailing_obs,
                obstype=obstype,
            )

            # Update gap
            gap.gapfill_technique = "debias era5 gapfill"
            gap.gapfill_values = filled_gap_series

            filled_gaps_series = pd.concat([filled_gaps_series, filled_gap_series])

        filled_gaps_series.name = obstype
        return filled_gaps_series


# =============================================================================
# Find gaps and missing values
# =============================================================================


def _find_closes_occuring_date(refdt, series_of_dt, where="before"):
    if where == "before":
        diff = refdt - (series_of_dt[series_of_dt < refdt])
    elif where == "after":
        diff = (series_of_dt[series_of_dt > refdt]) - refdt

    if diff.empty:
        # no occurences before of after

        return np.nan
    else:
        return min(diff).total_seconds()


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

    gap_df = pd.DataFrame()
    gap_indices = []
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
            # fill the gaps df
            datetime_of_gap_records = consec_missing_groups.get_group(gap_idx).index
            gap_df = pd.concat(
                [
                    gap_df,
                    pd.DataFrame(
                        data=[
                            [
                                datetime_of_gap_records.min(),
                                datetime_of_gap_records.max(),
                            ]
                        ],
                        index=[station],
                        columns=["start_gap", "end_gap"],
                    ),
                ]
            )

            logger.debug(
                f"Data gap from {datetime_of_gap_records.min()} --> {datetime_of_gap_records.max()} found for {station}."
            )
            gap_indices.extend(
                list(
                    zip(
                        [station] * datetime_of_gap_records.shape[0],
                        datetime_of_gap_records,
                    )
                )
            )

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

    df = df.sort_index()

    return df, missing_timestamp_series, gap_df
