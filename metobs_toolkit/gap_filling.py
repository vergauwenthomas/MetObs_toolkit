#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 17:05:26 2023

@author: thoverga
"""
import numpy as np
import pandas as pd
from datetime import timedelta

from metobs_toolkit.df_helpers import (
    remove_outliers_from_obs,
    init_multiindexdf,
    format_outliersdf_to_doubleidx,
)

# def interpolate_missing_obs(missing_ob, obsdf, outliersdf, obstype, method):
#     # 1 get trailing and leading obs



# =============================================================================
# Gap fillers
# =============================================================================


def interpolate_gap(
    gap, obsdf, outliersdf, dataset_res, obstype, method, max_consec_fill
):
    outliersdf = format_outliersdf_to_doubleidx(outliersdf)

    # 1 get trailing and leading + exploded index
    gap.update_gaps_indx_in_obs_space(obsdf, outliersdf, dataset_res)
    gap.update_leading_trailing_obs(obsdf, outliersdf, obs_only=True)

    # initiate return value when no interpolation can be performed
    empty_interp = pd.Series(data=np.nan, index=gap.exp_gap_idx.droplevel("name"))
    empty_interp.name = obstype

    # 2 check if there is a trailing and leading gap
    if gap.startgap == gap.leading_timestamp:
        message =f"No leading timestamp found for gap {gap}"
        print(message)
        gap.gapfill_errormessage[obstype]=message
        return empty_interp

    if gap.endgap == gap.trailing_timestamp:
        message = f"No trailing timestamp found for gap {gap}"
        print(message)
        gap.gapfill_errormessage[obstype]=message
        return empty_interp


    # 3. Get leading and trailing val
    if not bool(gap.leading_val):
        # empty dict --> no value in the obs
        message =f"No cadidate for leading {obstype} observation found for {gap}"
        print(message)
        gap.gapfill_errormessage[obstype]=message
        return empty_interp

    if not bool(gap.trailing_val):
        # empty dict --> no value in the obs
        message =f"No cadidate for trailing {obstype} observation found for {gap}"
        print(message)
        gap.gapfill_errormessage[obstype]=message
        return empty_interp


    leading_dt = gap.leading_timestamp
    leading_val = gap.leading_val[obstype]
    trailing_dt = gap.trailing_timestamp
    trailing_val = gap.trailing_val[obstype]


    # # 3 check both leading and trailing are in obs, and look for alternative leading/trailing if the original is an outlier.
    # sta_obs = obsdf.xs(gap.name, level="name")

    # # leading
    # if gap.leading_timestamp in sta_obs.index:
    #     # leading found in obs
    #     leading_dt = gap.leading_timestamp
    #     leading_val = sta_obs.loc[gap.leading_timestamp, obstype]
    # else:
    #     # look for last observation before leading timestamp
    #     delta_dt = (
    #         gap.leading_timestamp - sta_obs.index[sta_obs.index < gap.leading_timestamp]
    #     )
    #     if delta_dt.empty:
    #         print(
    #             f"No cadidate for leading {obstype} observation found for {gap.name} with gap: {gap.startgap} --> {gap.endgap}"
    #         )
    #         return empty_interp

    #     leading_dt = gap.leading_timestamp - delta_dt.min()
    #     leading_val = sta_obs.loc[leading_dt, obstype]

    # # trailing
    # if gap.trailing_timestamp in sta_obs.index:
    #     # leading found in obs
    #     trailing_dt = gap.trailing_timestamp
    #     trailing_val = sta_obs.loc[gap.trailing_timestamp, obstype]
    # else:
    #     # look for last observation before leading timestamp
    #     delta_dt = (
    #         sta_obs.index[sta_obs.index > gap.trailing_timestamp]
    #         - gap.trailing_timestamp
    #     )
    #     if delta_dt.empty:
    #         print(
    #             f"No cadidate for trailing {obstype} observation found for {gap.name} with gap: {gap.startgap} --> {gap.endgap}"
    #         )
    #         return empty_interp
    #     # TODO: settings restrictions on maximum delta_dt ??
    #     trailing_dt = gap.trailing_timestamp + delta_dt.min()
    #     trailing_val = sta_obs.loc[trailing_dt, obstype]

    # Make interpolation series
    gaps_series = pd.Series(data=np.nan, index=gap.exp_gap_idx.droplevel("name"))
    gaps_series = pd.concat(
        [
            gaps_series,
            pd.Series(
                index=[leading_dt, trailing_dt], data=[leading_val, trailing_val]
            ),
        ]
    )
    gaps_series = gaps_series.sort_index()

    # Interpolate series
    gaps_series.interpolate(
        method=method,
        limit=max_consec_fill,  # Maximum number of consecutive NaNs to fill. Must be greater than 0.
        limit_area="inside",
        inplace=True,
    )

    # Subset only gap indixes
    gaps_fill_series = gaps_series[gap.exp_gap_idx.droplevel("name")]
    gaps_fill_series.name = obstype

    # update gapfill info
    gap.gapfill_info = gaps_series.to_frame()


    return gaps_fill_series


# =============================================================================
# Debiasing period
# =============================================================================


def get_sample_size(sample_duration_hours, sta):
    stares = sta.metadf["dataset_resolution"].squeeze()
    sample_size = timedelta(hours=sample_duration_hours) / stares
    return int(sample_size)


def create_leading_trailing_debias_periods(
    station, gap, debias_period_settings, obstype
):
    # Get samplesizes
    debias_pref_sample_size_leading = get_sample_size(
        debias_period_settings["prefered_leading_sample_duration_hours"], station
    )
    debias_pref_sample_size_trailing = get_sample_size(
        debias_period_settings["prefered_trailing_sample_duration_hours"], station
    )
    debias_min_sample_size_leading = get_sample_size(
        debias_period_settings["minimum_leading_sample_duration_hours"], station
    )
    debias_min_sample_size_trailing = get_sample_size(
        debias_period_settings["minimum_trailing_sample_duration_hours"], station
    )

    # get all observations that can be used for debias training
    obs = station.df

    # remove blacklist
    # TODO
    obs = remove_outliers_from_obs(
        obs, format_outliersdf_to_doubleidx(station.outliersdf)
    )

    # add whitelist
    # TODO

    # only datetimes are relevant
    obs = obs.reset_index()
    obs = obs[["name", "datetime", obstype]]

    # Select all leading and all trailing obs
    leading_period = obs[obs["datetime"] < gap.startgap]
    trailing_period = obs[obs["datetime"] > gap.endgap]
    print(f'   {leading_period.shape[0]} leading records, {trailing_period.shape[0]} trailing records.')

    # some derived integers
    poss_shrinkage_leading = leading_period.shape[0] - debias_min_sample_size_leading
    poss_shrinkage_trailing = trailing_period.shape[0] - debias_min_sample_size_trailing
    poss_extention_leading = leading_period.shape[0] - debias_pref_sample_size_leading
    poss_extention_trailing = (
        trailing_period.shape[0] - debias_pref_sample_size_trailing
    )

    # check if desired sample sizes for leading and trailing are possible
    if (leading_period.shape[0] >= debias_pref_sample_size_leading) & (
        trailing_period.shape[0] >= debias_pref_sample_size_trailing
    ):
        print(f"leading and trailing periods are both available for debiassing.")
        # both periods are oke
        leading_df = leading_period[-debias_pref_sample_size_leading:]
        trailing_df = trailing_period[:debias_pref_sample_size_trailing]

    elif (leading_period.shape[0] <= debias_pref_sample_size_leading) & (
        trailing_period.shape[0] >= debias_pref_sample_size_trailing
    ):
        print(
            f"leading periods for debiassing does not have a preferable size. Try translation/shrinkage ..."
        )

        # leading period to small, trailing period is OK

        missing_records = debias_pref_sample_size_leading - leading_period.shape[0]

        # 1 if the leading period is smaller thatn the minimum leading size --> return default
        if poss_shrinkage_leading < 0:
            leading_df = init_multiindexdf()
            trailing_df = init_multiindexdf()  # this might be to strict
            print(
                "The available leading debias samplesize is smaller than the minimum. A translation/shrinking is not possible."
            )

        # 2 Try translation without shrinkage

        elif missing_records <= poss_extention_trailing:
            # translation without shrinkage is possible
            translation_trailing = missing_records

            leading_df = leading_period
            trailing_df = trailing_period[
                0 : (debias_pref_sample_size_trailing + translation_trailing)
            ]

            print(
                f"A translation of {translation_trailing} records is done towards the trailing period. (n_leading + n_trailing is conserved: {leading_df.shape[0] + trailing_df.shape[0]}"
            )

        # 3. Try if a translation is within the limits of shrinkage
        elif (missing_records - poss_extention_trailing) <= poss_shrinkage_leading:
            translation_trailing = poss_extention_trailing

            leading_df = leading_period
            trailing_df = trailing_period[
                0 : debias_pref_sample_size_trailing + translation_trailing
            ]
            print(
                f"A translation of {translation_trailing} records is done towards the trailing period. Since there was not engough translation space for the trailing obs, the condition n_leading + n_trailing is NOT conserved: {leading_df.shape[0] + trailing_df.shape[0]}. \
                  Both leading and trailing sizes still achieves minimal size restrictions."
            )
        # 4. If all else fails, it is not possible to make a leading period
        else:
            print(
                "The available leading samplesize can not reach minimal size restrictions."
            )
            # no translation is possible, even with shrinking
            leading_df = init_multiindexdf()
            trailing_df = init_multiindexdf()  # this might be to strict

    elif (leading_period.shape[0] >= debias_pref_sample_size_leading) & (
        trailing_period.shape[0] <= debias_pref_sample_size_trailing
    ):
        # leading period is ok, trailing period is to short
        print(
            f"trailing periods for debiassing does not have a preferable size. Try translation/shrinkage ..."
        )
        missing_records = debias_pref_sample_size_trailing - trailing_period.shape[0]

        # 1 if the trailing period is smaller thatn the minimum trailing size --> return default
        if poss_shrinkage_trailing < 0:
            leading_df = init_multiindexdf()  # might be to strict
            trailing_df = init_multiindexdf()
            print(
                "The available trailing debias samplesize is smaller than the minimum. A translation/shrinking is not possible."
            )
            # return

        # 2 Try translation without shrinkage
        elif missing_records <= poss_extention_leading:
            # translation without shrinkage is possible
            translation_leading = missing_records

            leading_df = leading_period[
                -(debias_pref_sample_size_leading + translation_leading) :
            ]
            trailing_df = trailing_period
            print(
                f"A translation of {translation_leading} records is done towards the leading period. (n_leading + n_trailing is conserved: {leading_df.shape[0] + trailing_df.shape[0]}"
            )

        # 3. Try if a translation is within the limits of shrinkage
        elif (missing_records - poss_extention_leading) <= poss_shrinkage_trailing:
            translation_leading = poss_extention_leading

            leading_df = leading_period[
                -(debias_pref_sample_size_leading + translation_leading)
            ]
            trailing_df = trailing_period
            print(
                f"A translation of {translation_leading} records is done towards the leading period. Since there was not engough translation space for the leading obs, the condition n_leading + n_trailing is NOT conserved: {leading_df.shape[0] + trailing_df.shape[0]}. \
                  Both leading and trailing sizes still achieves minimal size restrictions."
            )
        # 4. If all else fails, it is not possible to make a trailing period
        else:
            # no translation is possible, even with shrinking
            print(
                "The available trailing samplesize can not reach minimal size restrictions."
            )
            leading_df = init_multiindexdf()  # this might be to strict
            trailing_df = init_multiindexdf()

    else:
        # Both leading and trailing periods are not to small

        # 1 does both (leading and trailing) still acchieves the minimal size condition for shrinking?
        if (poss_shrinkage_leading >= 0) & (poss_shrinkage_trailing >= 0):
            print(
                "Both leading and trailing periods do not have a prefered size, but still meet the minimal conditions."
            )
            leading_df = leading_period
            trailing_df = trailing_period

        else:
            print(
                "Both leading and trailing periods do not have a prefered size, and eighter of them does NOT meet minimal condition."
            )
            # either one of the periods does not reach minimal condition, so return default
            leading_df = init_multiindexdf()
            trailing_df = init_multiindexdf()

    # convert to multiindex
    if not leading_df.empty:
        leading_df = leading_df.set_index(["name", "datetime"])
    if not trailing_df.empty:
        trailing_df = trailing_df.set_index(["name", "datetime"])

    return leading_df, trailing_df


def get_time_specific_biases(model, obs, obstype, period):
    diff = model - obs
    diff = diff.reset_index().set_index("datetime")
    diff["hours"] = diff.index.hour
    diff["minutes"] = diff.index.minute
    diff["seconds"] = diff.index.second

    biases = diff.groupby(["name", "hours", "minutes", "seconds"])[obstype].mean()
    biases.name = obstype + "_bias_" + period

    biases = biases.reset_index()
    return biases


def make_era_bias_correction(leading_model, trailing_model,
                             gap_model, leading_obs, trailing_obs,
                             obstype):
    error_message = ''
    # 1. get lead timestamp biases
    lead_biases = get_time_specific_biases(
        model=leading_model, obs=leading_obs, obstype=obstype, period="lead"
    )

    # 2. get trailing timestamp biases
    trail_biases = get_time_specific_biases(
        model=trailing_model, obs=trailing_obs, obstype=obstype, period="trail"
    )

    # 3. apply bias correction on modeldata in gap

    # linear interpolation of bias along the gap method:
    gap_model["trail_weight"] = np.linspace(0.0, 1.0, gap_model.shape[0])
    gap_model["lead_weight"] = 1.0 - gap_model["trail_weight"]

    # aggregate to timestamps
    gap_model["hours"] = gap_model.index.get_level_values("datetime").hour
    gap_model["minutes"] = gap_model.index.get_level_values("datetime").minute
    gap_model["seconds"] = gap_model.index.get_level_values("datetime").second


    # testing
    gap_model = gap_model.reset_index()

    gap_model = gap_model.merge(
        right=lead_biases[["hours", "minutes", "seconds", obstype + "_bias_lead"]],
        how="left",
        on=["hours", "minutes", "seconds"],
    )

    gap_model = gap_model.merge(
         right=trail_biases[["hours", "minutes", "seconds", obstype + "_bias_trail"]],
         how="left",
         on=["hours", "minutes", "seconds"],
     )

    gap_model = gap_model.set_index(['name', 'datetime'])

    # Idea: if BOTH leadin and trailing (hourly) biases is available, than use
    # use the debias corection (even if it is for a part of the gap!).
    # If either one or both are missing, than no bias correction is applied
    no_debias = gap_model[(gap_model[obstype + '_bias_lead'].isnull()) |
                          (gap_model[obstype + '_bias_trail'].isnull())].index
    error_message =f'WARNING!, No debias possible for these gap records: {no_debias},the gap will be filled by model data without bias correction. '
    print(error_message)


    # set weights to zero if not debias correction can be applied on that record
    gap_model.loc[no_debias, obstype+'_bias_trail'] = 0.
    gap_model.loc[no_debias, obstype+'_bias_lead'] = 0.


    # 5. compute the debiased fill value
    # leave this dataframe for debugging
    gap_model[obstype + "_fill"] = gap_model[obstype] - (
        (gap_model["lead_weight"] * gap_model[obstype + "_bias_lead"])
        + (gap_model["trail_weight"] * gap_model[obstype + "_bias_trail"])
    )

    # 6. make returen
    returnseries = gap_model[obstype + "_fill"]
    returnseries.name = obstype
    return returnseries, gap_model, error_message
