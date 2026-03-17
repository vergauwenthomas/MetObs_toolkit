from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Tuple
import numpy as np
import pandas as pd

logger = logging.getLogger("<metobs_toolkit>")


if TYPE_CHECKING:
    from ..buddywrapsensor import BuddyWrapSensor

# Import constants from buddywrapsensor
from ..buddywrapsensor import BC_NO_BUDDIES, BC_PASSED, BC_FLAGGED, BC_NOT_TESTED


def buddy_test_a_station(
    centerwrapsensor: BuddyWrapSensor,
    buddygroupname: str,
    widedf: pd.DataFrame,
    min_sample_size: int,
    min_sample_spread: float,
    outlier_threshold: float,
    iteration: int,
    check_type: str = 'spatial_check',
    use_z_robust_method: bool = True,
) -> Tuple[pd.MultiIndex, BuddyWrapSensor]:

    #TODO update docstring
    """Find outliers in a buddy group and update station flags/details.
    
    This function tests whether the center station is an outlier compared to
    its buddy stations using z-score analysis. The z-score is computed using
    the mean and standard deviation of the buddy stations only (the center
    station's values are excluded from the sample distribution).
    
    Parameters
    ----------
    centerwrapsensor : BuddyWrapSensor
        The wrapped sensor at the center of the buddy group to be tested.
    buddygroupname : str
        The name of the buddy group to use.
    widedf : pd.DataFrame
        Wide-format DataFrame with stations as columns and timestamps as index.
    min_sample_size : int
        Minimum number of valid buddy samples required for z-score calculation.
    min_sample_spread : float
        when use_z_robust_method is True, this is the equal to the minimum MAD to use (avoids division by near-zero).
        when use_z_robust_method is False, this is the standard deviation.
    outlier_threshold : float
        Z-score threshold above which a record is flagged as outlier.
    iteration : int
        The current iteration number.
    check_type : str, optional
        The type of check being performed ('spatial_check', 'safetynet_check:groupname').
        Default is 'spatial_check'.
        
    Returns
    -------
    pd.MultiIndex
        MultiIndex with levels ('name', 'datetime') containing only the outlier
        records for the center station.
    """
    
    # Get buddies (excluding center station)
    buddies = centerwrapsensor.get_buddies(groupname=buddygroupname)
    center_name = centerwrapsensor.name
    
    # Subset to buddies only (for sample distribution) and center station
    buddydf = widedf[buddies]
    center_series = widedf[center_name]
    
    # Count valid buddy samples per timestamp (center station NOT included)
    buddy_sample_sizes = buddydf.notna().sum(axis=1)
    
    # Mark timestamps where center station has no data as NOT_TESTED
    
    #Edgecaase: If a station has fewer records than others, they are present as NaN in widedf 
    #But these timestamps do not ex
    no_data = pd.Series(BC_NOT_TESTED, index=center_series[center_series.isna()].index)
    centerwrapsensor.add_flags(
        iteration=iteration,
        flag_series=no_data,
        column_name=check_type
    )
    
    
    # Find timestamps where center station has data
    center_has_data = center_series.notna()
    #TODO: pass the flag BC_NOT_TESTED for timestamps where center has no data
    
    # Separate timestamps by sample size condition (only where center has data)
    sufficient_samples_mask = (buddy_sample_sizes >= min_sample_size) & center_has_data
    insufficient_samples_mask = (buddy_sample_sizes < min_sample_size) & center_has_data
    
    timestamps_with_sufficient = widedf.index[sufficient_samples_mask]
    timestamps_insufficient = widedf.index[insufficient_samples_mask]
    
    # ---- Handle timestamps with insufficient buddy samples (BC_NO_BUDDIES) ----
    if not timestamps_insufficient.empty:
        # Create flags for NO_BUDDIES
        no_buddies_flags = pd.Series(BC_NO_BUDDIES, index=timestamps_insufficient)
        centerwrapsensor.add_flags(
            iteration=iteration,
            flag_series=no_buddies_flags,
            column_name=check_type
        )
        
        # Create detail messages
        no_buddies_details = pd.Series(
            [f"Insufficient buddy sample size (n={int(buddy_sample_sizes.loc[ts])}, "
             f"required={min_sample_size}) in {buddygroupname} buddy group "
             f"centered on {center_name}"
             for ts in timestamps_insufficient],
            index=timestamps_insufficient
        )
        if check_type == 'spatial_check':
            centerwrapsensor.add_spatial_details(
                iteration=iteration,
                detail_series=no_buddies_details
            )
        else:
            centerwrapsensor.add_safetynet_details(iteration=iteration,
                                                   safetynetname=buddygroupname,
                                                   detail_series=no_buddies_details)
    
    # ---- Handle timestamps with sufficient samples ----
    if timestamps_with_sufficient.empty:
        # No timestamps to process, return empty MultiIndex
        return (pd.MultiIndex.from_tuples([], names=['name', 'datetime']), centerwrapsensor)
    
    # Filter to rows with enough valid buddy samples
    buddydf_filtered = buddydf.loc[timestamps_with_sufficient]
    center_filtered = center_series.loc[timestamps_with_sufficient]
    buddy_sample_sizes_filtered = buddy_sample_sizes.loc[timestamps_with_sufficient]
    
    # Compute z-scores for center station using buddy distribution
    if use_z_robust_method:
        results_df = _compute_robust_z_scores(
            buddydf=buddydf_filtered,
            center_values=center_filtered,
            min_mad=min_sample_spread,
            outlier_threshold=outlier_threshold
        )
        
    else:
        results_df = _compute_center_z_scores(
            buddydf=buddydf_filtered,
            center_values=center_filtered,
            min_std=min_sample_spread,
            outlier_threshold=outlier_threshold
        )
    
    # Separate flagged (outliers) and passed
    outlier_timestamps = results_df.index[results_df['flagged']]
    passed_timestamps = results_df.index[~results_df['flagged']]
    
    # ---- Update PASSED flags and details ----
    if not passed_timestamps.empty:
        passed_flags = pd.Series(BC_PASSED, index=passed_timestamps)
        centerwrapsensor.add_flags(
            iteration=iteration,
            flag_series=passed_flags,
            column_name=check_type
        )
        
        # Create detail messages for passed
        passed_details = f"Passed {buddygroupname} check: " + results_df.loc[passed_timestamps, 'details']
        
        if check_type == 'spatial_check':
            centerwrapsensor.add_spatial_details(
                iteration=iteration,
                detail_series=passed_details
            )
        else:
            centerwrapsensor.add_safetynet_details(iteration=iteration,
                                                   safetynetname=buddygroupname,
                                                   detail_series=passed_details)
            
    # ---- Update FLAGGED (outlier) flags and details ----
    if not outlier_timestamps.empty:
        flagged_flags = pd.Series(BC_FLAGGED, index=outlier_timestamps)
        centerwrapsensor.add_flags(
            iteration=iteration,
            flag_series=flagged_flags,
            column_name=check_type
        )
        
        # Create detail messages for outliers
        outlier_details = f"Outlier in {buddygroupname} buddy group centered on {center_name}: " + results_df.loc[outlier_timestamps, 'details']
          
        
        
        if check_type == 'spatial_check':
            centerwrapsensor.add_spatial_details(
                iteration=iteration,
                detail_series=outlier_details
            )
        else:
            centerwrapsensor.add_safetynet_details(iteration=iteration,
                                                   safetynetname=buddygroupname,
                                                   detail_series=outlier_details)
    
    # ---- Return outliers as MultiIndex ----
    if not outlier_timestamps.empty:
        outlier_multiindex = pd.MultiIndex.from_arrays(
            [[center_name] * len(outlier_timestamps), outlier_timestamps],
            names=['name', 'datetime']
        )
        #Return the updated stations, this is needed when runned in multiprocessing
        return (outlier_multiindex, centerwrapsensor)
    else:
        return (pd.MultiIndex.from_tuples([], names=['name', 'datetime']), centerwrapsensor)



# ------------------------------------------
#    Statistical sample scoring
# ------------------------------------------

def _compute_robust_z_scores(
    buddydf: pd.DataFrame,
    center_values: pd.Series,
    min_mad: float,
    outlier_threshold: float
    ) -> pd.DataFrame:
    """Compute robust z-scores (MADFM-based) for outlier detection.

    Each centre-station value is compared against the median of its
    buddy values using the Median Absolute Deviation (MAD) as a
    robust measure of spread.

    Parameters
    ----------
    buddydf : pandas.DataFrame
        Wide DataFrame containing only the buddy stations' values
        (centre station excluded).  Rows are timestamps.
    center_values : pandas.Series
        Observations of the centre station, aligned to ``buddydf``'s index.
    min_mad : float
        Minimum MAD value to use in the denominator, preventing division
        by near-zero values.
    outlier_threshold : float
        Robust z-score above which a record is flagged as an outlier.

    Returns
    -------
    pandas.DataFrame
        DataFrame with the same index as ``buddydf`` and three columns:

        ``'z_score'``
            Computed robust z-score for each timestamp.
        ``'flagged'``
            Boolean; True when ``z_score > outlier_threshold``.
        ``'details'``
            Human-readable string summarising the z-score calculation.
    """
    buddy_not_na_counts = buddydf.notna().sum(axis=1)
    #Calculate MADFM (Median Absolute Deviation From Median)
    def MAD(x):
        "MEDIAN absolute deviation from median"
        return (x - x.median()).abs().median()
    
    mad_series = buddydf.apply(MAD, axis=1)
    # Replace std below minimum with the minimum (avoid division by near-zero)
    mad_series.loc[mad_series < min_mad] = np.float32(min_mad)
    
    # Calculate robust z-score for center station
    
    robust_z_scores = (center_values - buddydf.median(axis=1)).abs() / (1.4826 * mad_series)
    
    details = ('z (robust)=' + robust_z_scores.map('{:.2f}'.format) +
               ', threshold=' + str(outlier_threshold) + 
               ', n=' + buddy_not_na_counts.map(str) +
               ', MAD=' + mad_series.map('{:.2f}'.format) +
               ', median=' + buddydf.median(axis=1).map('{:.2f}'.format))
    
    # Build result DataFrame
    result_df = pd.DataFrame(
        index=buddydf.index,
        data={
            'z_score': robust_z_scores,
            'flagged': robust_z_scores > outlier_threshold,
            'details': details,
        }
    )
    return result_df
    

def _compute_center_z_scores(
    buddydf: pd.DataFrame,
    center_values: pd.Series,
    min_std: float,
    outlier_threshold: float
) -> pd.DataFrame:
    """Compute z-scores for center station using buddy distribution.
    
    The z-score is computed as the absolute deviation of the center station's
    value from the mean of the buddy stations, divided by the standard
    deviation of the buddy stations. The center station's values are NOT
    included in the mean/std calculation.
    
    Parameters
    ----------
    buddydf : pd.DataFrame
        DataFrame with buddy stations as columns (center station excluded).
    center_values : pd.Series
        Series with the center station's values to test.
    min_std : float
        Minimum standard deviation to use (avoids division by near-zero).
    outlier_threshold : float
        Z-score threshold above which a record is flagged as outlier.
        
    Returns
    -------
    pd.DataFrame
        DataFrame with columns: 'z_score', 'flagged', 'buddy_mean', 'buddy_std'.
    """
    # Compute mean and std from buddies only (center station excluded)
    buddy_mean_series = buddydf.mean(axis=1)
    buddy_std_series = buddydf.std(axis=1)
    buddy_not_na_counts = buddydf.notna().sum(axis=1)
    # Replace std below minimum with the minimum (avoid division by near-zero)
    buddy_std_series.loc[buddy_std_series < min_std] = np.float32(min_std)
    
    # Calculate z-score for center station
    z_scores = (center_values - buddy_mean_series).abs() / buddy_std_series
    
    
    details = ('z=' + z_scores.map('{:.2f}'.format) +
               ', threshold=' + str(outlier_threshold) + 
               ', n=' + buddy_not_na_counts.map(str) +
               ', mean=' + buddy_mean_series.map('{:.2f}'.format) +
               ', std=' + buddy_std_series.map('{:.2f}'.format))
    
    # Build result DataFrame
    result_df = pd.DataFrame(
        index=buddydf.index,
        data={
            'z_score': z_scores,
            'flagged': z_scores > outlier_threshold,
            'details': details,
        }
    )
    return result_df