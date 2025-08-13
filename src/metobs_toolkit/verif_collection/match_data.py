import pandas as pd


def match_obs_and_model(obsdf, modeldatadf, tolerance=pd.Timedelta('5min')):
    """
    Match observations and model data based on datetime, obstype, and name.

    Parameters
    ----------
    obsdf : pd.DataFrame
        DataFrame containing observation data with a DatetimeIndex.
    modeldatadf : pd.DataFrame
        DataFrame containing model data with a DatetimeIndex.
    tolerance : pd.Timedelta, optional
        Maximum allowed difference in time for matching, by default pd.Timedelta('5min').

    Returns
    -------
    pd.DataFrame
        Merged DataFrame with matched observations and model data.
    """
   
    modeldf_reset = modeldatadf.reset_index()
    obsdf_reset = obsdf.reset_index()

    # Sort by datetime for merge_asof
    modeldf_reset = modeldf_reset.sort_values('datetime')
    obsdf_reset = obsdf_reset.sort_values('datetime')

    # Merge using merge_asof on 'datetime', and exact match on 'obstype' and 'name'
    merged_df = pd.merge_asof(
            obsdf_reset,
            modeldf_reset,
            on='datetime',
            by=['obstype', 'name'],
            tolerance=tolerance,  # adjust tolerance as needed
            direction='nearest',
            suffixes=('_obs', '_model')
    )

    merged_df = merged_df.set_index(['datetime', 'obstype', 'name'])
    
    return merged_df