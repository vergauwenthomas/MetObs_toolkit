import pandas as pd
import numpy as np
from typing import Union, List, Dict


from metobs_toolkit.backend_collection.loggingmodule import log_entry
from metobs_toolkit.backend_collection.datetime_aggregates import (
    possible_time_aggregates,
    _get_time_derivates
)


@log_entry
def filter_verifdf(df: pd.DataFrame, filter_dict: Dict[str, str]) -> pd.DataFrame:
    """Filter the verification dataframe based on the provided filter dictionary."""
    for key, value in filter_dict.items():
        df = df[df[key] == value]
    return df


@log_entry
def merge_metadata(df: pd.DataFrame, metadf: pd.DataFrame, target_columns: List[str]) -> pd.DataFrame:
    """Merge metadata columns into the dataframe."""
    for col in target_columns:
        if col in df.columns:
            continue
        if col in metadf.columns:
            df = df.merge(metadf[[col]], how='left', left_on='name', right_index=True)
    return df


@log_entry
def merge_timederivatives(df: pd.DataFrame, target_columns: List[str]) -> pd.DataFrame:
    """Merge time derivative columns into the dataframe."""
    if any(col in possible_time_aggregates for col in target_columns):
        time_derivates = _get_time_derivates(df['datetime'])
        target_time_aggregates = [col for col in target_columns if col in possible_time_aggregates]
        df = df.merge(time_derivates[target_time_aggregates], how='left', left_on='datetime', right_index=True)
    return df


def construct_full_df(verifdf: pd.DataFrame,
                      metadf: pd.DataFrame,
                      filterby: Dict[str, str] = None,
                      needed_columns: List[str] = []) -> pd.DataFrame:
    
    df = verifdf.reset_index()

    # Merge metadata
    df = merge_metadata(df=df, metadf=metadf, target_columns=needed_columns)

    # Merge time derivatives
    df = merge_timederivatives(df=df, target_columns=needed_columns)
    
    # Validate required columns
    for col in needed_columns:
        if col not in df.columns:
            raise ValueError(f"Column '{col}' not found in verification dataframe, metadata, or possible_time_aggregates.")
    
    #filter the dataframe
    subdf = filter_verifdf(df, filterby)

    return subdf

def group_values(verifdf, metadf,
                groupby: List[str] = ['datetime'],
                filterby: Dict[str, str] = {'obstype': 'temp', 'label': 'ok'},
                agg_func=[np.mean]) -> pd.DataFrame:
  
    needed_columns = ['value_obs', 'value_model'] + groupby + list(filterby.keys())
    df = verifdf.reset_index()

    # Merge metadata
    df = merge_metadata(df=df, metadf=metadf, target_columns=needed_columns)

    # Merge time derivatives
    df = merge_timederivatives(df=df, target_columns=needed_columns)

    # Validate required columns
    for col in needed_columns:
        if col not in df.columns:
            raise ValueError(f"Column '{col}' not found in verification dataframe, metadata, or possible_time_aggregates.")

    # Apply filter
    df = filter_verifdf(df, filterby)
    
    relev_columns = ['value_obs', 'value_model'] + groupby

    #groupby
    groupdf = df[relev_columns].groupby(groupby).agg(agg_func)
    
    return groupdf



def construct_verifdf_with_ref(
    verifdf: pd.DataFrame,
    refstation: str,
    refvaluetype: Union['obs', 'fc'] = 'obs',
    tolerance: Union[str, pd.Timedelta] = pd.Timedelta('4min')) -> pd.DataFrame:

    
    if refvaluetype == 'obs':
        first_col = 'value_obs'
    elif refvaluetype == 'fc':
        first_col = 'value_model'
    else:
        raise ValueError("refvaluetype must be 'obs' or 'fc'")

    # Create a reference series
    refseries = verifdf[[first_col]].xs(refstation, level='name', drop_level=False)

    # Reset index for merging
    orig_idexes = list(verifdf.index.names)
    verifdf = verifdf.reset_index()
    refseries = refseries.reset_index()[['datetime', 'obstype', first_col]]

    # Merge verifdf with refseries
    verifdf = pd.merge_asof(right=verifdf,
                            left=refseries,
                            on='datetime',
                            by=['obstype'],
                            tolerance=tolerance,
                            direction='nearest',
                            suffixes=('', '_ref')
                            )
    
    verifdf = verifdf.set_index(orig_idexes)
    
    # Adjust values based on reference
    verifdf['value_obs'] -= verifdf[f'{first_col}_ref']
    verifdf['value_model'] -= verifdf[f'{first_col}_ref']
    
    #drop the reference column
    verifdf = verifdf.drop(columns=[f'{first_col}_ref'])
    
    return verifdf