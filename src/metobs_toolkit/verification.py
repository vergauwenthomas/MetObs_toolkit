import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Union, List, Dict
from sklearn.metrics import mean_squared_error, mean_absolute_error
from metobs_toolkit.backend_collection.datetime_aggregates import (
    possible_time_aggregates,
    _get_time_derivates
)
from metobs_toolkit.verif_collection.match_data import match_obs_and_model
from metobs_toolkit.backend_collection.loggingmodule import log_entry


class Verification:
    """
    A class for verification of model data against observations.
    
    Parameters
    ----------
    data_object : Dataset or Station
        The data object containing observations and model data.
    instantanious_tolerance : pd.Timedelta, optional
        Tolerance for matching observations and model data. Default is 5 minutes.
    """
    
    @log_entry
    def __init__(self, data_object: Union["Station", "Dataset"],
                 instantanious_tolerance: pd.Timedelta = pd.Timedelta('5min')):
        obsdf = data_object.df
        modeldatadf = data_object.modeldatadf

        if 'name' not in obsdf.index.names:
            obsdf = obsdf.assign(name=data_object.name)
        if 'name' not in modeldatadf.index.names:
            modeldatadf = modeldatadf.assign(name=data_object.name)

        self._verifdf = match_obs_and_model(
            obsdf=obsdf,
            modeldatadf=modeldatadf,
            tolerance=instantanious_tolerance
        )
        self.metadf = data_object.metadf

    @property
    def verifdf(self) -> pd.DataFrame:
        """Return the verification dataframe."""
        return self._verifdf

    @log_entry
    def traditional_scores(self,
                           groupby: List[str] = ['datetime'],
                           filter: Dict[str, str] = {'obstype': 'temp', 'label': 'ok'}) -> pd.DataFrame:
        """
        Calculate traditional verification scores (RMSE, MAE, Bias) grouped by specified columns.

        Parameters
        ----------
        groupby : list of str
            Columns to group by for score calculation.
        filter : dict
            Dictionary of column-value pairs to filter the data.

        Returns
        -------
        pd.DataFrame
            DataFrame with calculated scores for each group.
        """
        needed_columns = ['value_obs', 'value_model'] + groupby + list(filter.keys())
        df = self.verifdf.reset_index()

        # Merge metadata
        df = merge_metadata(df=df, metadf=self.metadf, target_columns=needed_columns)

        # Merge time derivatives
        df = merge_timederivatives(df=df, target_columns=needed_columns)

        # Validate required columns
        for col in needed_columns:
            if col not in df.columns:
                raise ValueError(f"Column '{col}' not found in verification dataframe, metadata, or possible_time_aggregates.")

        # Apply filter
        df = filter_verifdf(df, filter)

        # Drop non-overlap
        df = df.dropna(subset=['value_obs', 'value_model'])
        relev_columns = ['value_obs', 'value_model'] + groupby

        # Validate distinct categories
        distinct_cat = ['obstype', 'bandname', 'modelID']
        for cat in distinct_cat:
            if cat in df.columns and df[cat].nunique() > 1:
                raise ValueError(
                    f"Cannot calculate traditional scores for multiple categories in '{cat}'.\n"
                    f"Add it to the groupby or filter on it."
                )

        # Group by and calculate scores
        df = df[relev_columns].groupby(groupby).apply(
            lambda x: pd.Series({
                'rmse': np.sqrt(mean_squared_error(x['value_obs'], x['value_model'])),
                'mae': mean_absolute_error(x['value_obs'], x['value_model']),
                'modelbias': np.mean(x['value_model'] - x['value_obs']),
                'n': x['value_obs'].count()
            })
        )

        return df

    @log_entry
    def scatter_plot(self,
                     colorby: str = 'datetime',
                     filter: Dict[str, str] = {'obstype': 'temp'},
                     figsize: tuple = (10, 6),
                     bissectrice: bool = True,
                     **kwargs) -> plt.Axes:
        """
        Create a scatter plot of observations vs model values.

        Parameters
        ----------
        colorby : str
            Column to color the scatter points by.
        filter : dict
            Dictionary of column-value pairs to filter the data.
        figsize : tuple
            Figure size. Default is (10, 6).
        bissectrice : bool
            Whether to include the bissectrice line. Default is True.
        **kwargs
            Additional keyword arguments for the scatter plot.

        Returns
        -------
        plt.Axes
            The matplotlib Axes object for the plot.
        """
        needed_columns = ['value_obs', 'value_model', colorby] + list(filter.keys())
        df = self.verifdf.reset_index()

        # Merge metadata
        df = merge_metadata(df=df, metadf=self.metadf, target_columns=needed_columns)

        # Merge time derivatives
        df = merge_timederivatives(df=df, target_columns=needed_columns)

        # Validate required columns
        for col in needed_columns:
            if col not in df.columns:
                raise ValueError(f"Column '{col}' not found in verification dataframe, metadata, or possible_time_aggregates.")

        # Apply filter
        df = filter_verifdf(df, filter)

        # Create plot
        fig, ax = plt.subplots(figsize=figsize)
        scatter = ax.scatter(x=df['value_obs'], y=df['value_model'], c=df[colorby], **kwargs)

        if bissectrice:
            ax.plot(
                [df['value_obs'].min(), df['value_obs'].max()],
                [df['value_obs'].min(), df['value_obs'].max()],
                color='black', linestyle='--', zorder=2
            )

        #TODO: make axes reflect the band and obstype
        ax.set_xlabel('Observed Value')
        ax.set_ylabel('Model Value')
        ax.set_title('Scatter Plot of Observed vs Model Values')

        return ax

    def timeseries_plot(self):
        #TODO
        pass
        
    def timevariability_scores(self):
        #TODO
        pass

    def plot_scoringdf(self, scoringdf):
        #TODO
        pass
    
    
        
    
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


