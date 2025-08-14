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
from metobs_toolkit.verif_collection.plotting_helpers import (
    get_color,
    create_linestyle_map,
    default_colorscheme
)

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

        #verification dataframe
        self._verifdf = match_obs_and_model(
            obsdf=obsdf,
            modeldatadf=modeldatadf,
            tolerance=instantanious_tolerance
        )
        #Metadata
        self.metadf = data_object.metadf
        
        #Observation types
        self.obstypes = data_object.obstypes

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
                'RMSE': np.sqrt(mean_squared_error(x['value_obs'], x['value_model'])),
                'MAE': mean_absolute_error(x['value_obs'], x['value_model']),
                'modelbias': np.mean(x['value_model'] - x['value_obs']),
                'N_samples': x['value_obs'].count()
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

        #get obstype
        trg_obstype = self._get_obstype_from_df(df)
        


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
        ax.set_xlabel(f'Observed Value of {trg_obstype} in {trg_obstype.std_unit}')
        ax.set_ylabel(f'Model Value of {trg_obstype}')
        ax.set_title('Scatter Plot of Observed vs Model Values')

        return ax

    def timeseries_plot(self):
        #TODO
        pass
        
    def timevariability_scores(self):
        #TODO
        pass
    
    @log_entry
    def plot_scoringdf(self, scoresdf: pd.DataFrame,
                       figsize=(10, 6),
                       add_samplesize: bool = True,
                       custom_colorscheme: Dict[str, str] = {}) -> plt.Axes:
        
        colorscheme = default_colorscheme # Merge default and custom colorscheme
        colorscheme.update(custom_colorscheme)

        fig, ax = plt.subplots(figsize=figsize)
        
        cols_to_plot = [col for col in scoresdf.columns if col != 'N_samples']
        index_levels = scoresdf.index.nlevels
        # Prepare x-axis and grouping
        if (index_levels == 1) or (index_levels > 2):
                for col in cols_to_plot:
                        ax.plot(scoresdf.index, scoresdf[col],
                                label=col,
                                color=get_color(col),
                                linestyle='-')
                # No need for linestyle map

        elif index_levels == 2:
                
                #create linestyle map
                linestyle_map = create_linestyle_map(scoresdf.index.get_level_values(1))
                for col in cols_to_plot:
                        plotdf = scoresdf[col].unstack() #unstack last level
                        x = plotdf.index
                        for second_level in plotdf.columns:
                                ax.plot(x, plotdf[second_level],
                                        label=f"{col}@({second_level})",
                                        color=get_color(col),
                                        linestyle=linestyle_map[second_level]
                                        )
        if add_samplesize:
            # Plot 'N_samples' below with shared x-axis
            ax2 = ax.twinx()
            if 'N_samples' in scoresdf.columns:

                    if (index_levels == 1) or (index_levels > 2):
                            x= scoresdf.index
                            ax2.plot(x, scoresdf['N_samples'], color='gray', alpha=0.5, label='N', linestyle='dotted')
                    
                    elif index_levels == 2:
                            plotdf = scoresdf['N_samples'].unstack() #unstack last level
                            for val in plotdf.columns:
                                    ax2.plot(plotdf.index,
                                            plotdf[val],
                                            color='black',
                                            alpha=0.5,
                                            label=f"N_samples@({val})",
                                            linestyle=linestyle_map[val])
                    ax2.set_ylabel("Sample size")
                    ax2.legend(loc='lower right')

        ax.set_xlabel(scoresdf.index.names[0] if index_levels >= 1 else "Index")
        ax.set_ylabel("Score Value")
        ax.legend()
        ax.grid(True)
        return ax
    
    
    # ------------------------------------------
    #    Helper methods
    # ------------------------------------------
    def _get_obstype_from_df(self, df: pd.DataFrame) -> "Obstype":
        """
        Extract the obstype from the verification dataframe.

        Parameters
        ----------
        df : pd.DataFrame
            The verification dataframe.

        Returns
        -------
        Obstype
            The obstype object.
        """
        if 'obstype' in df.columns:
            obstype = df['obstype'].unique()
            if len(obstype) > 1:
                raise ValueError("Cannot extract a single obstype from multiple values.")
            else:
                target_obstype_str = obstype[0]
                if target_obstype_str in self.obstypes:
                    return self.obstypes[target_obstype_str]
                else:
                    raise ValueError(f"Obstype '{target_obstype_str}' not found in the verification.Obstypes.")
        else:
            raise ValueError("Column 'obstype' not found in the verification dataframe.")
        
    

        
    
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


