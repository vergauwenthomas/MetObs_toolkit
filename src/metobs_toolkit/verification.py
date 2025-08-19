from typing import Union, List, Dict, Callable, Literal, Optional
import pandas as pd
import matplotlib.pyplot as plt


from metobs_toolkit.verif_collection.match_data import match_obs_and_model

import metobs_toolkit.verif_collection.dataframe_func as df_func
import metobs_toolkit.verif_collection.scoring_metrics as metrics
import metobs_toolkit.verif_collection.agg_func as agg_mod
from metobs_toolkit.plot_collection.general_functions import (
    create_axes,
    set_legend,
    set_title,
    set_xlabel,
    set_ylabel,
    create_categorical_color_map,
    create_linestyle_map,
    format_datetime_axes,
)
from metobs_toolkit.backend_collection.loggingmodule import log_entry




class Verification:
    """
    A class for verification of model data against observations.
    """
    
    @log_entry
    def __init__(self, data_object: Union["Station", "Dataset"],
                 instantanious_tolerance: pd.Timedelta = pd.Timedelta('5min')):  # TYPO
        """Initialize the verification instance with matched obs and model data."""
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

    
    
    # ------------------------------------------
    #    DF construction
    # ------------------------------------------
    @property
    @log_entry
    def verifdf(self) -> pd.DataFrame:
        """Return the verification dataframe."""
        return self._verifdf

    @log_entry
    def traditional_scores(
        self,
        groupby: Optional[List[str]] = None,
        filterby: Optional[Dict[str, str]] = None,
    ) -> pd.DataFrame:
        """
        Calculate traditional verification scores (RMSE, MAE, Bias) grouped by specified columns.

        Parameters
        ----------
        groupby : list of str
            Columns to group by for score calculation.
        filterby : dict
            Column-value pairs to filter the data.

        Returns
        -------
        pd.DataFrame
            DataFrame with calculated scores for each group.
        """
        groupby = groupby or ['datetime']
        filterby = filterby or {'obstype': 'temp', 'label': 'ok'}
        needed_columns = ['value_obs', 'value_model'] + groupby + list(filterby.keys())
        df = df_func.construct_full_df(
            verifdf=self.verifdf,
            metadf=self.metadf,
            filterby=filterby,
            needed_columns=needed_columns
        )

        # Drop non-overlap
        df = df.dropna(subset=['value_obs', 'value_model'])
        relev_columns = ['value_obs', 'value_model'] + groupby

        # Validate distinct categories
        distinct_cat = ['obstype', 'bandname', 'modelID']
        for cat in distinct_cat:
            if cat in df.columns and df[cat].nunique() > 1:
                raise ValueError(
                    f"Cannot calculate traditional scores for multiple categories in '{cat}'. "
                    f"Add it to the groupby or filter on it."
                )

        # Group by and calculate scores
        df = df[relev_columns].groupby(groupby).apply(
            lambda x: pd.Series({
                'RMSE': metrics.rmse(x['value_obs'], x['value_model']),
                'MAE': metrics.mae(x['value_obs'], x['value_model']),
                'modelbias': metrics.bias(x['value_obs'], x['value_model']),
                'N_samples': metrics.N_samples(x['value_obs'], x['value_model'])
            })
        )
        return df

    @log_entry
    def group_values(
        self,
        groupby: Optional[List[str]] = None,
        filterby: Optional[Dict[str, str]] = None,
        agg_func: Optional[List[Callable]] = None,
    ) -> pd.DataFrame:
        """Group values and apply aggregation functions."""
        groupby = groupby or ['datetime']
        filterby = filterby or {'obstype': 'temp', 'label': 'ok'}
        agg_func = agg_func or [agg_mod.mean]
        return df_func.group_values(
            verifdf=self.verifdf,
            metadf=self.metadf,
            groupby=groupby,
            filterby=filterby,
            agg_func=agg_func
        )

    @log_entry
    def diurnal_cycle_values(
        self,
        obstype: str = 'temp',
        diurnal_group: str = 'name',
        refstation: Optional[str] = None,
        refvaluetype: Literal['obs', 'fc'] = 'obs',
        filterby: Optional[Dict[str, str]] = None,
        agg_func: Optional[List[Callable]] = None,
    ) -> pd.DataFrame:
        """Compute aggregated values for the diurnal cycle."""
        filterby = filterby or {'label': 'ok'}
        agg_func = agg_func or [agg_mod.mean, agg_mod.count_notna]

        if refstation is not None:
            verifdf = df_func.construct_verifdf_with_ref(
                verifdf=self.verifdf,
                refstation=refstation,
                refvaluetype=refvaluetype
            )
        else:
            verifdf = self.verifdf

        groupby = [diurnal_group, 'hour']
        filterby['obstype'] = obstype

        return df_func.group_values(
            verifdf=verifdf,
            metadf=self.metadf,
            groupby=groupby,
            filterby=filterby,
            agg_func=agg_func
        )

    # ------------------------------------------
    #    Plotting methods
    # ------------------------------------------
    
    @log_entry
    def plot_diurnal_cycle(
        self,
        diurnaldf: pd.DataFrame,
        plottype: Literal['obs', 'fc'] = 'obs',
        valuetype: str = 'mean',
        figsize: tuple = (15, 10),
    ) -> plt.Axes:
        """
        Plot the average diurnal cycle for observations or model values.

        Parameters
        ----------
        diurnaldf : pd.DataFrame
            Aggregated diurnal dataframe returned by diurnal_cycle_values.
        plottype : {'obs', 'fc'}
            Selects which values to plot.
        valuetype : str
            Column suffix to select from the MultiIndex column (e.g., 'mean').
        figsize : tuple
            Figure size.

        Returns
        -------
        plt.Axes
            The matplotlib axes used for the plot.
        """
        if plottype == 'obs':
            first_col = 'value_obs'
        elif plottype == 'fc':
            first_col = 'value_model'
        else:
            raise ValueError("plottype must be 'obs' or 'fc'")

        # get the name of the level(s) representing the groups
        group_name = list(set(diurnaldf.index.names) - {'hour'})

        # format values df for plotting
        plotdf = (
            diurnaldf[first_col]
            .stack(future_stack=True)
            .rename('value')
        )
        plotdf.index = plotdf.index.set_names('valuetype_col', level=-1)
        valueplotdf = plotdf.xs(valuetype, level='valuetype_col').unstack(level=group_name)

        # construct colordict
        colordict = create_categorical_color_map(valueplotdf.columns)

        # Create axes
        ax = create_axes(figsize=figsize)

        # plot each group
        for col in valueplotdf.columns:
            ax.plot(
                valueplotdf.index,
                valueplotdf[col],
                label=f"{col}@(Obs)" if plottype == 'obs' else f"{col}@(Model)",
                color=colordict[col],
                linestyle='-'
            )

        # styling
        ax = format_datetime_axes(ax=ax, set_diurnal_format=True)
        ax = set_legend(ax=ax)
        ax = set_ylabel(ax=ax, ylabel='Average values')
        ax = set_title(ax=ax, titlestr=f'Average Diurnal Cycle of {plottype}')
        return ax

    @log_entry
    def scatter_plot(
        self,
        colorby: str = 'datetime',
        filterby: Optional[Dict[str, str]] = None,
        figsize: tuple = (10, 6),
        bissectrice: bool = True,
        **kwargs
    ) -> plt.Axes:
        """
        Create a scatter plot of observations vs model values.

        Parameters
        ----------
        colorby : str
            Column to color the scatter points by.
        filterby : dict
            Column-value pairs to filter the data.
        figsize : tuple
            Figure size.
        bissectrice : bool
            Whether to include the bissectrice line.
        **kwargs
            Additional keyword arguments for the scatter plot.

        Returns
        -------
        plt.Axes
            The matplotlib Axes object for the plot.
        """
        filterby = filterby or {'obstype': 'temp'}
        needed_columns = ['value_obs', 'value_model', colorby] + list(filterby.keys())
        df = df_func.construct_full_df(
            verifdf=self.verifdf,
            metadf=self.metadf,
            filterby=filterby,
            needed_columns=needed_columns
        )
        
        #create color column
        colormap = create_categorical_color_map(df[colorby])
        mapped = df[colorby].map(colormap)
        df['color'] = mapped.where(mapped.notna(), None) #Nan to None
        
            
        # get obstype
        trg_obstype = self._get_obstype_from_df(df)

        # Create plot
        ax = create_axes(figsize=figsize)

        # main scatter
        ax.scatter(x=df['value_obs'], y=df['value_model'], c=df['color'], **kwargs)

        # add legend entries for each category in `colorby`
        # use empty plots so set_legend() can pick them up
        for cat in pd.unique(df[colorby]):
            ax.plot([], [], marker='o', linestyle='None', color=colormap[cat], label=str(cat))

        if bissectrice:
            ax.plot(
                [df['value_obs'].min(), df['value_obs'].max()],
                [df['value_obs'].min(), df['value_obs'].max()],
                color='black', linestyle='--', zorder=2
            )

        # TODO: make axes reflect the band and obstype
        ax = set_xlabel(ax=ax, xlabel=f'Observed Value of {trg_obstype} in {trg_obstype.std_unit}')
        ax = set_ylabel(ax=ax, ylabel=f'Model Value of {trg_obstype}')
        ax = set_title(ax=ax, titlestr='Scatter Plot of Observed vs Model Values')
        
        ax = set_legend(ax=ax)
        
        return ax

    @log_entry
    def timeseries_plot(
        self,
        colorby: str = 'name',
        filterby: Optional[Dict[str, str]] = None,
        figsize: tuple = (10, 6)
    ) -> plt.Axes:
        """
        Plot observed and model time series colored by a categorical column.

        Parameters
        ----------
        colorby : str
            Column name used to color/group the series.
        filterby : dict
            Column-value pairs to filter the data.
        figsize : tuple
            Figure size.

        Returns
        -------
        plt.Axes
            The matplotlib axes used for the plot.
        """
        filterby = filterby or {'obstype': 'temp', 'label': 'ok'}
        needed_columns = list(set(['value_obs', 'value_model', 'datetime'] + [colorby] + list(filterby.keys())))

        subdf = df_func.construct_full_df(
            verifdf=self.verifdf.reset_index(),
            metadf=self.metadf,
            filterby=filterby,
            needed_columns=needed_columns
        )

        # create axes
        ax = create_axes(figsize=figsize)

        # Create colorscheme
        colmap = create_categorical_color_map(subdf[colorby].unique())

        # plot observations
        plotdf = (
            subdf[['datetime', colorby, 'value_obs']]
            .dropna(subset=['value_obs'])
            .set_index(['datetime', colorby])
            .unstack()
        )
        plotdf.index = pd.to_datetime(plotdf.index)
        for col in plotdf['value_obs'].columns:
            ax.plot(
                plotdf.index,
                plotdf['value_obs'][col],
                label=f"{col}@(Obs)",
                color=colmap[col],
                linestyle='-'
            )

        # plot model values
        plotdf = (
            subdf[['datetime', colorby, 'value_model']]
            .dropna(subset=['datetime', colorby, 'value_model'])
            .set_index(['datetime', colorby])
            .unstack()
        )
        plotdf.index = pd.to_datetime(plotdf.index)
        for col in plotdf['value_model'].columns:
            ax.plot(
                plotdf.index,
                plotdf['value_model'][col],
                label=f"{col}@(Model)",
                color=colmap[col],
                linestyle=':'
            )

        ax = format_datetime_axes(ax=ax)
        ax = set_legend(ax=ax)
        return ax
    
    @log_entry  
    def timevariability_scores(self) -> None:
        # TODO
        pass
    
    @log_entry
    def scores_plot(
        self,
        scoresdf: pd.DataFrame,
        figsize: tuple = (10, 6),
        add_samplesize: bool = True,
        custom_colorscheme: Optional[Dict[str, str]] = None
    ) -> plt.Axes:
        """
        Plot one or more scores, optionally with sample size on a secondary axis.

        Parameters
        ----------
        scoresdf : pd.DataFrame
            Scores dataframe (optionally MultiIndexed).
        figsize : tuple
            Figure size.
        add_samplesize : bool
            If True, plot sample sizes on a twin axis.
        custom_colorscheme : dict
            Optional mapping from score name to color.

        Returns
        -------
        plt.Axes
            The matplotlib axes for the score plot.
        """
    

        ax = create_axes(figsize=figsize)
        
        
        cols_to_plot = [col for col in scoresdf.columns if col != 'N_samples']
        index_levels = scoresdf.index.nlevels
        
        colormap = create_categorical_color_map(
                catlist=cols_to_plot,
                user_color_defs=custom_colorscheme or {}
            )

        # Prepare x-axis and grouping
        if (index_levels == 1) or (index_levels > 2):
            for col in cols_to_plot:
                ax.plot(scoresdf.index, scoresdf[col],
                        label=col,
                        color=colormap(col),
                        linestyle='-')
                
        elif index_levels == 2:
            # create linestyle map
            linestyle_map = create_linestyle_map(scoresdf.index.get_level_values(1))
            for col in cols_to_plot:
                plotdf = scoresdf[col].unstack()  # unstack last level
                x = plotdf.index
                for second_level in plotdf.columns:
                    ax.plot(x, plotdf[second_level],
                            label=f"{col}@({second_level})",
                            color=colormap(col),
                            linestyle=linestyle_map[second_level])
        else:
            raise ValueError("scoresdf must have at least one index level.")

        if add_samplesize and ('N_samples' in scoresdf.columns):
            # Plot 'N_samples' below with shared x-axis
            ax2 = ax.twinx()
            if (index_levels == 1) or (index_levels > 2):
                x = scoresdf.index
                ax2.plot(x, scoresdf['N_samples'], color='gray', alpha=0.5, label='N', linestyle='dotted')
            elif index_levels == 2:
                plotdf = scoresdf['N_samples'].unstack()  # unstack last level
                for val in plotdf.columns:
                    ax2.plot(plotdf.index,
                             plotdf[val],
                             color='black',
                             alpha=0.5,
                             label=f"N_samples@({val})",
                             linestyle=linestyle_map[val])
            ax2.set_ylabel("Sample size")
            ax2.legend(loc='lower right')

        ax = set_xlabel(ax=ax, xlabel=scoresdf.index.names[0] if index_levels >= 1 else "Index")
        ax = set_ylabel(ax=ax, ylabel="Score Value")
        ax = set_legend(ax=ax)
        ax.grid(True)
        return ax
    
    # ------------------------------------------
    #    Helper methods
    # ------------------------------------------
    @log_entry
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
            








