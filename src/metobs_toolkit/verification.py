
import logging

import numpy as np
import pandas as pd
import sklearn.metrics as metrics
import matplotlib.pyplot as plt


import metobs_toolkit.plot_collection as plotting
from metobs_toolkit.backend_collection.errorclasses import (
    MetObsModelDataError,
    MetObsObstypeNotFound)
from metobs_toolkit.backend_collection.datetime_aggregates import _get_time_derivates
from metobs_toolkit.backend_collection.datetime_aggregates import (
    possible_time_aggregates,
    _get_time_derivates
)


class Verification:
    def __init__(self, obj):
        self.obj = obj

    def _get_meta_agg_groups(self):
        return list(set(self.obj.metadf.columns))

    def _get_time_agg_groups(self):
        return possible_time_aggregates

    def _get_all_agg_groups(self):
        all_groups = (['name', 'datetime'] 
                     + self._get_meta_agg_groups()
                     + self._get_time_agg_groups())
        return all_groups
    

    def _create_verifdf(self,
                target_obstype:str = 'temp',
                trg_modelID: str | None = None,
                trg_bandname: str | None = None,
                shift_tolerance = pd.Timedelta('5min')) -> pd.DataFrame:

        #TODO check of obstype is present

        #Construct a long series, with values. If the label is not 'ok', the
        #value is replaced by Nan
        obsdf = (self.obj.df
                    .xs(target_obstype, level='obstype', drop_level=True)
                    .rename(columns={'value': 'observation'})
                )
        
        #use only the 'ok' labels for verification
        obsdf.loc[obsdf['label'] != 'ok', 'observation'] = np.nan

        obsdf = obsdf[['observation']]


        #For stations
        if 'name' not in obsdf.index.names:
            obsdf['name'] = self.obj.name
            obsdf = obsdf.set_index('name', append=True)

        #model
        fcdf = self.obj.modeldatadf
        #1: filter on obstype
        #TODO check of obstype is present
        fcdf = fcdf.xs(target_obstype, level='obstype', drop_level=True)
        
        #2. filter on modelID
        if trg_modelID is None:
            # Check if modelID is unique
            model_ids = fcdf.index.get_level_values('modelID').unique()
            if len(model_ids) != 1:
                raise ValueError(f"Multiple modelID values found: {model_ids}. Please specify trg_modelID.")
            trg_modelID = model_ids[0]

        fcdf = fcdf.xs(trg_modelID, level='modelID', drop_level=True)
        
        #3. filter on bandname
        if trg_bandname is None:
            # Check if modelID is unique
            bandnames = fcdf.index.get_level_values('bandname').unique()
            if len(bandnames) != 1:
                raise ValueError(f"Multiple bandnames values found: {bandnames}. Please specify trg_bandname.")
            trg_bandname = bandnames[0]

        fcdf = fcdf.xs(trg_bandname, level='bandname', drop_level=True)

        
        #For stations
        if 'name' not in fcdf.index.names:
            fcdf['name'] = self.obj.name
            fcdf = fcdf.set_index('name', append=True)
        
        #Format the fcdf
        fcdf = fcdf.rename(columns={'value':'fc'})[['fc']]


        verifdf = pd.merge_asof(left=fcdf.reset_index(),
                                right=obsdf.reset_index(),
                                on='datetime',
                                by='name',
                                tolerance=shift_tolerance)

        verifdf = verifdf.set_index(['datetime', 'name'])
        return verifdf


    def get_scoring_df(
            self,
            target_obstype='temp',
            groupby=['name'],
            trg_modelID=None,
            trg_bandname=None,
            shift_tolerance = pd.Timedelta('5min')
            ):
        
        #trg_columns
        trg_columns = list(set([*groupby]))

        #Construct the verifdf
        verifdf = self._create_verifdf(
                    target_obstype=target_obstype,
                    shift_tolerance=shift_tolerance,
                    trg_modelID=trg_modelID,
                    trg_bandname=trg_bandname)
        verifdf = verifdf.reset_index()


        #Metamerge
        metadf = self.obj.metadf.reset_index()
        #subset to target columns
        metadf = metadf[list(set(trg_columns+['name']).intersection(set(metadf.columns)))]
        verifdf = pd.merge(left=verifdf,
                        right=metadf,
                        how='left',
                        on='name')


        #datetime agg merge
        dt_aggdf = _get_time_derivates(datetimes=pd.DatetimeIndex(verifdf['datetime']))
        dt_aggdf = dt_aggdf[list(set(trg_columns+['datetime']).intersection(set(dt_aggdf.columns)))]
        #add dt to verifdf
        verifdf = pd.concat([verifdf.set_index('datetime'), dt_aggdf], axis=1)


        #compute diff ?

        #Compute scores
        scorings = {}
        for groupname, groupdf in verifdf.groupby(groupby):
            scorings[groupname] = {}
            #Dropna
            groupdf = groupdf.dropna(subset=['observation', 'fc'])
            if groupdf.empty:
                
                scorings[groupname]['rmse'] = np.nan
                scorings[groupname]['mae'] = np.nan
                scorings[groupname]['bias'] = np.nan
                scorings[groupname]['samplesize'] = 0
                continue
            #RMSE
            scorings[groupname]['rmse'] = metrics.root_mean_squared_error(
                y_true=groupdf['observation'],
                y_pred=groupdf['fc'])
            #MAE
            scorings[groupname]['mae'] = metrics.mean_absolute_error(
                y_true=groupdf['observation'],
                y_pred=groupdf['fc'])
            #BIAS
            scorings[groupname]['bias'] = (groupdf['fc'] - groupdf['observation']).mean()
            #COUNT
            scorings[groupname]['samplesize'] = groupdf.shape[0]

        scoringdf = pd.DataFrame.from_dict(scorings, orient='index')
        scoringdf.index.set_names(groupby, inplace=True)


        #construct the df to plot
        plotdf = (scoringdf
                .stack(dropna=False)
                .rename('value')
                )
        plotdf.index = plotdf.index.set_names('score', level=-1)
        return plotdf
        
    
    def get_scores_time_related(
            self,
            target_obstype='temp',
            trg_modelID=None,
            trg_bandname=None):
        #TODO
        pass


    def make_scatter_plot(
            self,
            target_obstype='temp',
            timeinstancev: pd.Timestamp | None = None,
            trg_modelID=None,
            trg_bandname=None):
        pass
        #TODO


    def get_obstypes(self, trg_obstype:str,
                     trg_modelID=None,
                     trg_bandname=None): #passed to find_modeltimeseries
        
        #Get the obstype of the observations
        if hasattr(self.obj, 'obstypes'):
            #is a dataset
            obs_obstype = self.obj.obstypes[trg_obstype]
            for sta in self.obj.stations:
                is_found=False
                try:
                    fc_obstype = sta.get_modeltimeseries(
                                        obstype=trg_obstype,
                                        trg_modelID=trg_modelID,
                                        trg_bandname=trg_bandname).modelobstype
                except (MetObsModelDataError, MetObsObstypeNotFound):
                    continue

                is_found=True
                break

            if not is_found:
                raise MetObsModelDataError(f'No modeldata found for {trg_obstype} in {trg_modelID} with {trg_bandname}.')
        else:  
            # obj is a station
            obs_obstype = self.obj.get_sensor(trg_obstype).obstype
            fc_obstype = self.obj.get_modeltimeseries(
                                    obstype=trg_obstype,
                                    trg_modelID=trg_modelID,
                                    trg_bandname=trg_bandname).modelobstype

        return obs_obstype, fc_obstype
    
    def plot_scoring(
            self,
            scoringdf, 
            xaxis_level: str,
            color_level: str | None = None,
            to_plot_value: str | None = 'rmse',
            to_plot_level='score',
            colmap: dict | None = None,
            ax=None,
            figkwargs={},
            draw_h_line_at: None | float = None,
    ):
    
        #subset to relevant values
        if to_plot_value is None:
            pass
        else:
            scoringdf = scoringdf.xs(to_plot_value,
                                    level=to_plot_level,
                                    drop_level=False)
        
        if color_level is None:
            color_level = to_plot_level

        #Check if all levels are used
        targeted_lvls = set([xaxis_level, color_level, to_plot_level])
        targeted_lvls.discard(None)


        if len(scoringdf.index.names) > len(targeted_lvls):
            raise ValueError(f'to many index lvls present for a 2D plot: {scoringdf.index.names} > {targeted_lvls}.')

        #2 to few index lvls
        if len(scoringdf.index.names) < len(targeted_lvls):
            raise ValueError(f'to few index lvls present for a 2D plot: {scoringdf.index.names} < {targeted_lvls}.')

        # 3 not the same index levels
        if set(scoringdf.index.names) != set(targeted_lvls):
            raise ValueError(f'the index names do not match : {scoringdf.index.names} != {targeted_lvls}.')


        #line or scatter plot
        # Determine plot kind based on xaxis_level type
        xaxis_values = scoringdf.index.get_level_values(xaxis_level)
        if pd.api.types.is_datetime64_any_dtype(xaxis_values):
            #for datettime
            kind = 'line'
        elif xaxis_level in self._get_time_agg_groups():
            #for time aggregates
            kind = 'line'
        else:
            #No time related xaxis
            kind = 'scatter'

        # color scheme:
        if colmap is None:
            colmap = plotting.create_categorical_color_map(
                    catlist=scoringdf.index.get_level_values(color_level).unique()
                )
        
        #axes
        if ax is None:
            ax = plotting.create_axes(**figkwargs)
        
        #iterate over color group
        for colgroupid, colgroupdf in scoringdf.groupby(by=color_level):
            colgroupdf = colgroupdf.droplevel(color_level)


            if kind=='line':
                colgroupdf.plot(kind=kind,
                                    color=colmap[colgroupid],
                                    ax=ax,
                                    label=colgroupid,
                                )
            else:
                #scatter requires a dataframe (not a series)
                if isinstance(colgroupdf, pd.Series):
                    colgroupdf= colgroupdf.to_frame().reset_index()
                colgroupdf.plot(kind=kind,
                                x=xaxis_level,
                                y='value',
                                color=colmap[colgroupid],
                                ax=ax,
                                label=colgroupid)
        
        if isinstance(draw_h_line_at, (int, float, np.integer, np.floating)):
            ax.axhline(draw_h_line_at, color='black', linestyle=':', zorder=0.1)
        
        #Title
        if to_plot_level is None:
            plotting.set_title(ax, f"all values from the {to_plot_level}-level, grouped per {color_level} and {xaxis_level}.")
        else: 
            plotting.set_title(ax, f"{to_plot_value} values grouped per {color_level} and {xaxis_level}.")
        
        plotting.set_ylabel(ax, to_plot_value)
        plotting.set_xlabel(ax, xlabel=xaxis_level)
        if kind=='scatter':
            plt.setp(ax.get_xticklabels(), rotation=90, ha='right')
        
        plotting.set_legend(ax=ax)
        plt.tight_layout()
        return ax


    def plot_timeseries_of_obs_and_fc(self,
                                      target_obstype='temp',
                                      ax=None,
                                      colmap=None,
                                      trg_modelID=None,
                                      trg_bandname=None,
                                      figkwargs={}):
        

        # Get verifdf
        verifdf = self._create_verifdf(target_obstype=target_obstype,
                                       trg_bandname=trg_bandname,
                                       trg_modelID=trg_modelID,
                                      )

        #Get obstypes (for plotting details)
        obs_obstype, fc_obstype = self.get_obstypes(
                        trg_obstype=target_obstype,
                        trg_modelID=trg_modelID,
                        trg_bandname=trg_bandname)
                        
        
        # Prepare colormap (by stations)
        if colmap is None:
            colmap = plotting.create_categorical_color_map(
                catlist=verifdf.index.get_level_values('name').unique()
            )

        # Setup axes: main + residuals
        if ax is None:
            fig, (ax_main, ax_resid) = plt.subplots(
                2, 1, sharex=True,
                gridspec_kw={'height_ratios': [3, 1]},
                figsize=figkwargs.get('figsize', (10, 6))
            )
        else:
            ax_main = ax
            fig = ax_main.figure
            ax_resid = fig.add_axes([ax_main.get_position().x0, 
                                     ax_main.get_position().y0 - 0.18, 
                                     ax_main.get_position().width, 
                                     ax_main.get_position().height * 0.3],
                                    sharex=ax_main)
        
        # Plot timeseries: obs (solid), fc (dotted)
        for name, group in verifdf.groupby('name'):
            color = colmap[name]
            # Observations
            ax_main.plot(
                    group.index.get_level_values('datetime'),
                    group['observation'],
                    label=f"{name} obs",
                    color=color,
                    linestyle='-')
            # Forecasts
            ax_main.plot(group.index.get_level_values('datetime'),
                         group['fc'],
                         label=f"fc @ {name}",
                         color=color,
                         linestyle=':')

        plotting.set_ylabel(ax_main, obs_obstype._get_plot_y_label())
        # plotting._drop
        plotting.set_title(ax_main, f"Timeseries of {obs_obstype.name} (obs & fc)\n (originates from {fc_obstype.model_band})")

        # Residuals: fc - observation
        for name, group in verifdf.groupby('name'):
            color = colmap[name]
            resid = group['fc'] - group['observation']
            ax_resid.plot(
                group.index.get_level_values('datetime'),
                resid,
                label=name,
                color=color,
                linestyle='-')
            
        ax_resid.axhline(0, color='grey', linewidth=0.8, linestyle='--')
        plotting.set_ylabel(ax_resid,f'fc - obs (in {fc_obstype.std_unit})')
        plotting.format_datetime_axes(ax=ax_resid)
        plotting.set_legend(ax=ax_resid)
        # ax_resid.legend(loc='upper left', ncol=2, fontsize='x-small')
        plotting.set_title(ax_resid, 'Residuals')

        plt.tight_layout()
        return ax_main, ax_resid


    def plot_value_cycles(
            self,
            target_obstype:str='temp',
            reference_station: str | None = None, #Always observations as reference !!! 
            colorby:str = 'LCZ',
            xaxis:str='hour', #or datetime?
            colmap: dict|None = None,
            ax=None,
            figkwargs = {},
            trg_modelID=None,
            trg_bandname=None,
            output_is_obs = True,
            ):
        
        # ---- Filter datasources -------
        #Get obstypes (for plotting details)
        obs_obstype, fc_obstype = self.get_obstypes(
                        trg_obstype=target_obstype,
                        trg_bandname=trg_bandname,
                        trg_modelID=trg_modelID)

        if output_is_obs:
            df = (self.obj.df
                    .xs(target_obstype, level='obstype', drop_level=True)
                        )
            #use only the 'ok' labels for verification
            df.loc[df['label'] != 'ok', 'value'] = np.nan
        else:
            df = (self.obj.modeldatadf
                        .xs(target_obstype, level='obstype', drop_level=True)
                        #TODO: subset modelID
                        #TODO: subset bandname
                        .droplevel(['modelID', 'bandname'])
                        # .rename(columns={'value': 'fc'})
                        )
                
        #Get reference values
        if reference_station is not None:
            assert reference_station in df.index.get_level_values('name'), f'{reference_station} is not found in the name level.'
            refrecords = (self.obj.df
                        .xs(target_obstype, level='obstype', drop_level=True)
                        .xs(reference_station, level='name', drop_level=True)
                        .rename(columns={'value': 'ref_value'})
                        )
            
            #merge by datetimes
            df = pd.merge(left=df.reset_index(),
                                right=refrecords.reset_index(),
                                how='left',
                                on='datetime')

            #compute the difference
            df['value'] = df['value'] - df['ref_value']
            #Set index
            df = df.set_index(['datetime', 'name'])
        
        df = df[['value']]

        # ---- Combine meta and time derivatives ---- 


        trg_columns = list(set([colorby, xaxis]))
        #Metamerge
        #subset to target columns
        metadf = self.obj.metadf.reset_index()
        metadf = metadf[list(set(trg_columns+['name']).intersection(set(metadf.columns)))]
        df = pd.merge(left=df.reset_index(),
                        right=metadf,
                        how='left',
                        on='name')


        #datetime agg merge
        dt_aggdf = _get_time_derivates(datetimes=pd.DatetimeIndex(df['datetime']))
        dt_aggdf = dt_aggdf[list(set(trg_columns+['datetime']).intersection(set(dt_aggdf.columns)))]
        df = pd.concat([df.set_index('datetime'), dt_aggdf], axis=1)

        # --- Group and aggregate data to xaxis and color groups ---
        plotdf = df.reset_index()[list(set([colorby, xaxis, 'value']))].groupby([colorby, xaxis]).mean()

        # --- Data plotting -----
                

        # color scheme:
        if colmap is None:
            colmap = plotting.create_categorical_color_map(
                    catlist=plotdf.index.get_level_values(colorby).unique()
                )

        #axes
        if ax is None:
            ax = plotting.create_axes(**figkwargs)

        #iterate over color group
        for colgroupid, colgroupdf in plotdf.groupby(by=colorby):
            colgroupdf = colgroupdf.droplevel(colorby)

            colgroupdf['value'].plot(kind='line',
                            color=colmap[colgroupid],
                            ax=ax,
                            label=colgroupid,
                            )
            
        if output_is_obs:
            plotting.set_ylabel(ax,f'{obs_obstype._get_plot_y_label()})')
            if reference_station is None:
                plotting.set_title(ax, f"{obs_obstype} observational values grouped per {colorby}.")
            else:
                 plotting.set_title(ax, f"{obs_obstype} observational differences with {reference_station}, grouped per {colorby}.")
        else: 
            plotting.set_ylabel(ax,f'{fc_obstype._get_plot_y_label()})')

            if reference_station is None:
                plotting.set_title(ax, f"{fc_obstype} forecasted values grouped per {colorby}.")
            else:
                plotting.set_title(ax, f"{fc_obstype} forecasted differences with {reference_station}, grouped per {colorby}.")

        
        plotting.set_xlabel(ax, xlabel=xaxis)
        plotting.set_legend(ax=ax)
        plt.tight_layout()

        return ax



