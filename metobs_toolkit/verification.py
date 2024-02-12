#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 17:00:45 2024

@author: thoverga
"""


import os
import sys
import logging
import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import folium
import matplotlib
import branca.colormap as brcm
from folium import plugins as folium_plugins

from metobs_toolkit.obstypes import Obstype
from metobs_toolkit.analysis import _make_time_derivatives

from metobs_toolkit.df_helpers import get_seasons, _make_time_derivatives, metadf_to_gdf

# from metobs_toolkit.plotting_functions import make_folium_html_plot
# from metobs_toolkit.dataset import Dataset
# from metobs_verif.verification_methods import get_basic_scores_dict
# from metobs_verif.modeloutput import Modelfield
# import metobs_verif.plotting as plotting

logger = logging.getLogger(__name__)


class Verification:
    def __init__(self, modeldata, dataset):
        self.modeldata = modeldata
        self.dataset = dataset  # dataset or station

        self.verifdf = None
        self.loc_to_gp_distance = None

        # define plot defaults
        self.modelcolor = "blue"
        self.obscolor = "orange"
        self.biascolor = "#D5573B"
        self.rmsecolor = "#885053"
        self.maecolor = "#777DA7"
        self.corcolor = "#94C9A9"

        # self._construct_verifdf()
        self._check_compatibility()

    def __repr__(self):
        return f"Verification object of \n {self.dataset} \n ------------- And ---------------- \n {self.modeldata} "

    def __str__(self):
        return f"Verification object of \n {self.dataset} \n ------------- And ---------------- \n {self.modeldata} "

    def _check_compatibility(self):
        # test location and names are equal between dataset and modeldat
        assert not (self.modeldata.df.empty), f"The Modeldata is empty."
        assert not (self.dataset.df.empty), f"The Dataset is empty."
        assert not (
            self.modeldata.metadf.empty
        ), f"The metadf of the modeldata is empty"
        assert not (self.dataset.metadf.empty), f"The metadf of the modeldata is empty"

        startobs = self.dataset.df.index.get_level_values("datetime").min()
        endobs = self.dataset.df.index.get_level_values("datetime").max()

        startmod = self.modeldata.df.index.get_level_values("datetime").min()
        endmod = self.modeldata.df.index.get_level_values("datetime").max()
        assert (
            startobs < endmod
        ), f"Start of observations {startobs} is after the end of the modeldata {endmod}"
        assert (
            startmod < endobs
        ), f"Start of modeldata {startmod} is after the end of the modeldata {endobs}"

        assert set(self.modeldata.metadf.index) == set(
            self.dataset.metadf.index
        ), f"The stationnames are not equal between the Modeldata and the Dataset."
        assert (
            self.modeldata.metadf["lat"] == self.dataset.metadf["lat"]
        ).all(), "The coordinates of the stations are not equal between the observations and the modeldata."
        assert (
            self.modeldata.metadf["lon"] == self.dataset.metadf["lon"]
        ).all(), "The coordinates of the stations are not equal between the observations and the modeldata."

    # def get_obs_obstypes(self):
    #     #retunr list
    #     all_obstypes = self.dataset.obstypes
    #     return [val for key, val in all_obstypes.items() if key in self.dataset.df.columns]

    # def get_model_fields(self):
    #     return self.modeldata.get_fields()

    # def _get_model_field_obstypes_dict(self):
    #     modelfields = self.get_model_fields()
    #     return {field: self.modeldata.data.attrs['_obstypes'][field] for field in modelfields}

    # def _construct_verifdf(self):
    #     # Test the observations
    #     assert not self.dataset.metadf.empty, f'{self.obs} has an empty metadf attribute.'
    #     assert isinstance(self.dataset.metadf, type(gpd.GeoDataFrame())), 'the observations metadf attribute is not a GeoDataFrame.'

    #     modeldf = self.modeldata.extract_values_at_2d_fields(geoseries = self.dataset.metadf['geometry'])

    #     # Merge observations and modeldata
    #     obsdf = self.dataset.df.reset_index()
    #     #TODO add time tollerance???
    #     verifdf= obsdf.merge(modeldf,
    #                           how='inner',
    #                           left_on=['name', 'datetime'],
    #                           right_on=['name', 'timestamp'])
    #     if verifdf.empty:
    #         sys.exit('No overlap found for the observations and the model point extractions.')

    #     # create a loc_to_gp_distance series
    #     self.loc_to_gp_distance = pd.Series(dict(zip(verifdf['name'], verifdf['tollerance_distances'])))
    #     verifdf = verifdf.drop(columns=['tollerance_distances'])

    #     verifdf = verifdf.reset_index(drop=True).set_index(['name', 'datetime'])
    #     # sort obstypes
    #     colslist = list(self.dataset.df.columns)
    #     colslist.extend(self.get_model_fields())
    #     verifdf = verifdf[colslist]
    #     self.verifdf = verifdf

    # def _get_corresponding_model_var(self, obstype):

    #     if isinstance(obstype, str):
    #         obstypestr = obstype
    #     elif isinstance(obstype, Obstype):
    #         obstypestr = obstype.name
    #     else:
    #         sys.exit(f'{obstype} not a string or Obstype.')

    #     field_obstype_dict = self._get_model_field_obstypes_dict()
    #     model_variabels = [key for key, val in field_obstype_dict.items() if val.name == obstypestr]
    #     assert len(model_variabels) > 0, f'There are no fields found with {observation_obstype} as a unit: {field_obstype_dict}'
    #     return model_variabels

    # # def get_point_verif(self, observation_obstype, model_variabels=None):
    # #     # check if obstype exist
    # #     present_obstypes = [obstype.name for obstype in self.get_obs_obstypes()]
    # #     assert observation_obstype in present_obstypes, f'{observation_obstype} not in the knonw obstypes: {present_obstypes}'

    # #     if isinstance(model_variabels, str):
    # #         model_variabels = [model_variabels]
    # #     if model_variabels is None:
    # #         #get all fields in the same observationtype
    # #         model_variabels = self._get_corresponding_model_var(obstype=observation_obstype)

    # #     #Check if all model variables are known
    # #     if not np.all([var in self.verifdf.columns for var in model_variabels]):
    # #         sys.exit(f'No all {model_variabels} are found in the verificatin table columns: {self.verifdf.columns}')

    # #     for var in model_variabels:
    # #         assert var in self.verifdf.columns, f'{var} not found in the verification table columns: {self.verifdf.columns}'

    # #     # COmpute basic scores
    # #     scores = {}
    # #     obstype = self.dataset.obstypes[observation_obstype]

    # #     fig, axs = plotting._make_plot_score_grid_axes()
    # #     fig.suptitle(self.modeldata._get_duration_representation_str(), fontsize=16)
    # #     for var in model_variabels:
    # #         fig, axs = plotting._make_plot_score_grid_axes()
    # #         fig.suptitle(f'{var} point verif for {self.modeldata._get_duration_representation_str()}', fontsize=10)

    # #         scores[var] = get_basic_scores_dict(model=self.verifdf[var],
    # #                                             obs=self.verifdf[observation_obstype])
    # #         axs[0] = plotting.plot_table(data=pd.Series(scores[var]),
    # #                                    ax=axs[0])
    # #         # Make scatter
    # #         axs[1] = plotting.simple_scatter_plot(data_x = self.verifdf[observation_obstype].values,
    # #                                      data_y = self.verifdf[var].values,
    # #                                      ax=axs[1],
    # #                                      y_label=f'{var} in {obstype.std_unit}',
    # #                                      x_label=f'{observation_obstype} in {obstype.std_unit}')

    # #         # Make histogram cmpariosn
    # #         axs[2] = plotting.comparison_hist_plot(df=self.verifdf[[obstype.name, var]].reset_index(drop=True),
    # #                                             ax=axs[2],
    # #                                             col_map_dict = {obstype.name: self.obscolor,
    # #                                                             var: self.modelcolor},
    # #                                             xlabel = f'{obstype.name} in {obstype.std_unit}',
    # #                                             ylabel = 'Frequency',
    # #                                             title='',
    # #                                             bins='auto',
    # #                                             orientation='vertical')

    def get_verification_analysis(
        self, observation_obstype="temp", model_obstype="temp_sfx", groupby=[""]
    ):

        # Check if obstypes exists
        assert (
            observation_obstype in self.dataset.obstypes.keys()
        ), f"{observation_obstype} not found in the known observational obstypes: {self.dataset.obstypes}"
        assert (
            model_obstype in self.modeldata.obstypes.keys()
        ), f"{model_obstype} not found in the known model obstypes: {self.modeldata.obstypes}"

        # Check if both obstype have thes ame standard unit
        assert (
            self.dataset.obstypes[observation_obstype].get_standard_unit()
            == self.modeldata.obstypes[model_obstype].get_standard_unit()
        ), f"The standard units of {observation_obstype} is not equal to {model_obstype}"

        # if isinstance(model_variabels, str):
        #     model_variabels = [model_variabels]
        # if model_variabels is None:
        #     #get all fields in the same observationtype
        #     model_variabels = self._get_corresponding_model_var(obstype=observation_obstype)

        scoringdf = self._get_grouped_point_scoring_metrics(
            observation_obstype=self.dataset.obstypes[observation_obstype],
            model_obstype=self.modeldata.obstypes[model_obstype],
            groupby=groupby,
        )
        # return scoringdf
        # import matplotlib.pyplot as plt
        # return scoringdf
        plotdf = scoringdf[["bias", "RMSE", "MAE", "cor"]]
        plotdf.index = plotdf.index.droplevel("variabel")

        fig, ax = plt.subplots()
        ax = simple_multiline_plot(
            df=plotdf,
            ax=ax,
            col_map_dict={
                "bias": self.biascolor,
                "RMSE": self.rmsecolor,
                "MAE": self.maecolor,
                "cor": self.corcolor,
            },
            xlabel=f"{groupby}",
            ylabel=f"/ or {self.dataset.obstypes[observation_obstype].get_plot_y_label()}",
            title=f"Point verif scores for {self.dataset.obstypes[observation_obstype].name} and {self.modeldata.obstypes[model_obstype].get_bandname()}.",
            add_zero=True,
        )
        return ax

        # Find out which situation is applicable
        # cur_Obstype = self.dataset.obstypes[observation_obstype]

        # #aggregated plot (thus over all groups if defined)
        # add_scoring_var_plot = False
        # if len(set(scoringdf.index.names)) > 1:
        #     add_scoring_var_plot = True

        # fig, axdict = plotting._make_verif_grid(variables=model_variabels,
        #                                         add_extra_row=add_scoring_var_plot)

    #     # Plot table
    #     tabledf = scoringdf.reset_index().copy()
    #     numcols = ['bias','RMSE', 'MAE', 'cor']
    #     for col in numcols:
    #         tabledf[col] = tabledf[col].apply(lambda x: f'{x:.3f}')
    #     plotting.plot_table(data=tabledf,
    #                         ax=axdict['scoringtable'])

    #     for var in model_variabels:
    #         #scatter plots
    #         axdict[var]['scatter'] = plotting.simple_scatter_plot(data_x = self.verifdf[observation_obstype].values,
    #                                              data_y = self.verifdf[var].values,
    #                                              ax=axdict[var]['scatter'],
    #                                              y_label=f'{var} in {cur_Obstype.std_unit}',
    #                                              x_label=f'{cur_Obstype.name} in {cur_Obstype.std_unit}')

    #         #Histogram plots
    #         axdict[var]['hist'] = plotting.comparison_hist_plot(df=self.verifdf[[cur_Obstype.name, var]].reset_index(drop=True),
    #                                                ax=axdict[var]['hist'],
    #                                                col_map_dict = {cur_Obstype.name: self.obscolor,
    #                                                                var: self.modelcolor},
    #                                                xlabel = f'{cur_Obstype.name} in {cur_Obstype.std_unit}',
    #                                                ylabel = 'Frequency',
    #                                                title='',
    #                                                bins='auto',
    #                                                orientation='vertical')

    #     # 1: Situation 1 : Only variables in index --> no category/time evolution

    #     # 2: Situation 2 : only variables and datetime in index -->
    #     if len(set(scoringdf.index.names)) > 1 :
    #         print('time evolving situation')

    #         for var in model_variabels:
    #             plotdf = scoringdf.xs(var, level='variabel').drop(columns=['N_verifpoints'])

    #             axdict[var]['score_var'] = plotting.simple_multiline_plot(df=plotdf,
    #                                                 ax=axdict[var]['score_var'],
    #                                                 col_map_dict = {'bias': self.biascolor,
    #                                                                 'RMSE': self.rmsecolor,
    #                                                                 'MAE': self.maecolor,
    #                                                                 'cor': self.corcolor},
    #                                                 xlabel= f'{list(plotdf.index.names)}',
    #                                                 ylabel='',
    #                                                 add_zero=True)

    #             # return ax

    #     # # 3: Situation 3 : Variables and categorical levels in index --->
    #     # else:
    #     #     print('general situation')

    def _construct_mod_obs_df(self, obs_obstype, mod_obstype, interp=True):

        obsdf = self.dataset.df[[obs_obstype.name]]

        # TODO, without interpolation
        if interp:
            mod_df = self.modeldata.sample_data_as(target=obsdf)
            mod_df = mod_df[[mod_obstype.name]]
        else:
            sys.exit("not implemented yet")

        combdf = mod_df.merge(obsdf, how="outer", left_index=True, right_index=True)
        return combdf

    def _get_grouped_point_scoring_metrics(
        self, observation_obstype, model_obstype, groupby=[""]
    ):
        # check if obstype exist
        # present_obstypes = [obstype.name for obstype in self.get_obs_obstypes()]
        # assert observation_obstype in present_obstypes, f'{observation_obstype} not in the knonw obstypes: {present_obstypes}'

        # format groupby
        if groupby == [""]:
            groupby = None
        if groupby is not None:
            assert np.array(
                [
                    grp_id
                    in [
                        "minute",
                        "hour",
                        "month",
                        "year",
                        "day_of_year",
                        "week_of_year",
                        "season",
                        "datetime",
                        "name",
                        "lcz",
                    ]
                    for grp_id in groupby
                ]
            ).all(), f"Unknonw groupid in {groupby}."

        df = self._construct_mod_obs_df(
            obs_obstype=observation_obstype, mod_obstype=model_obstype, interp=True
        )
        df = df.reset_index()
        # add time aggregated columns
        df = _make_time_derivatives(df=df, required="", get_all=True)

        # get lcz for stations
        lcz_mapper = self.dataset.metadf["lcz"].to_dict()
        df["lcz"] = df["name"].map(lcz_mapper)

        # rename validate
        # df = df.rename(columns={'datetime': 'validate'})

        # subset the dataframe
        relevant_columns = [observation_obstype.name, model_obstype.name]
        if groupby is not None:
            relevant_columns.extend(groupby)

        df = df[relevant_columns]

        # Compute scores per group
        scoringlist = []
        if groupby is not None:
            # for var in model_variabels:
            for idx, group in df.groupby(groupby):
                groupscores = get_basic_scores_dict(
                    model=group[model_obstype.name], obs=group[observation_obstype.name]
                )
                # convert scores to a dataframe
                groupscores.update({"group": idx})
                groupscores.update({"variabel": model_obstype.name})

                groupscoresdf = pd.Series(groupscores).to_frame().transpose()
                trg_index = ["group", "variabel"]
                groupscoresdf = groupscoresdf.set_index(trg_index)
                scoringlist.append(groupscoresdf)
        if groupby is None:
            # for var in model_variabels:
            groupscores = get_basic_scores_dict(
                model=df[model_obstype.name], obs=df[observation_obstype.name]
            )

            # convert scores to a dataframe
            groupscores.update({"variabel": model_obstype.name})
            groupscoresdf = pd.Series(groupscores).to_frame().transpose()
            groupscoresdf = groupscoresdf.set_index("variabel")
            scoringlist.append(groupscoresdf)

        if not bool(scoringlist):
            sys.exit(f"No groups could be made for {groupby}")
        scoringdf = pd.concat(scoringlist)
        return scoringdf

    def make_interactive_diff_plot(
        self,
        observation_obstype,
        model_obstype,
        # save=True,
        outputfile,
        # starttime=None,
        # endtime=None,
        agg=None,
        vmin=None,
        vmax=None,
        mpl_cmap_name="viridis",
        radius=13,
        fill_alpha=0.6,
        max_fps=4,
        outlier_col="red",
        ok_col="black",
        gap_col="orange",
        fill_col="yellow",
    ):
        """Make interactive geospatial plot with time evolution.

        This function uses the folium package to make an interactive geospatial
        plot to illustrate the time evolution.



        Parameters
        ----------
        obstype : str or metobs_toolkit.Obstype, optional
            The observation type to plot. The default is 'temp'.
        save : bool, optional
            If true, the figure will be saved as an html-file. The default is True.
        outputfile : str, optional
            The path of the output html-file. The figure will be saved here, if
            save is True. If outputfile is not given, and save is True, than
            the figure will be saved in the default outputfolder (if given).
            The default is None.
        starttime : datetime.datetime, optional
             Specifiy the start datetime for the plot. If None is given it will
             use the start datetime of the dataset, defaults to None.
        endtime : datetime.datetime, optional
             Specifiy the end datetime for the plot. If None is given it will
             use the end datetime of the dataset, defaults to None.
        vmin : numeric, optional
            The value corresponding with the minimum color. If None, the
            minimum of the presented observations is used. The default is None.
        vmax : numeric, optional
            The value corresponding with the maximum color. If None, the
            maximum of the presented observations is used. The default is None.
        mpl_cmap_name : str, optional
            The name of the matplotlib colormap to use. The default is 'viridis'.
        radius : int, optional
            The radius (in pixels) of the scatters. The default is 13.
        fill_alpha : float ([0;1]), optional
            The alpha of the fill color for the scatters. The default is 0.6.
        max_fps : int (>0), optional
            The maximum allowd frames per second for the time evolution. The
            default is 4.
        outlier_col : str, optional
            The edge color of the scatters to identify an outliers. The default is 'red'.
        ok_col : str, optional
            The edge color of the scatters to identify an ok observation. The default is 'black'.
        gap_col : str, optional
            The edge color of the scatters to identify an missing/gap
            observation. The default is 'orange'.
        fill_col : str, optional
            The edge color of the scatters to identify a fillded observation.
            The default is 'yellow'.

        Returns
        -------
        m : folium.folium.map
            The interactive folium map.

        Note
        -------
        The figure will only appear when this is runned in notebooks. If you do
        not run this in a notebook, make shure to save the html file, and open it
        with a browser.

        """

        # Check if obstypes exists
        assert (
            observation_obstype in self.dataset.obstypes.keys()
        ), f"{observation_obstype} not found in the known observational obstypes: {self.dataset.obstypes}"
        assert (
            model_obstype in self.modeldata.obstypes.keys()
        ), f"{model_obstype} not found in the known model obstypes: {self.modeldata.obstypes}"

        # Check if both obstype have thes ame standard unit
        assert (
            self.dataset.obstypes[observation_obstype].get_standard_unit()
            == self.modeldata.obstypes[model_obstype].get_standard_unit()
        ), f"The standard units of {observation_obstype} is not equal to {model_obstype}"

        Obstype_obs = self.dataset.obstypes[observation_obstype]
        Obstype_mod = self.modeldata.obstypes[model_obstype]

        # Construct differences timeseries
        combdf = self._construct_mod_obs_df(
            obs_obstype=Obstype_obs, mod_obstype=Obstype_mod, interp=True
        )

        combdf["diff"] = combdf[Obstype_obs.name] - combdf[Obstype_mod.name]
        combdf = combdf.dropna()
        # combdf = combdf[['diff']]
        combdf = combdf.reset_index()
        # Aggregate
        if agg == None:
            combdf = combdf
            combdf = combdf[["name", "datetime", "diff"]]
        else:
            combdf = _make_time_derivatives(df=combdf, required="", get_all=True)
            combdf = combdf.groupby(["name", agg])[["diff"]].mean()
            combdf = combdf.reset_index()

        # return combdf

        # add metadf

        combdf = combdf.merge(
            self.dataset.metadf, how="left", left_on="name", right_index=True
        )

        combdf = combdf.reset_index()

        # to gdf
        combgdf = metadf_to_gdf(combdf, crs=4326)

        # make time estimation
        est_seconds = combgdf.shape[0] / 2411.5  # normal laptop
        print(
            f'The figure will take approximatly (laptop) {"{:.1f}".format(est_seconds)} seconds to make.'
        )
        logger.info(
            f'The figure will take approximatly (laptop) {"{:.1f}".format(est_seconds)} seconds to make.'
        )

        # Making the figure
        m = make_folium_html_plot(
            gdf=combgdf,
            variable_column="diff",
            var_display_name=Obstype_obs.name,
            var_unit=Obstype_obs.get_standard_unit(),
            # label_column="label",
            # label_col_map=label_col_map,
            vmin=vmin,
            vmax=vmax,
            radius=radius,
            fill_alpha=fill_alpha,
            mpl_cmap_name=mpl_cmap_name,
            max_fps=int(max_fps),
        )
        # if save:
        logger.info(f"Saving the htlm figure at {outputfile}")
        m.save(outputfile)
        return m


def get_basic_scores_dict(model, obs):
    tot_scores = {}
    tot_scores["bias"] = _calc_bias(model, obs)
    tot_scores["N_verifpoints"] = (model - obs).dropna().shape[0]
    tot_scores["RMSE"] = _calc_rmse(model, obs)
    tot_scores["MAE"] = _calc_mae(model, obs)
    tot_scores["cor"] = _calc_cor(model, obs)
    return tot_scores


def _calc_bias(model, obs):
    return (model - obs).mean(skipna=True)


def _calc_rmse(model, obs):
    return np.sqrt(np.mean((model - obs) ** 2))


def _calc_mae(model, obs):
    return (model - obs).abs().mean(skipna=True)


def _calc_cor(model, obs):
    return model.corr(obs)


# =============================================================================
# plotters
# =============================================================================
# =============================================================================


def _make_handle(labelname, color, lw=1, linestyle="solid"):
    return Line2D([0], [0], color=color, label=labelname, lw=lw, linestyle=linestyle)


def _make_verif_grid(variables, add_extra_row=False, figsize=(16, 10)):

    fig = plt.figure(figsize=figsize)
    # fig.tight_layout()
    if add_extra_row:
        rowmax = 15
    else:
        rowmax = 10

    colmax = len(variables)
    spec = fig.add_gridspec(rowmax, colmax)

    axdict = {}
    ax_int_map = dict(zip(variables, range((len(variables)))))

    # scoring table
    axdict["scoringtable"] = fig.add_subplot(spec[0, :])

    # scatter plots
    for variable, i in ax_int_map.items():

        axdict[variable] = {}
        axdict[variable]["scatter"] = fig.add_subplot(spec[2:6, i])
        axdict[variable]["hist"] = fig.add_subplot(spec[6:10:, i])

        # add variation score axes
        if add_extra_row:
            axdict[variable]["score_var"] = fig.add_subplot(spec[10:, i])

    return fig, axdict


# def _make_plot_score_grid_axes(figsize = (16,10)):

#     fig = plt.figure(figsize=figsize)
#     # fig.tight_layout()

#     spec = fig.add_gridspec(31, 10)


#     ax_scores = fig.add_subplot(spec[1, 0:10])
#     ax_scatter = fig.add_subplot(spec[4:15, :])
#     ax_hist = fig.add_subplot(spec[16:, :])

#     return fig, [ax_scores, ax_scatter, ax_hist]


# # =============================================================================
# Plotting functions
# =============================================================================
def simple_scatter_plot(data_x, data_y, ax, y_label=None, x_label=None):

    plotdf = pd.DataFrame({"x": data_x, "y": data_y})

    ax = plotdf.plot.scatter(x="x", y="y", ax=ax)

    # add bissectrice
    minobs, maxobs = plotdf["x"].min(), plotdf["x"].max()
    ax.axline(
        (minobs, minobs),
        (maxobs, maxobs),
        linestyle="--",
        linewidth=3,
        color="black",
        zorder=0,
    )

    # set axes labels
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)

    return ax


def plot_table(data, ax, fontsize=8):

    pd.plotting.table(
        ax=ax,
        data=data,
        loc="center",
        cellLoc="center",
        fontsize=fontsize,
    )
    ax.set_axis_off()
    return ax


def comparison_hist_plot(
    df, ax, col_map_dict, xlabel, ylabel, title, bins="auto", orientation="vertical"
):

    # make plot
    cols = [col_map_dict[typ] for typ in df.columns]
    ax.hist(
        x=df,
        color=cols,
        bins=bins,
        histtype="bar",
        orientation=orientation,
        align="right",
    )
    # add legend
    ax.legend(handles=[_make_handle(nam, col) for nam, col in col_map_dict.items()])

    # add styling attrs
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    return ax


def simple_multiline_plot(df, ax, col_map_dict, xlabel, ylabel, title, add_zero=False):

    collist = [col_map_dict[col] for col in df.columns]
    ax = df.plot(ax=ax, color=collist)

    if add_zero:
        ax.axhline(y=0, color="black", linestyle="--")

    # add styling attrs
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    return ax


def make_plot(
    dxr,
    ax,
    title=None,
    grid=False,
    land=True,
    coastline=True,
    contour=False,
    levels=10,
    cbar_kwargs={},
    **kwargs,
):
    if contour:
        dxr.plot.contourf(ax=ax, levels=levels, cbar_kwargs=cbar_kwargs, **kwargs)

    else:
        dxr.plot(ax=ax, cbar_kwargs=cbar_kwargs, **kwargs)
    if land:
        ax.add_feature(cfeature.LAND)
        ax.add_feature(cfeature.BORDERS)
    if coastline:
        ax.add_feature(cfeature.COASTLINE)
    if grid:
        ax.gridlines(draw_labels=True, dms=True, x_inline=False, y_inline=False)
    ax.set_title(title)

    return ax


def _get_init_mapcenter(gdf):
    center = gdf.dissolve().centroid.iloc[0]
    return [center.y, center.x]


def make_folium_html_plot(
    gdf,
    variable_column,
    var_display_name,
    var_unit,
    # label_column,
    # label_col_map,
    vmin=None,
    vmax=None,
    radius=13,
    fill_alpha=0.6,
    mpl_cmap_name="viridis",
    max_fps=4,
    dt_disp_fmt="%Y-%m-%d %H:%M",
):

    # create a map
    m = folium.Map(
        location=_get_init_mapcenter(gdf),
        tiles="cartodbpositron",
        zoom_start=10,
        attr="<a href=https://github.com/vergauwenthomas/MetObs_toolkit </a>",
    )

    # add extra tiles
    folium.TileLayer("OpenStreetMap", overlay=False, name="OSM").add_to(m)
    # RIP free Stamen tiles
    # folium.TileLayer("Stamen Terrain", overlay=False, name='Terrain', show=False).add_to(m)
    # folium.TileLayer("stamentoner", overlay=False, name='Toner', show=False).add_to(m)

    # Coloring
    if vmin is None:
        vmin = gdf[variable_column].min()
    if vmax is None:
        vmax = gdf[variable_column].max()

    # Create colormap to display on the map
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
    mapper = matplotlib.cm.ScalarMappable(
        norm=norm, cmap=matplotlib.colormaps[mpl_cmap_name]
    )
    colormap = brcm.LinearColormap(
        colors=mapper.cmap.colors,
        index=None,
        vmin=vmin,
        vmax=vmax,
        caption=f"{var_display_name} ({var_unit}) colorbar",
    )

    # linear colorscale for values
    def map_value_to_hex(series, vmin, vmax, cmapname="viridis"):
        norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
        mapper = matplotlib.cm.ScalarMappable(
            norm=norm, cmap=matplotlib.colormaps[cmapname]
        )

        return series.apply(lambda x: str(matplotlib.colors.to_hex(mapper.to_rgba(x))))

    gdf["value_color"] = map_value_to_hex(
        gdf[variable_column], vmin, vmax, cmapname=mpl_cmap_name
    )

    # # check if all labels are defined
    # if (
    #     len(
    #         [
    #             lab
    #             for lab in gdf[label_column].unique()
    #             if lab not in label_col_map.keys()
    #         ]
    #     )
    #     > 0
    # ):
    #     sys.exit(
    #         f'Unmapped labels found: {[lab for lab in gdf["label"].unique() if lab not in label_col_map.keys()]}'
    #     )

    # gdf["label_color"] = gdf[label_column].map(label_col_map)

    # Serialize Data to Features
    def make_scater_feature(row):
        dtstring = pd.to_datetime(
            f'2100-01-01 {str(row["hour"]).zfill(2)}:00'
        ).strftime(dt_disp_fmt)
        coords = [[row["geometry"].x, row["geometry"].y]]
        popup_str = f" <b>{row['name']}</b>  <br> {'{:.1f}'.format(row[variable_column])} {var_unit} <br> {'dummy'}"

        features_instance = {
            "type": "Feature",
            "geometry": {
                "type": "MultiPoint",
                "coordinates": coords,
            },
            "properties": {
                "times": [dtstring],
                "popup": popup_str,
                "tooltip": f'{row["name"]}',
                "id": "geenidee",
                "icon": "circle",
                "iconstyle": {
                    "fillColor": row["value_color"],
                    "fillOpacity": fill_alpha,
                    "stroke": "false",
                    "radius": radius,
                    # "color": row["label_color"],
                },
            },
        }
        return features_instance

    features = gdf.apply(make_scater_feature, axis=1).to_list()

    # Add data to the map
    folium_plugins.TimestampedGeoJson(
        {
            "type": "FeatureCollection",
            "features": features,
        },
        period="PT1H",
        duration="PT1H",
        add_last_point=False,
        auto_play=False,
        loop=False,
        max_speed=max_fps,  # fps
        loop_button=True,
        date_options="YYYY/MM/DD HH:mm:ss",
        time_slider_drag_update=True,
    ).add_to(m)

    m.add_child(colormap)
    # add control
    folium.LayerControl().add_to(m)

    return m
