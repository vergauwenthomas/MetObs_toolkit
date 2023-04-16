#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 11:26:52 2022

@author: thoverga
"""



import pandas as pd
import math
import numpy as np
import geopandas as gpd
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec

from metobs_toolkit.geometry_functions import find_largest_extent
from mpl_toolkits.axes_grid1 import make_axes_locatable



def geospatial_plot(plotdf, variable, timeinstance, title, legend, vmin, vmax,
                    plotsettings, categorical_fields, static_fields,
                    display_name_mapper, world_boundaries_map ):

    #Load default plot settings
    default_settings=plotsettings['spatial_geo']

    #subset to obstype
    plotdf = plotdf[[variable, 'geometry']]

    #Subset to the stations that have coordinates
    ignored_stations = plotdf[plotdf['geometry'].isnull()]
    plotdf = plotdf[~plotdf['geometry'].isnull()]
    if plotdf.empty:
        # logger.error(f'No coordinate data found, geoplot can not be made. Plotdf: {plotdf}')
        print(f'No coordinate data found, geoplot can not be made. Plotdf: {plotdf}')
        return

    if not ignored_stations.empty:
        # logger.error(f'No coordinate found for following stations: {ignored_stations.index.to_list()}, these will be ignored in the geo-plot!')
        print(f'No coordinate found for following stations: {ignored_stations.index.to_list()}, these will be ignored in the geo-plot!')



    #make color scheme for field
    if variable in categorical_fields:
        is_categorical=True
        if variable == 'lcz':
            #use all available LCZ categories
            use_quantiles=False
        else:
            use_quantiles=True
    else:
        is_categorical=False
        use_quantiles=False


    #if observations extend is contained by default exten, use default else use obs extend
    use_extent=find_largest_extent(geodf=gpd.GeoDataFrame(plotdf),
                                   extentlist=default_settings['extent'])


    #Style attributes
    if isinstance(title, type(None)):
        if variable in static_fields:
            title = display_name_mapper[variable]
        else:
            dtstring = datetime.strftime(timeinstance, default_settings['fmt'])
            title = display_name_mapper[variable] + ' at ' + dtstring

    ax = _spatial_plot(gdf=plotdf,
                      variable=variable,
                      legend=legend,
                      use_quantiles=use_quantiles,
                      is_categorical=is_categorical,
                      k_quantiles=default_settings['n_for_categorical'],
                      cmap = default_settings['cmap'],
                      world_boundaries_map=world_boundaries_map,
                      figsize=default_settings['figsize'],
                      extent=use_extent,
                      title=title,
                      vmin=vmin,
                      vmax=vmax
                      )
    return ax

def _spatial_plot(gdf, variable, legend, use_quantiles, is_categorical, k_quantiles,
                 cmap, world_boundaries_map, figsize, extent, title, vmin, vmax):

    #TODO: docstring + beter positionion of the lengends + fix size of pies

    gdf = gpd.GeoDataFrame(gdf)
    fig, ax = plt.subplots(1, 1, figsize=figsize)


    # Make color scheme
    if use_quantiles:
        # maybe better to use evenly spaced intervals rather than quantiles?
        scheme = 'equalinterval'
    else:
        scheme = None
        if (isinstance(vmin, type(None)) | isinstance(vmax, type(None))):
            vmin=gdf[variable].min()
            vmax=gdf[variable].max()

    if is_categorical:
        legend_kwds={'loc': 'best'}
        vmin=None
        vmax=None
        cax=None
    else:
        legend_kwds=None
        divider = make_axes_locatable(ax)

        cax = divider.append_axes("right", size="5%", pad=0.1)

    #world map as underlayer
    world_boundaries = gpd.read_file(world_boundaries_map)
    world_boundaries.plot(ax=ax)


    # add observations as scatters
    gdf.plot(

        column=variable,
        scheme=scheme,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        # edgecolor='white',
        # linewidth=0.5,
        # scale='NUMBER OF PERSONS KILLED',
        # limits=(8, 24),
        categorical=is_categorical,
        legend=legend,
        # legend_var='scale',
        # legend_kwargs={'loc': 'upper left', 'markeredgecolor': 'black'},
        # legend_values=[2, 1], legend_labels=['2 Fatalities', '1 Fatality'],
        ax=ax,
        cax=cax,
        legend_kwds=legend_kwds
    )

    #set extent
    ax.set_xlim(left=extent[0], right=extent[2])
    ax.set_ylim(bottom=extent[1], top=extent[3])


    ax.set_title(title)

    return ax


def timeseries_plot(mergedf, obstype, title, xlabel, ylabel, colorby,
                    show_legend, show_outliers, plot_settings, gap_settings,
                    qc_info_settings
                    ):

    # plot_settings = plot_settings['time_series']

    # init figure
    fig, ax = plt.subplots(figsize=plot_settings['time_series']['figsize'])


    # subbset and cleanup data

    mergedf = mergedf[~mergedf.index.duplicated()]
    init_idx = mergedf.index

    if not show_outliers:
        if obstype+'_final_label' in mergedf.columns:
            # remove outliers from plotting data
            mergedf = mergedf[mergedf[obstype+'_final_label'] == 'ok']
            mergedf = mergedf[[obstype, obstype+'_final_label']]
        else:
            # no final label available, so plot all data
            mergedf = mergedf[[obstype]]




    if colorby=='label':
        # iterate over label groups
        col_mapper = _all_possible_labels_colormapper(plot_settings, qc_info_settings, gap_settings) #get color mapper
        outl_groups = mergedf.groupby(obstype+'_final_label')
        legenddict={}
        for outl_label, groupdf in outl_groups:
            outl_color = col_mapper[outl_label]

            #plot data
            if outl_label=='ok': #ok data as lines

                # add init_idx andf fill with nans (to avoid matplotlib interpolation)
                fill_idx = init_idx.to_frame().drop(groupdf.index)
                groupdf = pd.concat([groupdf, fill_idx])
                groupdf = groupdf.drop(columns=['name', 'datetime'], errors='ignore')
                groupdf.sort_index()
                plotdf = groupdf.reset_index().pivot(index='datetime', columns='name', values=obstype) #long to wide

                ax=plotdf.plot(kind='line', color=outl_color, ax=ax, legend=False,
                                zorder=plot_settings['time_series']['linezorder'],
                                linewidth=plot_settings['time_series']['linewidth'])


            elif outl_label in list(gap_settings['gaps_fill_info']['label'].items()): #fill gaps as dashed lines

                fill_idx = init_idx.to_frame().drop(groupdf.index)
                groupdf = pd.concat([groupdf, fill_idx])
                groupdf = groupdf.drop(columns=['name', 'datetime'], errors='ignore')
                groupdf.sort_index()
                plotdf = groupdf.reset_index().pivot(index='datetime', columns='name', values=obstype) #long to wide

                ax=plotdf.plot(kind='line', style='--', color=outl_color, ax=ax, legend=False,
                                zorder=plot_settings['dashedzorder'],
                                linewidth=plot_settings['linewidth'])


            else: #outliers as scatters
                plotdf = groupdf[obstype]
                plotdf.index = plotdf.index.droplevel('name')
                plotdf = plotdf.reset_index()
                ax=plotdf.plot(kind='scatter', x='datetime', y=obstype,
                               ax=ax, color=outl_color, legend=False,
                               zorder=plot_settings['time_series']['scatterzorder'],
                               s=plot_settings['time_series']['scattersize'])

            legenddict[outl_label] = outl_color

        # make legend
        if show_legend:
            custom_handles = []
            # add ok label at the top
            if 'ok' in legenddict.keys():
                custom_handles.append(Line2D([0], [0],
                                             color=legenddict['ok'],
                                             label='ok', lw=4))
                del legenddict['ok'] #remove ok key

            for gapfillmethod, label in gap_settings['gaps_fill_info']['label'].items():
                if label in legenddict.keys():
                    custom_handles.append(Line2D([0], [0],
                                                color=legenddict[gap_settings['gap_fill_info']['label']],
                                                label=f'gap filled ({gapfillmethod})', lw=4, linestyle='--'))
                    del legenddict[label] #remove key


            custom_scatters = [Line2D([0], [0], marker='o', color='w',
                                   markerfacecolor=col, label=lab, lw=1) for lab, col in legenddict.items()]

            custom_handles.extend(custom_scatters)
            ax.legend(handles=custom_handles)



    elif colorby=='name':
        plotdf = mergedf.reset_index().pivot(index='datetime', columns='name', values=obstype)
        ax=plotdf.plot(kind='line', legend=show_legend, ax=ax)


    #Set title
    ax.set_title(title)
    # ax.legend().set_title('')

    #Set x and y labels
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    return ax




def _make_pie_from_freqs(freq_dict, colormapper, ax,
                         plot_settings, title=None):

    # To dataframe
    stats = pd.Series(freq_dict, name='freq').to_frame()

    # make color mapper
    stats['color'] = stats.index.map(colormapper)

    if (stats['freq'] == 0.0).all():
        print('No occurences in sample.')
        #add a 100% no occurences to it, so it can be plotted
        no_oc_df = pd.DataFrame(index=['No occurences'],
                     data={'freq': [100.0],
                           'color':[plot_settings['color_mapper']['ok']]})
        stats = pd.concat([stats, no_oc_df])


    # Make pie and legend
    patches, text = ax.pie(stats['freq'], colors=stats['color'],
                           radius=plot_settings['pie_charts']['radius_big'])
    ax.legend(handles=patches, labels = [f'{l}, {s:0.1f}%' for l, s in zip(stats.index.to_list(),
                                                                   stats['freq'].to_list())],
              loc = plot_settings['pie_charts']['anchor_legend_big'])

    # add subtitle
    if not isinstance(title, type(None)):
        ax.set_title(title)

    return ax

def _outl_value_to_colormapper(plot_settings, qc_check_info):
    """ Make color mapper for the outlier LABELVALUES to colors. """
    color_defenitions = plot_settings['color_mapper']
    outl_name_mapper = {val['outlier_flag']: key for key, val in qc_check_info.items()}
    outl_col_mapper = {outl_type: color_defenitions[outl_name_mapper[outl_type]] for outl_type in outl_name_mapper.keys()}
    return outl_col_mapper


def _all_possible_labels_colormapper(plot_settings,qc_info_settings, gap_settings):
    """ Make color mapper for all LABELVALUES to colors. """
    color_defenitions = plot_settings['color_mapper']

    mapper = dict()

    #get QC outlier labels

    outl_col_mapper = _outl_value_to_colormapper(plot_settings=plot_settings, qc_check_info=qc_info_settings)
    mapper.update(outl_col_mapper)
    #get 'ok' and 'not checked'
    mapper['ok'] = color_defenitions['ok']
    mapper['not checked'] = color_defenitions['not checked']

    # update gap and missing timestamp labels
    mapper[gap_settings['gaps_info']['gap']['outlier_flag']] = color_defenitions['gap']
    mapper[gap_settings['gaps_info']['missing_timestamp']['outlier_flag']] = color_defenitions['missing_timestamp']

    #add gapfill
    for method, label in gap_settings['gaps_fill_info']['label'].items():
        mapper[label] = color_defenitions[method]

    return mapper


def qc_stats_pie(final_stats, outlier_stats, specific_stats,
                 plot_settings, qc_check_info):

    color_defenitions = plot_settings['color_mapper']
    # Define layout

    fig = plt.figure(figsize=plot_settings['pie_charts']['figsize'])
    fig.tight_layout()
    spec = fig.add_gridspec(4, 4, wspace=10)

    ax_thl = fig.add_subplot(spec[0,:2]) #top half left
    ax_thr = fig.add_subplot(spec[0,2:]) #top half right




    # 1. Make the finale label pieplot
    # make color mapper
    final_col_mapper = {
        'ok': color_defenitions['ok'],
        'QC outliers': color_defenitions['outlier'],
        'missing (gaps)': color_defenitions['gap'],
        'missing (individual)': color_defenitions['missing_timestamp']
        }

    _make_pie_from_freqs(final_stats, final_col_mapper,
                         ax_thl,
                         plot_settings,
                         'Final label frequencies')


    # 2. Make QC overview pie
    # make color mapper
    outl_col_mapper = _outl_value_to_colormapper(plot_settings, qc_check_info)

    _make_pie_from_freqs(outlier_stats, outl_col_mapper,
                         ax_thr,
                         plot_settings,
                         'Outlier performance')



    # 3. Make a specific pie for each indvidual QC + gap + missing

    spec_col_mapper = {
        'ok': color_defenitions['ok'],
        'not checked': color_defenitions['not checked'],
        'outlier': color_defenitions['outlier'],
        'gap': color_defenitions['gap'],
        'missing timestamp': color_defenitions['missing_timestamp']
        }

    i=0
    for checkname, stats in specific_stats.items():
        ax = fig.add_subplot(spec[math.floor(i/4)+1,i%4])
        _make_pie_from_freqs(stats, spec_col_mapper,
                             ax,
                             plot_settings,
                             checkname)
        i += 1



    plt.show()

    return

