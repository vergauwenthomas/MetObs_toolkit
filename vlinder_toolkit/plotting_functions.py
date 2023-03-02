#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 11:26:52 2022

@author: thoverga
"""
import pandas as pd
import math
import geopandas as gpd
from datetime import datetime
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec

from .geometry_functions import find_largest_extent
from .settings import Settings
from mpl_toolkits.axes_grid1 import make_axes_locatable


def geospatial_plot(plotdf, variable, timeinstance, static_fields, categorical_fields, title, legend, vmin, vmax):
    
    #Load default plot settings
    default_settings=Settings.plot_settings['spatial_geo']
    
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
            title = Settings.display_name_mapper[variable]
        else:
            dtstring = datetime.strftime(timeinstance, default_settings['fmt'])
            title = Settings.display_name_mapper[variable] + ' at ' + dtstring
    
    ax = spatial_plot(gdf=plotdf,
                      variable=variable,
                      legend=legend,
                      use_quantiles=use_quantiles,
                      is_categorical=is_categorical,
                      k_quantiles=default_settings['n_for_categorical'],
                      cmap = default_settings['cmap'],
                      world_boundaries_map=Settings.world_boundary_map,
                      figsize=default_settings['figsize'],
                      extent=use_extent,
                      title=title,
                      vmin=vmin,
                      vmax=vmax
                      )
    return ax



def timeseries_plot(dtseries, title, xlabel, ylabel,figsize):
    fig, ax = plt.subplots(figsize=figsize) 
    
    ax=dtseries.plot(ax=ax)
    
    ax.set_title(title)
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    return ax


def timeseries_comp_plot(plotdf, title, xlabel, ylabel, figsize):
        
    fig, ax = plt.subplots(figsize=figsize) 
    
    
    plotdf.plot(ax=ax, title=title)
    
    #titles and axis
    # ax.set_title=title
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    
    #legend
    ax.legend().set_title('')
    handles, labels = ax.get_legend_handles_labels()


    ax.legend(bbox_to_anchor=(1.0, -0.05), ncol=7)
    return ax



def spatial_plot(gdf, variable, legend, use_quantiles, is_categorical, k_quantiles,
                 cmap, world_boundaries_map, figsize, extent, title, vmin, vmax):
    
    gdf = gpd.GeoDataFrame(gdf)
    #create figure object
    # ax = plt.subplot(111)
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



def qc_stats_pie(valid_records, final_labels_df, qc_stats, figsize, title, obstype='temp'):
     
    valid_records_without_na = valid_records[valid_records[obstype].notna()]
    num_valid_records_without_na = len(valid_records_without_na)
    num_not_valid_observations = len(valid_records) - len(valid_records_without_na)
    if ((final_labels_df[obstype+'_final_label'] == 'missing timestamp (gap)').any()):
        num_gaps = final_labels_df[obstype+'_final_label'].value_counts()['missing timestamp (gap)']
    else:
        num_gaps = 0
    num_outliers = len(final_labels_df) - num_gaps
    
    colors = {'ok': 'green', 'not checked': 'orange', 'outlier': 'red'}
    names_1 = ['ok', 'outliers', 'gaps', 'not valid observations']
    colors_1 = ["green", "orange", "blue", "red"]
    values_1 = [num_valid_records_without_na, num_outliers, num_gaps, num_not_valid_observations]
    
    fig = plt.figure(figsize=(10,10))
    fig.tight_layout()
    
    spec = fig.add_gridspec(3, 4, wspace=10)
    
    ax_1 = fig.add_subplot(spec[0,:2])
    ax_2 = fig.add_subplot(spec[0,2:])
    
    patches1, texts1 = ax_1.pie(values_1, radius=1.5, colors = colors_1)
    percents_1 = [item/sum(values_1) * 100 for item in values_1]
    ax_1.legend(patches1, labels = [f'{l}, {s:0.1f}%' for l, s in zip(names_1, percents_1)], loc = (-0.25, 0.75))

    for i in range(len(qc_stats.columns)):
        data = list(qc_stats.iloc[:,i].values)
        names = list(qc_stats.index)
        ax = fig.add_subplot(spec[math.floor(i/4)+1,i%4])
        patches, texts = ax.pie(data, colors=list(colors.values()), radius = 8)
        ax.set_title(qc_stats.columns[i], pad=60, fontweight ="bold")
        ax.legend(patches, labels = [f'{l}, {s:0.1f}%' for l, s in zip(names, data)], loc = (0.25, 1))
       
    names_2 = list(final_labels_df[obstype+'_final_label'].value_counts().index)
    values_2 = list(final_labels_df[obstype+'_final_label'].value_counts().values)
    
    percents_2 = [item/sum(values_2) * 100 for item in values_2]
    patches2, texts2 = ax_2.pie(values_2, radius=1.5)
    ax_2.legend(patches2, labels = [f'{l}, {s:0.1f}%' for l, s in zip(names_2, percents_2)], loc = (-1, 0.5))
    plt.show()
    
    #colors = {'ok': 'green', 'not checked': 'orange', 'outlier': 'red'}
    #fig, ax = plt.subplots(2)
    
    #ax[0] = qc_stats.plot(kind='pie', subplots=True,
     #               colors = list(colors.values()), 
      #              autopct='%1.f%%', startangle=270, fontsize=10,
       #             layout=(3,3), figsize=(10,10),
        #            title = qc_stats.columns.to_list(),
         #           legend=False, ylabel='')
    
    #names = list(final_labels_df[obstype+'_final_label'].value_counts().index)
    #values = list(final_labels_df[obstype+'_final_label'].value_counts().values)
    #print(final_labels_df[obstype+'_final_label'].value_counts())
    # Creating plot
    #ax[1].pie(values, labels = names)
    #plt.show()
    #Set title
    #fig = ax[0][0][0][0].get_figure()
    #fig.suptitle(title)
    
    #Set legend
    
    return fig