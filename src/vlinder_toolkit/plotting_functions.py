#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 11:26:52 2022

@author: thoverga
"""
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
# import geoplot as gplt
import mapclassify as mc
from mpl_toolkits.axes_grid1 import make_axes_locatable

def timeseries_plot(dtseries, title, xlabel, ylabel,figsize):
    fig, ax = plt.subplots(figsize=figsize) 
    
    ax=dtseries.plot(ax=ax)
    
    ax.set_title(title)
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    return ax


def timeseries_comp_plot(plotdf, title, xlabel, ylabel, figsize):
        
    fig, ax = plt.subplots(figsize=figsize) 
    
    
    plotdf.plot(ax=ax)
    
    
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend().set_title('')
    return ax

# def spatial_plot(gdf, variable, legend, proj_type, use_quantiles, k_quantiles,
#                  cmap, world_boundaries_map, figsize, extent, title, vmin, vmax):
    
#     # get projection of plot
#     if proj_type == 'Orthographic':
#         proj=gplt.crs.Orthographic()
#     elif proj_type == 'AlbersEqualArea':
#         proj=gplt.crs.AlbersEqualArea()
#     else:
#         print("The projection: ", proj_type, " is not yet implemented!")
#         return
    
    
#     # Make color scheme
#     if use_quantiles:
#         #maybe better to use evenly spaced intervals rather than quantiles?
#         scheme = mc.Quantiles(gdf[variable],
#                               k=k_quantiles)
#     else:
#         scheme = None
#         if (isinstance(vmin, type(None)) | isinstance(vmax, type(None))):
#             vmin=gdf[variable].min()
#             vmax=gdf[variable].max()
            
#         norm=Normalize(vmin=vmin,
#                        vmax=vmax)
    
    
#     #create figure object
#     ax = plt.subplot(111, projection=proj)
    
#     # add observations as scatters
#     ax = gplt.pointplot(
#         df=gdf,
#         projection=proj,
#         hue=variable, 
#         scheme=scheme,
#         cmap=cmap,
#         norm=norm,
#         # edgecolor='white', 
#         # linewidth=0.5,
#         figsize=figsize,
#         # scale='NUMBER OF PERSONS KILLED',
#         # limits=(8, 24),
#         legend=legend,
#         # legend_var='scale',
#         # legend_kwargs={'loc': 'upper left', 'markeredgecolor': 'black'},
#         # legend_values=[2, 1], legend_labels=['2 Fatalities', '1 Fatality'],
#         ax=ax,
#     )
    
#     world_boundaries = read_file(world_boundaries_map)
#     ax = gplt.polyplot(world_boundaries, ax=ax, extent=extent)
    
    
#     ax.set_xlim(left=extent[0], right=extent[2])
#     ax.set_ylim(bottom=extent[1], top=extent[3])
    
#     ax.set_title(title)
#     return ax

def spatial_plot(gdf, variable, legend, use_quantiles, is_categorical, k_quantiles,
                 cmap, world_boundaries_map, figsize, extent, title, vmin, vmax):
    
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