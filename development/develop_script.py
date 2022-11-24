#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:25:02 2022

@author: thoverga
"""

# import vlinder_toolkit
import os, sys
import pandas as pd
from pathlib import Path


lib_folder = Path(__file__).resolve().parents[1]
sys.path.append(str(lib_folder))
print(str(lib_folder))

from src import vlinder_toolkit

# % Import

testdatafile = os.path.join(str(lib_folder), 'tests', 'test_data',  'testdata_okt.csv')

static_data =  os.path.join(str(lib_folder), 'static_data', 'vlinder_metadata.csv')

lcz_map = os.path.join(str(lib_folder), 'physiograpy', 'lcz_filter_v1.tif')


#% Setup dataset
settings = vlinder_toolkit.Settings()
settings.update_settings(input_data_file=testdatafile,
                          input_metadata_file=static_data,
                          geotiff_lcz_file=lcz_map
                          )



dataset_hourly = vlinder_toolkit.Dataset()
dataset_5min = vlinder_toolkit.Dataset()
dataset_hourly.import_data_from_file(coarsen_timeres=True)
dataset_5min.import_data_from_file()
dataset_5min.apply_quality_control(obstype='temp', show_qc_info=False)
vlinder75 = dataset_5min.get_stations(['vlinder75'])
vlinder75_data = vlinder75['vlinder75'].temp

#dataset_hourly.apply_quality_control(obstype='temp')

#%%
from datetime import datetime

starttime=datetime(2022, 10,4) # 2022/09/04 00:00:00
endtime=datetime(2022,10,7,12,45)

stationnames = ['vlinder04', 'vlinder05']
plotdf = dataset.df[dataset.df['name'].isin(stationnames)]


# test = plotdf.loc['20221004':'20221005']

def datetime_subsetting(df, starttime, endtime):
    
    stand_format = '%Y-%m-%d %H:%M:%S'
    
    if isinstance(starttime, type(None)):
        startstring = None #will select from the beginning of the df
    else:
        startstring = starttime.strftime(stand_format)
    if isinstance(endtime, type(None)):
        endstring = None
    else: 
        endstring = endtime.strftime(stand_format)

    return df[startstring: endstring]


plotdf = datetime_subsetting(plotdf, None, None)


print(plotdf.index.min())
print(plotdf.index.max())



# plotdf = plotdf.loc[(plotdf >= starttime)]
# plotdf = plotdf.loc[(plotdf <= endtime)]


#%%

dataset.apply_quality_control(obstype='temp')


sta1 = dataset.get_stations(['vlinder09'])
data_sta1 = sta1['vlinder09'].df()
data_sta1_qc = sta1['vlinder09'].qc_labels_df['temp']
sta1.make_plot()


dataset.make_geo_plot()
#%%
import matplotlib.pyplot as plt
import geopandas as gpd
import mapclassify as mc

from mpl_toolkits.axes_grid1 import make_axes_locatable

variable='lcz'


default_settings=settings.plot_settings['spatial_geo']
timeinstance=dataset.df.index.min()

subdf = dataset.df.loc[timeinstance]
#create geodf
gdf = gpd.GeoDataFrame(subdf,
                       geometry=gpd.points_from_xy(subdf['lon'], subdf['lat']))

gdf = gdf[[variable, 'geometry']]


from matplotlib.colors import Normalize

def spatial_plot(gdf, variable, legend, use_quantiles, is_categorical, k_quantiles,
                 cmap, world_boundaries_map, figsize, extent, title, vmin, vmax):
    
    
    
    
    # Make color scheme
    if use_quantiles:
        # maybe better to use evenly spaced intervals rather than quantiles?
        scheme = 'equalinterval'
    else:
        scheme = None
        if (isinstance(vmin, type(None)) | isinstance(vmax, type(None))):
            vmin=gdf[variable].min()
            vmax=gdf[variable].max()
            
        
    
    #create figure object
    # ax = plt.subplot(111)
    fig, ax = plt.subplots(1, 1, figsize=figsize)

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
        # legend_kwds={'label': "colorbar label"}
    )
    
    #set extent
    ax.set_xlim(left=extent[0], right=extent[2])
    ax.set_ylim(bottom=extent[1], top=extent[3])

    
    ax.set_title(title)
    
    plt.show()
    return 

ax = spatial_plot(gdf=gdf,
                  variable=variable,
                  legend=True,
                  use_quantiles=True,
                  is_categorical=True,
                  k_quantiles=default_settings['n_for_categorical'],
                  cmap = default_settings['cmap'],
                  world_boundaries_map=settings.world_boundary_map,
                  figsize=default_settings['figsize'],
                  extent=default_settings['extent'],
                  title='hoi',
                  vmin=None,
                  vmax=None,
                  )

# #%

# #% Setup dataset
#%%
 # world_boundaries = gpd.read_file(world_boundaries_map)
# ax = world_boundaries.plot(world_boundaries, ax=ax, extent=extent)



# dataset = vlinder_toolkit.Dataset()
# dataset.import_data_from_file(coarsen_timeres=True)
#%% apply QC

