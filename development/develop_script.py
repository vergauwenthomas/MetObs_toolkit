#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:25:02 2022
@author: thoverga
"""

#%%

# import metobs_toolkit
#
import os
import sys
from pathlib import Path
import pandas as pd


lib_folder = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(lib_folder))

# add testdata paths
from tests.push_test.test_data_paths import testdata




tmp_pickle=os.path.join(lib_folder, 'development', 'tmp', 'dev_pickle.pkl')

import metobs_toolkit




#%%

# # use_dataset = 'debug_wide'
use_dataset = 'Congo_single_station'

dataset = metobs_toolkit.Dataset()


dataset.update_settings(output_folder=None,
                        input_data_file=testdata[use_dataset]['datafile'],
                        input_metadata_file=testdata[use_dataset]['metadatafile'],
                        data_template_file=testdata[use_dataset]['template'],
                        metadata_template_file=testdata[use_dataset]['template'],
                        )


dataset.import_data_from_file(**testdata[use_dataset]['kwargs'])

dataset.coarsen_time_resolution(freq = testdata[use_dataset]['coarsen'])
# dataset.apply_quality_control()
#%%
dataset.make_geo_plot(boundbox=[7.0, -14, 47.3, 14])

#%%
# from datetime import datetime
# lon_min, lat_min, lon_max, lat_max
# extentlist = [2.260609, 49.25, 6.118359, 52.350618]




#%%
import geopandas as gpd
import matplotlib.pyplot as plt

# import geopandas as gpd
from shapely.geometry import box


def gpd_to_extent_box(geodf):
    return box(*geodf.total_bounds)


def extent_list_to_box(extentlist):
    return box(*extentlist)

def box_to_extent_list(bbox):
    return list(bbox.bounds)



def find_extend_of_geodf(geodf, lat_size=0.5, lon_size=0.5):
    """ This function will construct a bounding box for the plot.

        If the geodf contains more than one point, the bounding box is
        defined as the spatial span of the points.

        If the geodf contains only one point, a minimal span of lat_size,
        lon_size is created with the point at the centroid.
    """

    geodf_extent_box = gpd_to_extent_box(geodf)

    if geodf_extent_box.area != 0.0:
        # multiple stations can span the zoombox
        return geodf_extent_box

    # else: on station
    center_x, center_y = geodf_extent_box.centroid.x, geodf_extent_box.centroid.y

    minx, maxx = center_x - (lon_size/2.), center_x + (lon_size/2.)
    miny, maxy = center_y - (lat_size/2.), center_y + (lat_size/2.)

    print( minx, maxx, miny, maxy)

    return box(min([minx, maxx]), min([miny, maxy]), max([minx, maxx]), max([miny, maxy]))




def find_largest_extent(geodf, extentlist):
    # get extent of geodf as a box
    # geodf_extent_box = gpd_to_extent_box(geodf)
    geodf_extent_box = find_extend_of_geodf(geodf)

    # get extendbox of extendlist
    default_extent_box = extent_list_to_box(extentlist)

    # Check if default covers the geodf (default is belgium)
    if default_extent_box.covers(geodf_extent_box):
        return extentlist


    return box_to_extent_list(geodf_extent_box)



# geodf = dataset.metadf
# extentlist =[2.260609, 49.25, 6.118359, 52.350618]
# extent = find_largest_extent(geodf, extentlist)

# gdf = gpd.GeoDataFrame(geodf)
# fig, ax = plt.subplots(1, 1, figsize=(10,10))

# geodf.plot(ax=ax)


# # set extent
# ax.set_xlim(left=extent[0], right=extent[2])
# ax.set_ylim(bottom=extent[1], top=extent[3])
