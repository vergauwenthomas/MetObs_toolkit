#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 09:13:01 2022

@author: thoverga
"""
# import geopandas as gpd
from shapely.geometry import box


def gpd_to_extent_box(geodf):
    return box(*geodf.total_bounds)


def extent_list_to_box(extentlist):
    return box(*extentlist)

def box_to_extent_list(bbox):
    return list(bbox.bounds)



def find_extend_of_geodf(geodf, lat_size=1., lon_size=1.):
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

    return box(min([minx, maxx]), min([miny, maxy]), max([minx, maxx]), max([miny, maxy]))




def find_plot_extent(geodf, user_bounds, default_extentlist):
    # test if user_bounds is valid and can be used
    if bool(user_bounds):
        user_bounds = [float(x) for x in user_bounds]

        return user_bounds


    # get extent of geodf as a box
    # geodf_extent_box = gpd_to_extent_box(geodf)
    geodf_extent_box = find_extend_of_geodf(geodf)

    # get extendbox of extendlist
    default_extent_box = extent_list_to_box(default_extentlist)

    # Check if default covers the geodf (default is belgium)
    if default_extent_box.covers(geodf_extent_box):
        return default_extentlist


    return box_to_extent_list(geodf_extent_box)




