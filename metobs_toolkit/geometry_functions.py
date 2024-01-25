#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 09:13:01 2022

@author: thoverga
"""
# import geopandas as gpd
from shapely.geometry import box


def gpd_to_extent_box(geodf):
    """Convert GeoDataFrame to a box with coordinates of the bounds."""
    return box(*geodf.total_bounds)


def extent_list_to_box(extentlist):
    """Convert list of coordinates to a shapely box."""
    return box(*extentlist)


def box_to_extent_list(bbox):
    """Convert shapely box to a list of the bound coordinates."""
    return list(bbox.bounds)


def find_extend_of_geodf(geodf, lat_size=1.0, lon_size=1.0):
    """Construct a bounding box for the plot.

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

    minx, maxx = center_x - (lon_size / 2.0), center_x + (lon_size / 2.0)
    miny, maxy = center_y - (lat_size / 2.0), center_y + (lat_size / 2.0)

    return box(
        min([minx, maxx]), min([miny, maxy]), max([minx, maxx]), max([miny, maxy])
    )


def find_plot_extent(geodf, user_bounds, default_extentlist):
    """Find the most suitable plot bounds for spatial plot.

    If the user_bounds are valid, these are used. Else the bounds of the goedf
    computed. If these bounds are not contained by the default (Belgium) bounds
    than the geodf extend is used else the default.
    Parameters
    ----------
    geodf : geopandas.geoDataFrame
        The geometry dataframe containing all the stations to plot.
    user_bounds : list
        List of bound coordinates.
    default_extentlist : list
        List of default bounds (Belgium).

    Returns
    -------
    list
        A list of bounds for the spatial plot.

    """
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
