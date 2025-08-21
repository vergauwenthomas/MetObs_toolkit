#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Geometry utility functions for spatial plotting and bounding box calculations.

Created on Fri Oct 21 09:13:01 2022

@author: thoverga
"""

import logging
from typing import List, Any
from shapely.geometry import box

# import geopandas as gpd

from metobs_toolkit.backend_collection.loggingmodule import log_entry

logger = logging.getLogger("<metobs_toolkit>")


@log_entry
def gpd_to_extent_box(geodf: Any) -> box:
    """
    Convert a GeoDataFrame to a shapely box with coordinates of the bounds.

    Parameters
    ----------
    geodf : GeoDataFrame
        The GeoDataFrame whose total bounds will be used to create the box.

    Returns
    -------
    shapely.geometry.box
        A shapely box representing the bounds of the GeoDataFrame.
    """
    return box(*geodf.total_bounds)


@log_entry
def extent_list_to_box(extentlist: List[float]) -> box:
    """
    Convert a list of coordinates to a shapely box.

    Parameters
    ----------
    extentlist : list of float
        List of coordinates [minx, miny, maxx, maxy] to define the box.

    Returns
    -------
    shapely.geometry.box
        A shapely box defined by the given coordinates.
    """
    return box(*extentlist)


@log_entry
def box_to_extent_list(bbox: box) -> List[float]:
    """
    Convert a shapely box to a list of its bound coordinates.

    Parameters
    ----------
    bbox : shapely.geometry.box
        The shapely box to convert.

    Returns
    -------
    list of float
        List of coordinates [minx, miny, maxx, maxy] representing the box bounds.
    """
    return list(bbox.bounds)


@log_entry
def find_extent_of_geodf(
    geodf: Any, lat_size: float = 1.0, lon_size: float = 1.0
) -> box:
    """
    Construct a bounding box for the plot based on the GeoDataFrame.

    If the GeoDataFrame contains more than one point, the bounding box is
    defined as the spatial span of the points. If the GeoDataFrame contains
    only one point, a minimal span of `lat_size` and `lon_size` is created
    with the point at the centroid.

    Parameters
    ----------
    geodf : GeoDataFrame
        The GeoDataFrame containing the spatial points.
    lat_size : float, optional
        The minimum latitude span for a single point, by default 1.0.
    lon_size : float, optional
        The minimum longitude span for a single point, by default 1.0.

    Returns
    -------
    shapely.geometry.box
        A shapely box representing the bounding box for the plot.
    """
    geodf_extent_box = gpd_to_extent_box(geodf)

    if geodf_extent_box.area != 0.0:
        # Multiple stations can span the zoom box
        return geodf_extent_box

    # Else: one station
    center_x, center_y = geodf_extent_box.centroid.x, geodf_extent_box.centroid.y

    minx, maxx = center_x - (lon_size / 2.0), center_x + (lon_size / 2.0)
    miny, maxy = center_y - (lat_size / 2.0), center_y + (lat_size / 2.0)

    return box(
        min([minx, maxx]), min([miny, maxy]), max([minx, maxx]), max([miny, maxy])
    )


@log_entry
def find_plot_extent(
    geodf: Any, user_bounds: List[float], default_extentlist: List[float]
) -> List[float]:
    """
    Find the most suitable plot bounds for a spatial plot.

    If `user_bounds` are valid, these are used. Otherwise, the bounds of the
    GeoDataFrame are computed. If these bounds are not contained by the default
    (e.g., Belgium) bounds, then the GeoDataFrame extent is used; otherwise, the default.

    Parameters
    ----------
    geodf : GeoDataFrame
        The geometry dataframe containing all the stations to plot.
    user_bounds : list of float
        List of bound coordinates [minx, miny, maxx, maxy] provided by the user.
    default_extentlist : list of float
        List of default bounds [minx, miny, maxx, maxy] (e.g., Belgium).

    Returns
    -------
    list of float
        A list of bounds for the spatial plot.
    """
    # Test if user_bounds is valid and can be used
    if bool(user_bounds):
        user_bounds = [float(x) for x in user_bounds]
        return user_bounds

    # Get extent of geodf as a box
    geodf_extent_box = find_extent_of_geodf(geodf)

    # Get extent box of extentlist
    default_extent_box = extent_list_to_box(default_extentlist)

    # Check if default covers the geodf (default is Belgium)
    if default_extent_box.covers(geodf_extent_box):
        return default_extentlist

    return box_to_extent_list(geodf_extent_box)
