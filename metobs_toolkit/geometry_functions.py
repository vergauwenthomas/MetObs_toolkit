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


def find_largest_extent(geodf, extentlist):
    # get extent of geodf as a box
    geodf_extent_box = gpd_to_extent_box(geodf)

    # get extendbox of extendlist
    default_extent_box = extent_list_to_box(extentlist)

    # Check if default covers the geodf
    if default_extent_box.covers(geodf_extent_box):
        return extentlist

    return geodf.total_bounds
