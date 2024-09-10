#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions that are used for GEE interactions.

@author: thoverga
"""

import sys
import logging
from time import sleep
import pytz
import pandas as pd
import ee

from metobs_toolkit.df_helpers import init_multiindexdf

logger = logging.getLogger(__name__)


def connect_to_gee(**kwargs):
    """
    Setup authentication for the use of the GEE Python API.

    For a fresh kernel, without stored credentials, a prompt/browser window
    will appear with further instructions for the authentication.


    Parameters
    ----------
    **kwargs : Kwargs passed to ee.Authenticate()
        Kwargs are only used by the user, for resetting the gee connection. See
        the Note below.

    Returns
    -------
    None.

    Note
    ------
    Upon calling, this function assumes you have a Google developers account,
    and a project with the Google Earth Engine API enabled.
    See the * Using Google Earth Engine * page for more info.

    Note
    ------
    During the Authentication, you will be asked if you want a read-only scope.
    A read-only scope is sufficient when the data is transferred directly to your
    machine (small data transfers), but will not be sufficient when extracting
    large amounts of data (typical for extracting Modeldata). This is because
    modeldata is written directly to your Google Drive, and therefore
    the read-only scope is insufficient.

    Note
    ------
    Due to several reasons, an EEExeption may be thrown. This is
    likely because of an invalid credential file. To fix this, you
    can update your credential file, and specify a specific authentication method.
    We found that the "notebook" authentication method works best for most users.

    Here is an example on how to update the credentials:

    .. code-block:: python

        import metobs_toolkit

        metobs_toolkit.connect_to_gee(force=True, #create new credentials
                                      auth_mode='notebook', # 'notebook', 'localhost', 'gcloud' (requires gcloud installed) or 'colab' (works only in colab)
                                      )

    """

    if bool(kwargs):  # kwargs are always passed by user, so reinitialize
        ee.Authenticate(**kwargs)
        ee.Initialize()
        return

    if not ee.data._credentials:  # check if ee connection is initialized
        ee.Authenticate()
        ee.Initialize()
    return


# =============================================================================
# Object convertors
# =============================================================================


def _datetime_to_gee_datetime(datetime):
    # ASSUME DATETIME IN UTC!!!
    # utcdt = datetime.astimezone(pytz.utc)#This will c
    logger.debug(datetime.replace(tzinfo=None))
    return ee.Date(datetime.replace(tzinfo=None))


def get_ee_obj(model, target_bands=[], force_mosaic=None):
    """Get an image from a GEE object."""

    # get the dataset
    if model.is_image:
        obj = ee.Image(model.location).select(model.band_of_use)
    else:
        obj = ee.ImageCollection(model.location)

        # filter the bands
        if bool(target_bands):
            obj = obj.select(ee.List(target_bands))
        else:
            obj = obj.select(model.band_of_use)

        # mosaic over the collectiontion (== convert imagecollection to image)
        if force_mosaic is None:
            # use attribute of the model
            if model._is_mosaic:
                obj = obj.mosaic()
        elif force_mosaic:
            obj = obj.mosaic()
        else:
            obj = obj

    return obj


# =============================================================================
# Helpers
# =============================================================================


def _is_image(geeobj):
    return isinstance(geeobj, ee.image.Image)


def _is_imagecollection(geeobj):
    return isinstance(geeobj, ee.imagecollection.ImageCollection)


def _validate_metadf(metadf):
    """Test if metadf is valid for GEE extraction.

    Returns True if metadata is suitable for gee extraction.

    :param metadf: metadata dataframe
    :type metadf: pd.DataFrame
    :return: True if oke, else False
    :rtype: Bool

    """
    if metadf.empty:
        return False
    if metadf["geometry"].x.isnull().values.all():
        return False
    if metadf["geometry"].y.isnull().values.all():
        return False
    try:
        # Just testing if it can be converted
        metadf = metadf.to_crs("epsg:4326")
    except:
        return False

    return True


def _addDate(image):
    """Add the image datetime as a band."""
    img_date = ee.Date(image.date())
    img_date = ee.Number.parse(img_date.format("YYYYMMddHHmmss"))
    return image.addBands(ee.Image(img_date).rename("datetime"))


def _df_to_features_point_collection(df):
    """Convert a dataframe to a featurecollections row-wise."""
    features = []
    for index, row in df.reset_index().iterrows():
        #     construct the geometry from dataframe
        poi_geometry = ee.Geometry.Point([row["lon"], row["lat"]])
        #     construct the attributes (properties) for each point
        poi_properties = poi_properties = {"name": row["name"]}
        #     construct feature combining geometry and properties
        poi_feature = ee.Feature(poi_geometry, poi_properties)
        features.append(poi_feature)

    return ee.FeatureCollection(features)


def _df_to_features_buffer_collection(df, bufferradius):
    """Convert a dataframe to a featurecollections row-wise."""
    features = []
    for index, row in df.reset_index().iterrows():
        #     construct the geometry from dataframe
        poi_geometry = ee.Geometry.Point([row["lon"], row["lat"]]).buffer(
            distance=bufferradius
        )
        #     construct the attributes (properties) for each point
        poi_properties = poi_properties = {"name": row["name"]}
        #     construct feature combining geometry and properties
        poi_feature = ee.Feature(poi_geometry, poi_properties)
        features.append(poi_feature)

    return ee.FeatureCollection(features)


def coordinates_available(metadf, latcol="lat", loncol="lon"):
    """Test if all coordinates are available."""
    if metadf[latcol].isnull().all():
        logger.warning("No coordinates are found!")
        return False
    if metadf[loncol].isnull().all():
        logger.warning("No coordinates are found!")
        return False
    return True


def _estimate_data_size(metadf, startdt, enddt, time_res, n_bands=1):
    datatimerange = pd.date_range(start=startdt, end=enddt, freq=time_res)

    return metadf.shape[0] * len(datatimerange) * n_bands


# =============================================================================
# Docstring test
# =============================================================================
if __name__ == "__main__":
    from metobs_toolkit.doctest_fmt import setup_and_run_doctest

    setup_and_run_doctest()
