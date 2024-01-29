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


# =============================================================================
# Object convertors
# =============================================================================


def _datetime_to_gee_datetime(datetime):
    # covert to UTC!
    # utcdt = datetime.astimezone(pytz.utc) #this will assume datetime in lt !!!
    # logger.debug(utcdt.replace(tzinfo=None))
    return ee.Date(datetime.replace(tzinfo=None))


def get_ee_obj(trg_location, is_image, is_imagecollection, band=None):
    """Get an image from a GEE object."""
    if is_image:
        obj = ee.Image(trg_location)
    elif is_imagecollection:
        if isinstance(band, type(None)):
            obj = ee.ImageCollection(trg_location)
        else:
            obj = ee.ImageCollection(trg_location).select(band)

    else:
        sys.exit("Map type is not an Image or Imagecollection.")
    return obj


def coords_to_geometry(lat=[], lon=[], proj="EPSG:4326"):
    """Convert coordinates to GEE geometries."""
    if len(lat) == 1:
        return ee.Geometry.Point(coords=[lon[0], lat[0]], proj=proj)
    else:
        return ee.Geometry.MultiPoint(list(zip(lon, lat)), proj=proj)


# =============================================================================
# Helpers
# =============================================================================


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
    for index, row in df.iterrows():
        #     construct the geometry from dataframe
        poi_geometry = ee.Geometry.Point([row["lon"], row["lat"]])
        #     construct the attributes (properties) for each point
        poi_properties = {"feature_idx": ee.String(index)}
        #     construct feature combining geometry and properties
        poi_feature = ee.Feature(poi_geometry, poi_properties)
        features.append(poi_feature)

    return ee.FeatureCollection(features)


def _df_to_features_buffer_collection(df, bufferradius):
    """Convert a dataframe to a featurecollections row-wise."""
    features = []
    for index, row in df.iterrows():
        #     construct the geometry from dataframe
        poi_geometry = ee.Geometry.Point([row["lon"], row["lat"]]).buffer(
            distance=bufferradius
        )
        #     construct the attributes (properties) for each point
        poi_properties = poi_properties = {"feature_idx": index}
        #     construct feature combining geometry and properties
        poi_feature = ee.Feature(poi_geometry, poi_properties)
        features.append(poi_feature)

    return ee.FeatureCollection(features)


def _estimate_data_size(metadf, startdt, enddt, time_res, n_bands=1):
    datatimerange = pd.date_range(start=startdt, end=enddt, freq=time_res)

    return metadf.shape[0] * len(datatimerange) * n_bands


# =============================================================================
# Data extractors
# =============================================================================


def extract_pointvalues(
    metadf, scale, trg_gee_loc, band_of_use, is_imagecollection, is_image
):
    """
    TODO: update this docstring
    Extract values for point locations from a GEE dataset.

    The pointlocations are defined in a dataframe by EPSG:4326 lat lon coordinates.

    A dataframe with the extracted values is returned.
    The values are mapped to human classes if the dataset value type is labeld as categorical.

    Parameters
    ----------
    metadf : pd.DataFrame
        dataframe containing coordinates and a column "name", representing the name for each location.
    mapinfo : Dict
        The information about the GEE dataset.
    output_column_name : String
        Column name for the extracted values.
    latcolname : String, optional
        Columnname of latitude values. The default is 'lat'.
    loncolname : String, optional
        Columnname of longitude values. The default is 'lon'.

    Returns
    -------
    pd.DataFrame
        A dataframe with name as index, all columns from the metadf + extracted extracted values column.

    """

    # =============================================================================
    # df to featurecollection
    # =============================================================================

    ee_fc = _df_to_features_point_collection(metadf)

    # =============================================================================
    # extract raster values
    # =============================================================================
    raster = get_ee_obj(
        trg_location=trg_gee_loc,
        is_image=is_image,
        is_imagecollection=is_imagecollection,
        band=band_of_use,
    )

    if is_imagecollection:

        def rasterExtraction(image):
            feature = image.sampleRegions(
                collection=ee_fc,  # feature collection here
                scale=scale,  # Cell size of raster
            )
            return feature

        results = raster.map(rasterExtraction).flatten().getInfo()
    elif is_image:
        # raster = get_ee_obj(mapinfo, mapinfo["band_of_use"])  # dataset
        results = raster.sampleRegions(
            collection=ee_fc, scale=scale  # feature collection here
        ).getInfo()
    else:
        sys.exit("gee dataset is neighter image nor imagecollection.")

    # extract properties
    if not bool(results["features"]):
        # no data retrieved
        logger.warning(f"Something went wrong, gee did not return any data: {results}")
        logger.info(
            f"(Could it be that (one) these coordinates are not on the map: {metadf}?)"
        )
        return pd.DataFrame()

    # =============================================================================
    # to dataframe
    # =============================================================================

    properties = [x["properties"] for x in results["features"]]
    df = pd.DataFrame(properties)

    return df


def extract_buffer_frequencies(
    metadf, scale, trg_gee_loc, band_of_use, is_imagecollection, is_image, bufferradius
):
    """
    update docstring
    Extract buffer fractions from a GEE categorical dataset.

    The pointlocations are defined in a dataframe by EPSG:4326 lat lon coordinates.

    A dataframe with the extracted values is returned.
    The values are mapped to human classes if the dataset value type is labeld as categorical.

    Parameters
    ----------
    metadf : pd.DataFrame
        dataframe containing coordinates and a column "name", representing the name for each location.
    mapinfo : Dict
        The information about the GEE dataset.
    latcolname : String, optional
        Columnname of latitude values. The default is 'lat'.
    loncolname : String, optional
        Columnname of longitude values. The default is 'lon'.

    Returns
    -------
    pd.DataFrame
        A dataframe with name as index, all columns from the metadf + extracted extracted values column.

    """

    # =============================================================================
    # df to featurecollection
    # =============================================================================

    ee_fc = _df_to_features_buffer_collection(metadf, bufferradius)

    # =============================================================================
    # extract raster values
    # =============================================================================
    raster = get_ee_obj(
        trg_location=trg_gee_loc,
        is_image=is_image,
        is_imagecollection=is_imagecollection,
        band=band_of_use,
    )

    def rasterExtraction(image):
        feature = image.reduceRegions(
            reducer=ee.Reducer.frequencyHistogram(),
            collection=ee_fc,  # feature collection here
            scale=scale,  # Cell size of raster
        )
        return feature

    results = raster.map(rasterExtraction).flatten().getInfo()

    # =============================================================================
    # to dataframe
    # =============================================================================

    freqs = {
        staprop["properties"]["feature_idx"]: staprop["properties"]["histogram"]
        for staprop in results["features"]
    }
    freqsdf = pd.DataFrame(freqs)

    # format frequency df
    freqsdf = freqsdf.transpose().fillna(0)
    freqsdf.index.name = "feature_idx"

    # normalize freqs
    freqsdf = freqsdf.div(freqsdf.sum(axis=1), axis=0)

    return freqsdf


def gee_extract_timeseries(
    metadf,
    bandnames,
    startdt,
    enddt,
    scale,
    timeres,
    trg_gee_loc,
    is_imagecollection,
    is_image,
    gdrive_filename,
):
    """TODO: update this docstring
    Extract timeseries data at the stations location from a GEE dataset.

    Extract a timeseries, for a given obstype, for point locations from a GEE
    dataset. The pointlocations are defined in a dataframe by EPSG:4326 lat lon
    coordinates.

    The startdate is included, the enddate is excluded.

    A multi-index dataframe with the timeseries is returned

    Parameters
    ----------
    metadf : pd.DataFrame
        dataframe containing coordinates and a column "name", representing the name for each location.
    band_mapper : dict
        the name of the band to extract data from as keys, the default name of
        the corresponding obstype as values.
    mapinfo : Dict
        The information about the GEE dataset.
    startdt : datetime obj
        Start datetime for timeseries (included).
    enddt : datetime obj
        End datetime for timeseries (excluded).
    latcolname : String, optional
        Columnname of latitude values. The default is 'lat'.
    loncolname : String, optional
        Columnname of longitude values. The default is 'lon'.

    Returns
    -------
    pd.DataFrame
        A dataframe with name - datetime multiindex, all columns from the metadf + extracted timeseries
        column with the same name as the obstypes.

    """

    use_drive = False
    _est_data_size = _estimate_data_size(
        metadf=metadf,
        startdt=startdt,
        enddt=enddt,
        time_res=timeres,
        n_bands=len(bandnames),
    )

    if _est_data_size > 4000:
        print(
            "THE DATA AMOUT IS TO LAREGE FOR INTERACTIVE SESSION, THE DATA WILL BE EXPORTED TO YOUR GOOGLE DRIVE!"
        )
        logger.info(
            "THE DATA AMOUT IS TO LAREGE FOR INTERACTIVE SESSION, THE DATA WILL BE EXPORTED TO YOUR GOOGLE DRIVE!"
        )

        use_drive = True
    # =============================================================================
    # df to featurecollection
    # =============================================================================

    ee_fc = _df_to_features_point_collection(metadf)

    # =============================================================================
    # extract raster values
    # =============================================================================

    def rasterExtraction(image):
        feature = image.sampleRegions(
            collection=ee_fc,  # feature collection here
            scale=scale,  # Cell size of raster
        )
        return feature

    # Because the daterange is maxdate exclusive, add the time resolution to the enddt
    enddt = enddt + pd.Timedelta(timeres)

    # construct raster obj
    raster = get_ee_obj(
        trg_location=trg_gee_loc,
        is_image=is_image,
        is_imagecollection=is_imagecollection,
        band=bandnames,
    )

    # filter out timeseries
    results = (
        raster.filter(
            ee.Filter.date(
                _datetime_to_gee_datetime(startdt), _datetime_to_gee_datetime(enddt)
            )
        )
        .map(_addDate)
        .map(rasterExtraction)
        .flatten()
    )

    if not use_drive:
        results = results.getInfo()

        # =============================================================================
        # to dataframe
        # =============================================================================

        # extract properties
        properties = [x["properties"] for x in results["features"]]
        df = pd.DataFrame(properties)

        if df.empty:
            sys.exit("ERROR: the returned timeseries from GEE are empty.")

        # df = format_df(df, band_mapper)
        return df

    else:
        _filename = gdrive_filename
        _drivefolder = "era5_timeseries"

        print(
            f"The timeseries will be writen to your Drive in {_drivefolder}/{_filename} "
        )
        logger.info(
            f"The timeseries will be writen to your Drive in {_drivefolder}/{_filename} "
        )

        data_columns = ["datetime", "feature_idx"]
        data_columns.extend(bandnames)

        task = ee.batch.Export.table.toDrive(
            collection=results,
            description="extracting_era5",
            folder=_drivefolder,
            fileNamePrefix=str(gdrive_filename),
            fileFormat="CSV",
            selectors=data_columns,
        )

        task.start()
        logger.info("The google server is handling your request ...")
        sleep(3)
        finished = False
        while finished is False:
            if task.status()["state"] == "READY":
                logger.info("Awaitening execution ...")
                sleep(4)
            elif task.status()["state"] == "RUNNING":
                logger.info("Running ...")
                sleep(4)
            else:
                logger.info("finished")
                finished = True

        doc_folder_id = task.status()["destination_uris"][0]
        print("The data is transfered! Open the following link in your browser: \n\n")
        print(f"{doc_folder_id} \n\n")
        print(
            "To upload the data to the model, use the Modeldata.set_model_from_csv() method"
        )

        return init_multiindexdf()
