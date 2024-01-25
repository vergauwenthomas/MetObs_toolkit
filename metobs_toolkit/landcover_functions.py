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
#  Connection functions
# =============================================================================


def connect_to_gee():
    """Authenticate to GEE if needed."""
    if not ee.data._credentials:  # check if ee connection is initialized
        ee.Authenticate()
        ee.Initialize()
    return


# =============================================================================
# Top level functions (can be called by dataset)
# =============================================================================


def lcz_extractor(metadf, mapinfo):
    """Extract LCZ for all stations in the metadf."""
    # make return in case something went wrong
    default_return = pd.Series(
        index=metadf.index, data="Location_unknown", name="lcz", dtype=object
    )
    # test if metadata is suitable
    if not _validate_metadf(metadf):
        logger.warning(f"Metadf is not suitable for GEE extractiond: {metadf}")
        return default_return

    relevant_metadf = metadf.reset_index()[["name", "lat", "lon"]]

    lcz_df = extract_pointvalues(
        metadf=relevant_metadf, mapinfo=mapinfo, output_column_name="lcz"
    )
    if lcz_df.empty:
        return pd.Series(dtype=object)
    return lcz_df["lcz"]  # return series


def lc_fractions_extractor(metadf, mapinfo, buffer, agg):
    """Get landcover fractions for all buffers from GEE."""
    # make return in case something went wrong
    default_return = (pd.DataFrame(index=metadf.index), buffer)

    # test if metadata is suitable
    if not _validate_metadf(metadf):
        logger.warning(f"Metadf is not suitable for GEE extractiond: {metadf}")
        return default_return

    relevant_metadf = metadf.reset_index()[["name", "lat", "lon"]]

    freqs_df = extract_buffer_frequencies(
        metadf=relevant_metadf, mapinfo=mapinfo, bufferradius=buffer
    )

    # apply aggregation if required
    if agg:
        logger.info(f"Using aggregation scheme: {mapinfo['aggregation']}")
        agg_df = pd.DataFrame()
        for agg_name, agg_classes in mapinfo["aggregation"].items():
            present_agg_classes = [
                str(num) for num in agg_classes if str(num) in freqs_df.columns
            ]
            agg_df[agg_name] = freqs_df[present_agg_classes].sum(axis=1)

        return agg_df, buffer

    else:
        # map numeric classes to human
        mapper = {
            str(num): human for num, human in mapinfo["categorical_mapper"].items()
        }
        freqs_df = freqs_df.rename(columns=mapper)

        return freqs_df, buffer


def height_extractor(metadf, mapinfo):
    """Get altitude for all stations from GEE."""
    # make return in case something went wrong
    default_return = pd.Series(
        index=metadf.index, data="Location_unknown", name="altitude", dtype=object
    )

    # test if metadata is suitable
    if not _validate_metadf(metadf):
        logger.warning(f"Metadf is not suitable for GEE extractiond: {metadf}")
        return default_return

    relevant_metadf = metadf.reset_index()[["name", "lat", "lon"]]

    altitude_df = extract_pointvalues(
        metadf=relevant_metadf, mapinfo=mapinfo, output_column_name="altitude"
    )
    return altitude_df["altitude"]  # return series


# =============================================================================
# Object convertors
# =============================================================================


def _datetime_to_gee_datetime(datetime):
    # covert to UTC!
    utcdt = datetime.astimezone(pytz.utc)
    logger.debug(utcdt.replace(tzinfo=None))
    return ee.Date(utcdt.replace(tzinfo=None))


def get_ee_obj(mapinfo, band=None):
    """Get an image from a GEE object."""
    if mapinfo["is_image"]:
        obj = ee.Image(mapinfo["location"])
    elif mapinfo["is_imagecollection"]:
        if isinstance(band, type(None)):
            obj = ee.ImageCollection(mapinfo["location"])
        else:
            obj = ee.ImageCollection(mapinfo["location"]).select(band)

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
# Data extractors
# =============================================================================


def extract_pointvalues(metadf, mapinfo, output_column_name):
    """Extract values for point locations from a GEE dataset.

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
    scale = mapinfo["scale"]

    # test if coordiantes are available
    if not coordinates_available(metadf, "lat", "lon"):
        return pd.DataFrame()

    # =============================================================================
    # df to featurecollection
    # =============================================================================

    ee_fc = _df_to_features_point_collection(metadf)

    # =============================================================================
    # extract raster values
    # =============================================================================

    raster = get_ee_obj(mapinfo, mapinfo["band_of_use"])  # dataset
    if mapinfo["is_imagecollection"]:

        def rasterExtraction(image):
            feature = image.sampleRegions(
                collection=ee_fc,  # feature collection here
                scale=scale,  # Cell size of raster
            )
            return feature

        results = raster.map(rasterExtraction).flatten().getInfo()
    elif mapinfo["is_image"]:
        raster = get_ee_obj(mapinfo, mapinfo["band_of_use"])  # dataset
        results = raster.sampleRegions(
            collection=ee_fc, scale=scale  # feature collection here
        ).getInfo()
    else:
        sys.exit(
            f'gee dataset {mapinfo["location"]} is neighter image nor imagecollection.'
        )

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

    # map to human space if categorical
    if mapinfo["value_type"] == "categorical":
        df[mapinfo["band_of_use"]] = df[mapinfo["band_of_use"]].map(
            mapinfo["categorical_mapper"]
        )

    # rename to values to toolkit space
    df = df.rename(columns={mapinfo["band_of_use"]: output_column_name})

    # #format index
    df = df.set_index(["name"])

    return df


def extract_buffer_frequencies(metadf, mapinfo, bufferradius):
    """Extract buffer fractions from a GEE categorical dataset.

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
    scale = mapinfo["scale"]

    # test if coordiantes are available
    if not coordinates_available(metadf, "lat", "lon"):
        return pd.DataFrame()

    # test if map is categorical
    if not mapinfo["value_type"] == "categorical":
        logger.warning(
            "Extract buffer frequencies is only implemented for categorical datasets!"
        )
        return pd.DataFrame()

    # =============================================================================
    # df to featurecollection
    # =============================================================================

    ee_fc = _df_to_features_buffer_collection(metadf, bufferradius)

    # =============================================================================
    # extract raster values
    # =============================================================================

    def rasterExtraction(image):
        feature = image.reduceRegions(
            reducer=ee.Reducer.frequencyHistogram(),
            collection=ee_fc,  # feature collection here
            scale=scale,  # Cell size of raster
        )
        return feature

    raster = get_ee_obj(mapinfo, mapinfo["band_of_use"])  # dataset
    results = raster.map(rasterExtraction).flatten().getInfo()

    # =============================================================================
    # to dataframe
    # =============================================================================

    freqs = {
        staprop["properties"]["name"]: staprop["properties"]["histogram"]
        for staprop in results["features"]
    }
    freqsdf = pd.DataFrame(freqs)

    # format frequency df
    freqsdf = freqsdf.transpose().fillna(0)
    freqsdf.index.name = "name"

    # normalize freqs
    freqsdf = freqsdf.div(freqsdf.sum(axis=1), axis=0)

    return freqsdf


def gee_extract_timeseries(
    metadf, band_mapper, mapinfo, startdt, enddt, latcolname="lat", loncolname="lon"
):
    """Extract timeseries data at the stations location from a GEE dataset.

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
    scale = mapinfo["scale"]
    bandnames = list(band_mapper.keys())

    # test if coordiantes are available
    if not coordinates_available(metadf, latcolname, loncolname):
        return pd.DataFrame()

    use_drive = False
    _est_data_size = _estimate_data_size(
        metadf=metadf,
        startdt=startdt,
        enddt=enddt,
        time_res=mapinfo["time_res"],
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
    enddt = enddt + pd.Timedelta(mapinfo["time_res"])

    raster = get_ee_obj(mapinfo, bandnames)  # dataset
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

    def format_df(df, band_mapper):
        # format datetime
        df["datetime"] = pd.to_datetime(df["datetime"], format="%Y%m%d%H%M%S")
        # set timezone
        df["datetime"] = df["datetime"].dt.tz_localize("UTC")

        # format index
        df = df.set_index(["name", "datetime"])
        df = df.sort_index()

        # rename to values to toolkit space
        df = df.rename(columns=band_mapper)
        return df

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

        df = format_df(df, band_mapper)
        return df

    else:
        _filename = "era5_data"
        _drivefolder = "era5_timeseries"

        print(
            f"The timeseries will be writen to your Drive in {_drivefolder}/{_filename} "
        )
        logger.info(
            f"The timeseries will be writen to your Drive in {_drivefolder}/{_filename} "
        )

        data_columns = ["datetime", "name"]
        data_columns.extend(bandnames)

        task = ee.batch.Export.table.toDrive(
            collection=results,
            description="extracting_era5",
            folder=_drivefolder,
            fileNamePrefix=_filename,
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
