import logging
from typing import Union

import numpy as np
import pandas as pd
import geopandas as gpd

from metobs_toolkit.geedatasetmanagers import GEEStaticDatasetManager
from metobs_toolkit.gee_api import connect_to_gee
from metobs_toolkit.geedatasetmanagers import default_datasets as default_gee_datasets

logger = logging.getLogger("<metobs_toolkit>")


class Site:
    """
    Holds all metadata and location-based data of a site (station location).

    Parameters
    ----------
    stationname : str
        Name of the station.
    latitude : float
        Latitude of the station.
    longitude : float
        Longitude of the station.
    extradata : dict, optional
        Additional metadata for the site, by default {}.
    """

    def __init__(
        self, stationname: str, latitude: float, longitude: float, extradata: dict = {}
    ):
    
        # Set data
        self._stationname = stationname
        self._lat = float(latitude)
        self._lon = float(longitude)

        # additional data
        self._extradata = dict(extradata)

        # Data extracted from other sources
        self._geedata = {}  # example: "LCZ": 'LCZ-4'
        self._gee_buffered_fractions = {}  # example: {100: {pervious: 0.8, impervious: 0.2}}

    def __eq__(self, other):
        """Check equality with another Site object."""
        if not isinstance(other, Site):
            return False

        lat_equal = (self.lat == other.lat) or (
            pd.isnull(self.lat) and pd.isnull(other.lat)
        )
        lon_equal = (self.lon == other.lon) or (
            pd.isnull(self.lon) and pd.isnull(other.lon)
        )

        return (
            self.stationname == other.stationname
            and lat_equal
            and lon_equal
            and self.extradata == other.extradata
            and self._geedata == other._geedata
            and self._gee_buffered_fractions == other._gee_buffered_fractions
        )

    @property
    def stationname(self) -> str:
        """Return the station name."""
        return str(self._stationname)

    @property
    def lat(self) -> float:
        """Return the latitude."""
        return self._lat

    @property
    def lon(self) -> float:
        """Return the longitude."""
        return self._lon

    @property
    def altitude(self):
        """Return the altitude if available, else NaN."""
        if "altitude" in self.extradata.keys():
            # if altitude info was extracted from the metadatafile
            return self.extradata["altitude"]
        elif "altitude" in self._geedata.keys():
            # if altitude info was extracted from the GEE API
            return self._geedata["altitude"]
        else:
            return np.nan

    @property
    def LCZ(self):
        """Return the LCZ if available, else NaN."""
        if "LCZ" in self.extradata.keys():
            # if altitude info was extracted from the metadatafile
            return self.extradata["LCZ"]
        elif "LCZ" in self._geedata.keys():
            # if altitude info was extracted from the GEE API
            return self._geedata["LCZ"]
        else:
            return np.nan

    @property
    def buffered_fractions(self):
        """Return the buffered fractions dictionary."""
        return self._gee_buffered_fractions

    @property
    def extradata(self):
        """Return the extra metadata dictionary."""
        return self._extradata

    @property
    def metadf(self) -> pd.DataFrame:
        """Return a DataFrame with all metadata and geometry."""
        metadf = pd.DataFrame(
            data={
                "lat": self.lat,
                "lon": self.lon,
                **self._geedata,  # unfold all gee extracted data
                **self.extradata,
            },  # unfold all extra data
            index=pd.Index(data=[self.stationname], name="name"),
        )

        # add buffered fractions
        for bufradius, fracdict in self._gee_buffered_fractions.items():
            for covername, fraction in fracdict.items():
                metadf[f"{covername}_frac_{bufradius}m"] = fraction

        # Create geometry column (geopandasdataframe)
        metadf = gpd.GeoDataFrame(
            metadf, geometry=gpd.points_from_xy(metadf["lon"], metadf["lat"])
        )
        # add CRS
        metadf = metadf.set_crs("WGS84")
        return metadf

    # ------------------------------------------
    #   Setters
    # ------------------------------------------

    def set_geedata(self, dataname: str, value: Union[str, float]) -> None:
        """
        Set a value in the geedata dictionary.

        Parameters
        ----------
        dataname : str
            Name of the data field.
        value : str or float
            Value to set.
        """
        logger.debug("Entering set_geedata for %s", self)
        self._geedata[dataname] = value

    def set_gee_buffered_frac_data(self, buffer: float, data: dict) -> None:
        """
        Set buffered fraction data for a given buffer radius.

        Parameters
        ----------
        buffer : float
            Buffer radius.
        data : dict
            Dictionary of cover fractions.
        """
        logger.debug("Entering set_gee_buffered_frac_data for %s", self)
        # Replace NaNs in data with 0
        data = {k: (0 if pd.isna(v) else v) for k, v in data.items()}
        self._gee_buffered_fractions.update({buffer: data})

    def add_metadata(self, metadata: dict) -> None:
        """
        Add metadata to the site.

        Parameters
        ----------
        metadata : dict
            Dictionary of metadata to add.

        Raises
        ------
        ValueError
            If metadata is not a dictionary.
        """
        logger.debug("Entering add_metadata for %s", self)
        if not isinstance(metadata, dict):
            raise ValueError(f"metadata should be a dictionary, not {type(metadata)}")
        logger.debug(f"Adding metadata: {metadata}")
        self._extradata.update(dict(metadata))
        logger.info(f"Updated metadata: {self.extradata}")

    # ------------------------------------------
    #    Flaggers
    # ------------------------------------------

    def flag_has_altitude(self) -> bool:
        """
        Check if the site has altitude information.

        Returns
        -------
        bool
            True if altitude is available, False otherwise.
        """
        logger.debug("Entering flag_has_altitude for %s", self)
        return not pd.isnull(self.altitude)

    def flag_altitude_from_gee(self) -> bool:
        """
        Check if altitude information is from GEE.

        Returns
        -------
        bool
            True if altitude is from GEE, False otherwise.
        """
        logger.debug("Entering flag_altitude_from_gee for %s", self)
        return "altitude" in self._geedata.keys()

    def flag_has_LCZ(self) -> bool:
        """
        Check if the site has LCZ information.

        Returns
        -------
        bool
            True if LCZ is available, False otherwise.
        """
        logger.debug("Entering flag_has_LCZ for %s", self)
        return not pd.isnull(self.LCZ)

    def flag_LCZ_from_gee(self) -> bool:
        """
        Check if LCZ information is from GEE.

        Returns
        -------
        bool
            True if LCZ is from GEE, False otherwise.
        """
        logger.debug("Entering flag_LCZ_from_gee for %s", self)
        return "LCZ" in self._geedata.keys()

    def flag_has_landcoverfractions(self) -> bool:
        """
        Check if the site has land cover fractions.

        Returns
        -------
        bool
            True if land cover fractions are available, False otherwise.
        """
        logger.debug("Entering flag_has_landcoverfractions for %s", self)
        return bool(self.landcover_fractions)

    def flag_has_coordinates(self) -> bool:
        """
        Check if the site has valid coordinates.

        Returns
        -------
        bool
            True if both latitude and longitude are available, False otherwise.
        """
        logger.debug("Entering flag_has_coordinates for %s", self)
        return (not pd.isnull(self.lat)) and (not pd.isnull(self.lon))

    # ------------------------------------------
    #    Methods
    # ------------------------------------------

    def get_gee_point_metadata(
        self, geestaticdataset: GEEStaticDatasetManager, initialize_gee: bool = True
    ) -> Union[str, float]:
        """
        Extract static point metadata from a GEEStaticDatasetManager.

        Parameters
        ----------
        geestaticdataset : GEEStaticDatasetManager
            The GEE static dataset manager.
        initialize_gee : bool, optional
            Whether to initialize the GEE API, by default True.

        Returns
        -------
        str or float
            The extracted metadata value.

        Raises
        ------
        ValueError
            If geestaticdataset is not a GEEStaticDatasetManager or coordinates are unknown.
        """
        logger.debug("Entering get_gee_point_metadata for %s", self)
        # test if modeldata is static
        if not isinstance(geestaticdataset, GEEStaticDatasetManager):
            raise ValueError(
                f"geestaticdataset should be an isntance of GeeStaticDataset, not {type(geestaticdataset)}"
            )

        # test if coordinates are knonw
        if not self.flag_has_coordinates():
            raise ValueError(f"Coordinates of {self.stationname} are unknown")

        # initialize gee api
        if initialize_gee:
            connect_to_gee()

        geedf = geestaticdataset.extract_static_point_data(self.metadf)
        return geedf[geestaticdataset.name].iloc[0]

    def get_gee_point_buffer_fractions(
        self,
        geestaticdataset: GEEStaticDatasetManager,
        buffers: list = [100.0],
        aggregate: bool = False,
        initialize_gee: bool = True,
    ) -> dict:
        """
        Extract buffer fractions from a GEEStaticDatasetManager.

        Parameters
        ----------
        geestaticdataset : GEEStaticDatasetManager
            The GEE static dataset manager.
        buffers : list, optional
            List of buffer radii, by default [100.0].
        aggregate : bool, optional
            Whether to aggregate fractions, by default False.
        initialize_gee : bool, optional
            Whether to initialize the GEE API, by default True.

        Returns
        -------
        dict
            Dictionary of buffer fractions.

        Raises
        ------
        ValueError
            If geestaticdataset is not a GEEStaticDatasetManager or coordinates are unknown.
        """
        logger.debug("Entering get_gee_point_buffer_fractions for %s", self)
        # test if modeldata is static
        if not isinstance(geestaticdataset, GEEStaticDatasetManager):
            raise ValueError(
                f"geestaticdataset should be an isntance of GeeStaticDataset, not {type(geestaticdataset)}"
            )

        # test if coordinates are knonw
        if not self.flag_has_coordinates():
            raise ValueError(f"Coordinates of {self.stationname} are unknown")

        # initialize gee api
        if initialize_gee:
            connect_to_gee()

        bufferdict = {}
        for bufferradius in buffers:
            geedf = geestaticdataset.extract_static_buffer_frac_data(
                metadf=self.metadf, bufferradius=bufferradius, agg_bool=aggregate
            )
            fracdict = (
                geedf.reset_index()
                .set_index("buffer_radius")
                .drop(columns="name")
                .to_dict(orient="index")
            )
            bufferdict.update(fracdict)
        return bufferdict

    def get_info(self, printout: bool = True) -> Union[str, None]:
        """
        Retrieve and optionally print basic information about the site.

        Parameters
        ----------
        printout : bool, optional
            If True, print the information. If False, return as string.

        Returns
        -------
        str or None
            Information string if printout is False, else None.
        """
        logger.debug("Entering get_info for %s", self)
        infostr = f"Site of {self.stationname}:\n"

        infostr += f" -- metadata from file --\n"
        if self.flag_has_coordinates():
            infostr += (
                f"  * Coordinates ({self.lat}, {self.lon}) (latitude, longitude)\n"
            )
        else:
            infostr += "  * Coordinates are unknown\n"

        if self.flag_has_altitude() & (not self.flag_altitude_from_gee()):
            infostr += f"  * Altitude: {self.altitude} (m)\n"

        infostr += "  Extra metadata from the metadatafile:\n"
        for key, value in self.extradata.items():
            infostr += f"    * {key}: {value}\n"

        infostr += f" -- data extracted from GEE --\n"
        if self.flag_has_altitude() & (self.flag_altitude_from_gee()):
            infostr += f"  * Altitude: {self.altitude} (m)\n"

        if self.flag_has_LCZ() & (self.flag_LCZ_from_gee()):
            infostr += f"  * LCZ: {self.LCZ}\n"

        for key, value in self._geedata.items():
            infostr += f"  * {key}: {value}\n"

        if printout:
            print(infostr)
        else:
            return infostr
