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
    """Holds all metadata and location based data of a site (station location)"""

    def __init__(
        self, stationname: str, latitude: float, longitude: float, extradata: dict = {}
    ):

        logger.debug(
            f"Initializing Site of {stationname} with latitude={latitude}, longitude={longitude}, extradata={extradata}"
        )
        # Set data
        self._stationname = stationname
        self._lat = float(latitude)
        self._lon = float(longitude)

        # additional data
        self._extradata = dict(extradata)

        # Data extracted from other sources
        self._geedata = {}  # example: "lcz": 'LCZ-4
        self._gee_buffered_fractions = {}  # example: {100: {pervious: 0.8,
        #                                                impervious: 0.2}}

    def __eq__(self, other):
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
    def stationname(self):
        return str(self._stationname)

    @property
    def lat(self):
        return self._lat

    @property
    def lon(self):
        return self._lon

    @property
    def altitude(self):
        if "altitude" in self.extradata.keys():
            # if altitude info was extracted from the metadatafile
            return self.extradata["altitude"]
        elif "altitude" in self._geedata.keys():
            # if altitude info was extracted from the GEE API
            return self._geedata["altitude"]
        else:
            return np.nan

    @property
    def lcz(self):
        if "lcz" in self.extradata.keys():
            # if altitude info was extracted from the metadatafile
            return self.extradata["lcz"]
        elif "lcz" in self._geedata.keys():
            # if altitude info was extracted from the GEE API
            return self._geedata["lcz"]
        else:
            return np.nan

    @property
    def buffered_fractions(self):
        return self._gee_buffered_fractions

    @property
    def extradata(self):
        return self._extradata

    @property
    def metadf(self) -> pd.DataFrame:

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

    def set_geedata(self, dataname: str, value: Union[str, float]):
        self._geedata[dataname] = value

    def set_gee_buffered_frac_data(self, buffer, data: dict):
        # Replace NaNs in data with 0
        data = {k: (0 if pd.isna(v) else v) for k, v in data.items()}
        self._gee_buffered_fractions.update({buffer: data})

    def add_metadata(self, metadata: dict):
        if not isinstance(metadata, dict):
            raise ValueError(f"metadata should be a dictionary, not {type(metadata)}")
        logger.debug(f"Adding metadata: {metadata}")
        self._extradata.update(dict(metadata))
        logger.info(f"Updated metadata: {self.extradata}")

    # ------------------------------------------
    #    Flaggers
    # ------------------------------------------

    def flag_has_altitude(self):
        return not pd.isnull(self.altitude)

    def flag_altitude_from_gee(self):
        return "altitude" in self._geedata.keys()

    def flag_has_lcz(self):
        return not pd.isnull(self.lcz)

    def flag_lcz_from_gee(self):
        return "lcz" in self._geedata.keys()

    def flag_has_landcoverfractions(self):
        return bool(self.landcover_fractions)

    def flag_has_coordinates(self):
        return (not pd.isnull(self.lat)) and (not pd.isnull(self.lon))

    # ------------------------------------------
    #    Methods
    # ------------------------------------------

    def get_gee_point_metadata(
        self, geestaticdataset, initialize_gee: bool = True
    ) -> Union[str, float]:
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
        geestaticdataset,
        buffers=[100.0],
        aggregate=False,
        initialize_gee: bool = True,
    ) -> dict:

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

        if self.flag_has_lcz() & (self.flag_lcz_from_gee()):
            infostr += f"  * LCZ: {self.lcz}\n"

        for key, value in self._geedata.items():
            infostr += f"  * {key}: {value}\n"

        if printout:
            print(infostr)
        else:
            return infostr
