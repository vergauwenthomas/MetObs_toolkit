import logging
import copy
from typing import Union

from numpy import nan, float32
import pandas as pd
import geopandas as gpd

from metobs_toolkit.backend_collection.df_helpers import convert_to_numeric_series
from metobs_toolkit.geedatasetmanagers import GEEStaticDatasetManager, global_LCZ_map
from metobs_toolkit.gee_api import connect_to_gee
from metobs_toolkit.backend_collection.errorclasses import (
    MetObsAdditionError,
)
import metobs_toolkit.backend_collection.printing_collection as printing


from metobs_toolkit.backend_collection.loggingmodule import log_entry

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

        # Common/special physical properties
        self._LCZ = nan  # set by the setup
        self._altitude = nan  # set by the setup

        # additional data
        self._extradata = dict(extradata)

        # Data extracted from other sources
        self._geedata = {}  # example: "LCZ": 'LCZ-4'
        self._gee_buffered_fractions = (
            {}
        )  # example: {100: {pervious: 0.8, impervious: 0.2}}

        self.setup()

    def setup(self):
        # transfer altitude to attr
        if "altitude" in self.extradata.keys():
            logger.debug(
                f'Setting altitude({float(self.extradata["altitude"])}) attribute for site {self.stationname} from extradata'
            )
            self._altitude = float(self.extradata["altitude"])
        else:
            self._altitude = nan

        # transfer LCZ to attr
        if "LCZ" in self.extradata.keys():
            logger.debug(
                f"Setting LCZ attribute for site {self.stationname} from extradata"
            )
            self._LCZ = str(self.extradata["LCZ"])
        else:
            self._LCZ = nan

    def _id(self) -> str:
        """A physical unique id.

        In the __add__ methods, if the id of two instances differs, adding is
        a regular concatenation.
        """
        return f"{self.stationname}"

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
        alt_equal = (self.altitude == other.altitude) or (
            pd.isnull(self.altitude) and pd.isnull(other.altitude)
        )
        lcz_equal = (self.LCZ == other.LCZ) or (
            pd.isnull(self.LCZ) and pd.isnull(other.LCZ)
        )

        return (
            self.stationname == other.stationname
            and lat_equal
            and lon_equal
            and alt_equal
            and lcz_equal
            and self.extradata == other.extradata
            and self._geedata == other._geedata
            and self._gee_buffered_fractions == other._gee_buffered_fractions
        )

    def __add__(self, other: "Site") -> "Site":
        # We assume an outside merge on the same name, and same coordinates

        # lat, lon and name must be the same!
        lat_equal = (self.lat == other.lat) or (
            pd.isnull(self.lat) and pd.isnull(other.lat)
        )
        lon_equal = (self.lon == other.lon) or (
            pd.isnull(self.lon) and pd.isnull(other.lon)
        )

        if not (self._id() == other._id() and lat_equal and lon_equal):
            raise MetObsAdditionError(
                f"Could not sum {self} and {other}, since the name or coordinates are not equal."
            )

        new_alt = next(
            (x for x in [self.altitude, other.altitude] if not pd.isna(x)), nan
        )
        new_LCZ = next((x for x in [self.LCZ, other.LCZ] if not pd.isna(x)), nan)

        # Update metadata (dict-style)
        new_extradata = copy.copy(self.extradata)
        new_extradata.update(other.extradata)  # dictupdate

        new_geedata = copy.copy(self._geedata)
        new_geedata.update(other._geedata)

        new_gee_buffered = copy.copy(self._gee_buffered_fractions)
        new_gee_buffered.update(other._gee_buffered_fractions)

        # Create a new site
        newsite = Site(
            stationname=self.stationname,
            latitude=self.lat,
            longitude=self.lon,
            extradata=new_extradata,
        )
        newsite._geedata = new_geedata
        newsite._gee_buffered_fractions = new_gee_buffered

        newsite.set_altitude(new_alt)
        newsite.set_LCZ(new_LCZ)
        return newsite

    def __repr__(self):
        """Return a string representation for debugging."""
        return f"{type(self).__name__}(id={self._id()})"

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
    def altitude(self) -> Union[float, type(nan)]:
        """Return the altitude."""
        return self._altitude

    @property
    def LCZ(self) -> Union[str, type(nan)]:
        """Return the LCZ."""
        return self._LCZ

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
                "altitude": self.altitude,
                "LCZ": self.LCZ,
                **self._geedata,  # unfold all gee extracted data
                **self.extradata,
            },  # unfold all extra data
            index=pd.Index(data=[self.stationname], name="name"),
        )

        # Ensure lat, lon, and altitude columns are float64
        for col in ["lat", "lon", "altitude"]:
            if col in metadf.columns:
                metadf[col] = convert_to_numeric_series(
                    metadf[col], datadtype=float32
                ).values

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

    @log_entry
    def set_altitude(self, altitude: Union[float, type(nan)]) -> None:
        """
        Set the altitude attribute.

        Parameters
        ----------
        altitude : float
            Altitude value to set.
        """
        self._altitude = float(altitude)

    @log_entry
    def set_LCZ(self, LCZ: Union[str, type(nan)]) -> None:
        """
        Set the LCZ attribute.

        Parameters
        ----------
        LCZ : str
            LCZ value to set.
        """
        # test if label is compatible with LCZ from global_LCZ_map
        if pd.isna(LCZ):
            self._LCZ = nan
        else:
            all_possible_labels = list(global_LCZ_map.class_map.values())
            if LCZ not in all_possible_labels:
                raise ValueError(
                    f"LCZ should be one of {all_possible_labels}, not {LCZ}"
                )
            self._LCZ = str(LCZ)

    @log_entry
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
        self._geedata[dataname] = value

    @log_entry
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
        # Replace NaNs in data with 0
        data = {k: (0 if pd.isna(v) else v) for k, v in data.items()}
        self._gee_buffered_fractions.update({buffer: data})

    @log_entry
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
        if not isinstance(metadata, dict):
            raise ValueError(f"metadata should be a dictionary, not {type(metadata)}")
        logger.debug(f"Adding metadata: {metadata}")
        self._extradata.update(dict(metadata))
        logger.info(f"Updated metadata: {self.extradata}")

    # ------------------------------------------
    #    Flaggers
    # ------------------------------------------

    @log_entry
    def flag_has_altitude(self) -> bool:
        """
        Check if the site has altitude information.

        Returns
        -------
        bool
            True if altitude is available, False otherwise.
        """
        return not pd.isnull(self.altitude)

    @log_entry
    def flag_has_LCZ(self) -> bool:
        """
        Check if the site has LCZ information.

        Returns
        -------
        bool
            True if LCZ is available, False otherwise.
        """
        return not pd.isnull(self.LCZ)

    @log_entry
    def flag_has_landcoverfractions(self) -> bool:
        """
        Check if the site has land cover fractions.

        Returns
        -------
        bool
            True if land cover fractions are available, False otherwise.
        """
        return bool(self.buffered_fractions)

    @log_entry
    def flag_has_coordinates(self) -> bool:
        """
        Check if the site has valid coordinates.

        Returns
        -------
        bool
            True if both latitude and longitude are available, False otherwise.
        """
        return (not pd.isnull(self.lat)) and (not pd.isnull(self.lon))

    # ------------------------------------------
    #    Methods
    # ------------------------------------------

    @log_entry
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
        if geedf.empty:
            logger.warning(
                f"No data returned by GEE when point extraction on {self} for {geestaticdataset.name}"
            )
            return nan

        return geedf[geestaticdataset.name].iloc[0]

    @log_entry
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

    def _get_info_core(self, nident_root=1) -> str:
        """
        Generate a formatted string containing metadata information for parent objects.
        This method compiles various metadata details such as coordinates, altitude,
        land cover zone (LCZ), land cover fractions, and extra metadata into a formatted
        string. It is primarily used by the `get_info` methods of parent objects.
        Parameters
        ----------
        nident_root : int, optional
            The base indentation level for the formatted output, by default 1.
        Returns
        -------
        str
            A formatted string containing the metadata information.
        """

        infostr = ""
        # Coordinates
        if self.flag_has_coordinates():
            infostr += printing.print_fmt_line(
                f"Coordinates ({self.lat}, {self.lon}) (latitude, longitude)",
                nident_root,
            )
        else:
            infostr += printing.print_fmt_line("Coordinates are unknown", nident_root)

        # Altitude
        if self.flag_has_altitude():
            infostr += printing.print_fmt_line(
                f"Altitude: {self.altitude} (m)", nident_root
            )
        else:
            infostr += printing.print_fmt_line("Altitude is unknown", nident_root)

        # LCZ
        if self.flag_has_LCZ():
            infostr += printing.print_fmt_line(f"LCZ: {self.LCZ}", nident_root)
        else:
            infostr += printing.print_fmt_line("LCZ is unknown", nident_root)

        # Buffered fractions
        if self.flag_has_landcoverfractions():
            infostr += printing.print_fmt_line(
                "Land cover fractions are available", nident_root
            )
        else:
            infostr += printing.print_fmt_line(
                "Land cover fractions are unknown", nident_root
            )

        # Extra metadata
        if bool(self):
            infostr += printing.print_fmt_line(
                "Extra metadata from the metadata file:", nident_root
            )
            for key, value in self.extradata.items():
                infostr += printing.print_fmt_line(f"{key}: {value}", nident_root + 1)
        else:
            infostr += printing.print_fmt_line(
                "No extra metadata available", nident_root
            )

        return infostr

    @log_entry
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
        infostr = ""
        infostr += printing.print_fmt_title("General Info of Site")
        infostr += printing.print_fmt_line(f"Site of {self.stationname}:", 0)
        infostr += self._get_info_core(nident_root=1)

        if printout:
            print(infostr)
        else:
            return infostr
