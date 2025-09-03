import logging
import numpy as np
import pandas as pd

from metobs_toolkit.io_collection.filereaders import FileReader
from metobs_toolkit.template import Template
from metobs_toolkit.backend_collection.errorclasses import MetObsInconsistentStationName

from metobs_toolkit.backend_collection.loggingmodule import log_entry

logger = logging.getLogger("<metobs_toolkit>")


class MetaDataParser:
    """
    Parser class for metadata.

    This class handles the parsing, renaming, and formatting of metadata files
    according to a provided template.

    Parameters
    ----------
    metadatafilereader : FileReader
        The file reader object to read the metadata file.
    template : Template
        The template object that defines the metadata structure.
    """

    def __init__(self, metadatafilereader: FileReader, template: Template):
        """Initialize MetaDataParser."""
        self.filereader = metadatafilereader
        self.template = template

        # datadf with 'name' as index
        self.datadf = pd.DataFrame()  # metadata in formatted DataFrame style

    @log_entry
    def parse(self, **readkwargs) -> None:
        """
        Parse the metadata file and format it according to the template.

        Reads the raw metadata, applies blacklist renaming, sets the station name,
        renames columns to standard names, and subsets to mapped columns.

        Parameters
        ----------
        **readkwargs
            Additional keyword arguments passed to the file reader's read method.

        Returns
        -------
        None
        """
        # Read in the raw metadata
        raw_metadata = self.filereader.read(**readkwargs)
        self.template._metadata_template_compatibility_test(raw_metadata.columns)

        # 1. blacklist label handling
        metadf = self._rename_blacklabels(raw_metadata)

        # 2. name column handling (ID)
        metadf = self._set_name(metadf)

        # 3. Rename all columns to standards
        metadf = self._rename_raw_columns(metadf)

        # 4. Subset to mapped columns
        metadf = self._subset_to_mapped_columns(metadf)

        # set the metadf as data attribute
        self.datadf = metadf

    def _rename_blacklabels(self, rawdf: pd.DataFrame) -> pd.DataFrame:
        """
        Rename columns in the DataFrame using the blacklist mapping from the template.

        Parameters
        ----------
        rawdf : pd.DataFrame
            The raw metadata DataFrame.

        Returns
        -------
        pd.DataFrame
            DataFrame with columns renamed according to the blacklist.
        """
        blacklist_mapper = self.template._apply_blacklist(
            columns=rawdf.columns, on_data=False
        )
        rawdf.rename(columns=blacklist_mapper, inplace=True)
        return rawdf

    def _set_name(self, rawdf: pd.DataFrame) -> pd.DataFrame:
        """
        Add a 'name' column and set it as the index.

        Handles single-station and multi-station cases, ensures uniqueness,
        and sets the index to 'name'.

        Parameters
        ----------
        rawdf : pd.DataFrame
            The metadata DataFrame.

        Returns
        -------
        pd.DataFrame
            DataFrame with 'name' as the index.
        """
        if (self.template._is_data_single_station()) & (
            self.template.metadata_namemap["name"] is None
        ):
            rawdf["name"] = self.template.single_station_name
        else:
            rawdf.rename(columns=self.template._get_metadata_name_map(), inplace=True)

        # make sure the name column values are strings
        rawdf["name"] = rawdf["name"].astype(str)

        # test uniqueness
        if not rawdf["name"].is_unique:
            duplicates = rawdf[rawdf["name"].duplicated(keep=False)]
            logger.warning(
                f'Duplicate names found in metadata: {duplicates["name"].tolist()}'
            )
            rawdf = rawdf.drop_duplicates(subset="name", keep="first")

        # set index
        rawdf.set_index("name", inplace=True)

        return rawdf

    def _rename_raw_columns(self, rawdf: pd.DataFrame) -> pd.DataFrame:
        """
        Rename columns in the DataFrame to standard names as defined in the template.

        Parameters
        ----------
        rawdf : pd.DataFrame
            The metadata DataFrame.

        Returns
        -------
        pd.DataFrame
            DataFrame with columns renamed to standard names.
        """
        metacolmap = self.template._get_metadata_column_map()
        rawdf.rename(columns=metacolmap, inplace=True)
        return rawdf

    def _subset_to_mapped_columns(self, rawdf: pd.DataFrame) -> pd.DataFrame:
        """
        Subset the DataFrame to only include columns mapped in the template.

        Parameters
        ----------
        rawdf : pd.DataFrame
            The metadata DataFrame.

        Returns
        -------
        pd.DataFrame
            DataFrame containing only the relevant mapped columns.
        """
        logger.debug(
            f"Entering _subset_to_mapped_columns() of {self.__class__.__name__}"
        )
        metacolmap = self.template._get_metadata_column_map()
        relev_columns = list(metacolmap.values())

        # drop name --> that is the index
        relev_columns = [col for col in relev_columns if col != "name"]
        return rawdf[relev_columns]

    @log_entry
    def get_station_lon(self, stationname: str) -> float:
        """
        Get the longitude of a station from the metadata.

        Parameters
        ----------
        stationname : str
            The name of the station.

        Returns
        -------
        float
            Longitude of the station, or np.nan if not found.
        """
        if self._check_stationname_is_known(stationname=stationname):
            if "lon" in self.datadf.columns:
                return float(self.datadf.loc[stationname, "lon"])
            else:
                logger.warning(
                    f"No longitude is found for {stationname} in the metadata!"
                )
                return np.nan
        else:
            # station not in metadata
            return np.nan

    @log_entry
    def get_station_lat(self, stationname: str) -> float:
        """
        Get the latitude of a station from the metadata.

        Parameters
        ----------
        stationname : str
            The name of the station.

        Returns
        -------
        float
            Latitude of the station, or np.nan if not found.
        """
        if self._check_stationname_is_known(stationname=stationname):
            if "lat" in self.datadf.columns:
                return float(self.datadf.loc[stationname, "lat"])
            else:
                logger.warning(
                    f"No latitude is found for {stationname} in the metadata!"
                )
                return np.nan
        else:
            # station not in metadata
            return np.nan

    @log_entry
    def get_station_extra_metadata(self, stationname: str) -> dict:
        """
        Get extra metadata for a station, excluding latitude and longitude.

        Parameters
        ----------
        stationname : str
            The name of the station.

        Returns
        -------
        dict
            Dictionary of extra metadata for the station.
        """
        logger.debug(
            f"Entering get_station_extra_metadata() of {self.__class__.__name__}"
        )
        not_extra_columns = ["lat", "lon"]
        if self._check_stationname_is_known(stationname=stationname):
            extra_info = (
                self.datadf.loc[stationname]
                .drop(labels=not_extra_columns, errors="ignore")
                .to_dict()
            )
            return extra_info
        else:
            # station not in metadata
            return {}

    def _check_stationname_is_known(self, stationname: str) -> bool:
        """
        Check if the station name is present in the metadata.

        Parameters
        ----------
        stationname : str
            The name of the station.

        Returns
        -------
        bool
            True if the station name is known, False otherwise.
        """
        logger.debug(
            f"Entering _check_stationname_is_known() of {self.__class__.__name__}"
        )
        if stationname in self.datadf.index:
            return True
        else:
            logger.warning(f"{stationname} is not found in the metadata!")
            return False

    def _overwrite_name(self, target_single_name: str) -> None:
        """
        Overwrite the station name in single-station cases if needed.

        In single-station cases, the name can be defined in the template,
        but a (different) name can be present in the metadata file for the same station.
        This results in incompatible data-metadata.

        If the target name is present in the metadata, nothing is changed.
        If not present and multiple stations are in the metadata, an error is raised.
        If only one station is present, it is renamed to the target name.

        Parameters
        ----------
        target_single_name : str
            The target station name to enforce.

        Returns
        -------
        None

        Raises
        ------
        MetObsInconsistentStationName
            If multiple stations are present and the target name is not found.
        """
        # 1. target_single_name is included in the metadata
        if target_single_name in self.datadf.index:
            return
        else:
            # 2. not present in the metadata
            if self.datadf.shape[0] > 1:
                raise MetObsInconsistentStationName(
                    f"""
    The station name used in the single-station data is {target_single_name} (is a column in the data file,
    or defined in the template file). This station name is NOT present in the metadata, and multiple stations are
    present. Make sure that the single station name used by the data is the same
    as in the metadata file. (Or make sure that there is only one station present
    in the metadata file).
    """
                )
            else:
                logger.warning(
                    f"Due to a mismatch in single-station-name, the metadata name: {self.datadf.index[0]} --> {target_single_name} is renamed to be in line with the data file and template."
                )
                # rename
                self.datadf.index = pd.Index(name="name", data=[target_single_name])
                return

    # ------------------------------------------
    #    Getters
    # ------------------------------------------

    @log_entry
    def get_df(self) -> pd.DataFrame:
        """Return the parsed metadata DataFrame."""
        return self.datadf
