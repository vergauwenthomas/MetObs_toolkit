import logging

import pandas as pd
from typing import List, Union
import metobs_toolkit.backend_collection.printing_collection as printing

logger = logging.getLogger("<metobs_toolkit>")


class SensorWhiteSet:
    """Whitelist container for a single sensor (station-obstype combination).

    This class manages whitelisted timestamps for a specific sensor, allowing
    certain observations to be excluded from outlier detection in QC checks.

    Parameters
    ----------
    white_timestamps : list, optional
        List of datetime objects to whitelist. Default is empty list.
    all_timestamps : bool, optional
        If True, all timestamps are whitelisted. Default is False.

    Attributes
    ----------
    white_timestamps : list
        List of whitelisted datetime objects.
    all_timestamps : bool
        Whether all timestamps are whitelisted.
    """

    def __init__(
        self, white_timestamps: List = [], all_timestamps: bool = False
    ) -> None:

        if all_timestamps:
            assert (
                len(white_timestamps) == 0
            ), "If all_timestamps is True, white_timestamps must be empty."

        self.white_timestamps = white_timestamps
        self.all_timestamps = all_timestamps

        logger.debug(
            "Initialized SensorWhiteSet: all_timestamps=%s, n_timestamps=%s",
            all_timestamps,
            len(white_timestamps),
        )

    def __repr__(self) -> str:
        """Return a string representation for debugging."""
        if self.all_timestamps:
            return f"{type(self).__name__}(all_timestamps=True)"
        return f"{type(self).__name__}(n_timestamps={len(self.white_timestamps)})"

    def __str__(self) -> str:
        return self.__repr__()

    def has_whites(self) -> bool:
        """Check if any timestamps are whitelisted.

        Returns
        -------
        bool
            True if any timestamps are whitelisted, False otherwise.
        """
        if self.all_timestamps:
            return True
        if len(self.white_timestamps) > 0:
            return True
        return False

    def all_timestamps_white(self) -> bool:
        """Check if all timestamps are whitelisted.

        Returns
        -------
        bool
            True if all timestamps are whitelisted, False otherwise.
        """
        return self.all_timestamps

    def get_white_timestamps(self) -> pd.DatetimeIndex:
        """Get whitelisted timestamps as a DatetimeIndex.

        Returns
        -------
        pd.DatetimeIndex
            DatetimeIndex containing all whitelisted timestamps.
        """
        return pd.DatetimeIndex(data=self.white_timestamps, name="datetime")

    def catch_white_records(self, outliers_idx: pd.DatetimeIndex) -> pd.DatetimeIndex:
        """Remove whitelisted timestamps from outliers index.

        Filters the provided outliers index by removing any timestamps that are
        whitelisted in this SensorWhiteSet.

        Parameters
        ----------
        outliers_idx : pd.DatetimeIndex
            Index of outlier timestamps to filter.

        Returns
        -------
        pd.DatetimeIndex
            Filtered outliers index with whitelisted records removed.
        """
        logger.debug(
            "Filtering %s outliers with SensorWhiteSet (has_whites=%s)",
            len(outliers_idx),
            self.has_whites(),
        )

        if self.has_whites():

            if self.all_timestamps_white():
                # all timestamps are white, return empty index
                logger.debug(
                    "All timestamps whitelisted, removing all %s outliers",
                    len(outliers_idx),
                )
                outliers = pd.DatetimeIndex([], name="datetime")
            else:
                # Get the white timestamps
                white_records = self.get_white_timestamps()
                # Remove white records from outliers
                outliers = outliers_idx.difference(white_records)
                n_filtered = len(outliers_idx) - len(outliers)
                logger.debug(
                    "Filtered out %s whitelisted outliers, %s outliers remain",
                    n_filtered,
                    len(outliers),
                )
        else:
            # no whites
            outliers = outliers_idx
            logger.debug(
                "No whitelist applied, all %s outliers retained", len(outliers)
            )

        return outliers


class WhiteSet:
    """Whitelist container for multiple stations and observation types.

    This class manages a collection of whitelisted records across multiple stations
    and observation types. It uses a pandas Index or MultiIndex with optional levels
    for 'name' (station), 'obstype', and 'datetime' to define which records should
    be excluded from outlier detection in QC checks.

    Parameters
    ----------
    white_records : pd.Index, optional
        Index with levels 'name', 'obstype', and/or 'datetime' defining whitelisted
        records. Default is an empty Index.

    Attributes
    ----------
    white_records : pd.Index
        The whitelist index containing whitelisted record identifiers.

    Notes
    -----
    The white_records index must contain at least one of: 'name', 'obstype', or
    'datetime' as level names. If 'datetime' is not present, all timestamps for
    matching station/obstype combinations are whitelisted.
    """

    def __init__(self, white_records: pd.Index = pd.Index([])) -> None:
        self.white_records = white_records

        # Validate white_records structure
        self.test_white_records()

        if not self.is_empty():
            logger.debug(
                "Initialized WhiteSet: n_records=%s, levels=%s",
                len(white_records),
                list(white_records.names),
            )
        else:
            logger.debug("Initialized empty WhiteSet")

    def __repr__(self) -> str:
        if self.is_empty():
            return f"{type(self).__name__}(empty)"
        levels = list(self.white_records.names)
        n_records = len(self.white_records)
        return f"{type(self).__name__}(n_records={n_records}, levels={levels})"

    def __str__(self) -> str:
        return self.__repr__()

    def test_white_records(self) -> None:
        """Validate the structure and content of white_records index.

        Validates that white_records has a valid structure for use in QC methods.
        The index must contain at least one of the expected level names ('name',
        'obstype', or 'datetime') and may not contain any unexpected levels.

        Raises
        ------
        ValueError
            If white_records does not contain at least one of 'name', 'obstype', or
            'datetime' as index level names.
        ValueError
            If white_records contains unexpected index levels.
        """
        if self.is_empty():
            return

        logger.debug(
            "Validating WhiteSet structure with levels: %s", self.white_records.names
        )

        if not any(
            [
                idxname in self.white_records.names
                for idxname in ["name", "obstype", "datetime"]
            ]
        ):
            logger.debug("Validation failed: missing required index levels")
            raise ValueError(
                "white_records must contain at least one of the following index levels: 'name', 'obstype', 'datetime'"
            )
        if not all(
            [
                idxname in ["name", "obstype", "datetime"]
                for idxname in self.white_records.names
            ]
        ):
            logger.debug("Validation failed: unexpected index levels found")
            raise ValueError(
                "white_records contains unexpected index levels. Only 'name', 'obstype', and 'datetime' are allowed."
            )

        logger.debug("WhiteSet validation passed")

    def is_empty(self) -> bool:
        """Check if white_records is empty.

        Returns
        -------
        bool
            True if white_records is empty, False otherwise.
        """
        return self.white_records.empty

    def get_info(self, printout: bool = True) -> Union[str, None]:
        """
        Retrieve and optionally print detailed information about the WhiteSet.

        Parameters
        ----------
        printout : bool, optional
            If True, prints the information to the console. If False, returns
            the information as a string. Default is True.

        Returns
        -------
        str or None
            A string containing the WhiteSet information if `printout` is False.
            Otherwise, returns None.
        """
        infostr = ""
        infostr += printing.print_fmt_title("General info of WhiteSet")

        if self.is_empty():
            infostr += printing.print_fmt_section("Whitelist details")
            infostr += printing.print_fmt_line(
                "Empty WhiteSet (no whitelisted records)"
            )
        else:
            # Basic information
            infostr += printing.print_fmt_section("Whitelist details")
            n_records = len(self.white_records)
            infostr += printing.print_fmt_line(f"Total records: {n_records}")

            # Index levels information
            levels = [lvl for lvl in self.white_records.names if lvl is not None]
            infostr += printing.print_fmt_line(f"Index levels: {', '.join(levels)}")

            # Count unique values per level
            if "name" in levels:
                n_stations = self.white_records.get_level_values("name").nunique()
                stations = sorted(self.white_records.get_level_values("name").unique())
                infostr += printing.print_fmt_line(
                    f"Stations ({n_stations}): {', '.join(map(str, stations))}"
                )

            if "obstype" in levels:
                n_obstypes = self.white_records.get_level_values("obstype").nunique()
                obstypes = sorted(
                    self.white_records.get_level_values("obstype").unique()
                )
                infostr += printing.print_fmt_line(
                    f"Observation types ({n_obstypes}): {', '.join(map(str, obstypes))}"
                )

            if "datetime" in levels:
                n_times = self.white_records.get_level_values("datetime").nunique()
                infostr += printing.print_fmt_line(f"Unique timestamps: {n_times}")

                # Time range information
                datetimes = self.white_records.get_level_values("datetime")
                min_time = datetimes.min()
                max_time = datetimes.max()
                infostr += printing.print_fmt_line(
                    f"Time range: {min_time} to {max_time}", 1
                )
            else:
                infostr += printing.print_fmt_line(
                    "All timestamps whitelisted (no 'datetime' level)"
                )

        if printout:
            print(infostr)
        else:
            return infostr

    def create_sensorwhitelist(
        self, trg_station: str, trg_obstype: str
    ) -> SensorWhiteSet:
        """Create a sensor-specific whitelist for a station and observation type.

        Filters the white_records by station name and obstype to create a
        SensorWhiteSet containing only the relevant whitelisted timestamps.

        Parameters
        ----------
        trg_station : str
            Target station name to filter for.
        trg_obstype : str
            Target observation type to filter for.

        Returns
        -------
        SensorWhiteSet
            A SensorWhiteSet instance containing whitelisted timestamps for the
            specified station and obstype combination.

        Notes
        -----
        If the white_records index does not contain a 'datetime' level but does
        match the station/obstype, all timestamps are whitelisted for that sensor.
        """
        logger.debug(
            "Creating SensorWhiteSet for station='%s', obstype='%s'",
            trg_station,
            trg_obstype,
        )

        if self.is_empty():
            logger.debug("WhiteSet is empty, returning empty SensorWhiteSet")
            return SensorWhiteSet(white_timestamps=[], all_timestamps=False)

        # Filter white_records for the target station and obstype
        trg_whitelist = self.white_records

        if "name" in trg_whitelist.names:
            if trg_station in trg_whitelist.get_level_values("name"):
                trg_whitelist = trg_whitelist[
                    (trg_whitelist.get_level_values("name") == trg_station)
                ]
                logger.debug(
                    "Filtered whitelist by station name, %s records remain",
                    len(trg_whitelist),
                )
            else:
                # name is specified, but no matches in whitelist
                logger.debug(
                    "Station '%s' not found in whitelist, returning empty SensorWhiteSet",
                    trg_station,
                )
                return SensorWhiteSet(white_timestamps=[], all_timestamps=False)

        # filter on obstype if present
        if "obstype" in trg_whitelist.names:
            if trg_obstype in trg_whitelist.get_level_values("obstype"):
                trg_whitelist = trg_whitelist[
                    (trg_whitelist.get_level_values("obstype") == trg_obstype)
                ]
                logger.debug(
                    "Filtered whitelist by obstype, %s records remain",
                    len(trg_whitelist),
                )
            else:
                # obstype is specified, but no matches in whitelist
                logger.debug(
                    "Obstype '%s' not found in whitelist, returning empty SensorWhiteSet",
                    trg_obstype,
                )
                return SensorWhiteSet(white_timestamps=[], all_timestamps=False)

        if "datetime" in trg_whitelist.names:
            # Get the white datetimes
            white_datetimes = pd.DatetimeIndex(
                trg_whitelist.get_level_values("datetime").unique()
            )
            logger.debug(
                "Created SensorWhiteSet with %s unique timestamps", len(white_datetimes)
            )
            return SensorWhiteSet(
                white_timestamps=white_datetimes, all_timestamps=False
            )
        else:
            # if no datetime level is set, and name and/or obstype match, all timestamps are white
            logger.debug(
                "No datetime level in whitelist, all timestamps whitelisted for this sensor"
            )
            return SensorWhiteSet(white_timestamps=[], all_timestamps=True)
