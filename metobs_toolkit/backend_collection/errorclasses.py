class MetObsMissingFile(Exception):
    """Raise when a file is missing which is required."""

    pass


class MetObsSensorDataNotFound(Exception):
    """Raise when a user request a specific sensordata
    that does not exist for a station"""

    pass


class MetObsMetadataNotFound(Exception):
    """Raise when a user request a specific metadata
    that does not exist for a station"""

    pass


class MetObsObstypeNotFound(Exception):
    """Raise when a user request an obstype that is
    unknown by the instance."""

    pass


class MetObsStationNotFound(Exception):
    """Exception raised for errors when a station is not found."""

    pass


class MetObsWrongType(Exception):
    """Exception raised when a wrong type of variable is detected."""

    pass


class MetObsDataAlreadyPresent(Exception):
    """Raise when something wants to be set, but it is already available."""

    pass


class MetobsQualityControlError(Exception):
    """Exception raised for errors in the datasetbase."""

    pass


class MetObsModelDataError(Exception):
    """Exception raised when something is wrong or missing with Modeldata"""

    pass


class MetObsTimeSimplifyError(Exception):
    """Exception raised when something is wrong with time resampling/syncing/simplifying"""

    pass


class MetObsStationClassError(Exception):
    """Raises general errors related to the station class."""

    pass


class MetObsMissingArgument(Exception):
    """Raises general errors when a argument is required for a specic situation. (like metadata-only)"""

    pass
