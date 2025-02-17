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
