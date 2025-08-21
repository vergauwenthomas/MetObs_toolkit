import logging

from metobs_toolkit.backend_collection.loggingmodule import log_entry

logger = logging.getLogger("<metobs_toolkit>")


class MetObsMissingFile(Exception):
    """
    Exception raised when a required file is missing.
    """

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)


class MetObsSensorDataNotFound(Exception):
    """
    Exception raised when requested sensor data does not exist for a station.
    """

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)


class MetObsMetadataNotFound(Exception):
    """
    Exception raised when requested metadata does not exist for a station.
    """

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)


class MetObsObstypeNotFound(Exception):
    """
    Exception raised when a requested observation type is unknown.
    """

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)


class MetObsStationNotFound(Exception):
    """
    Exception raised when a station is not found.
    """

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)


class MetObsWrongType(Exception):
    """
    Exception raised when a variable is of the wrong type.
    """

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)


class MetObsDataAlreadyPresent(Exception):
    """
    Exception raised when data is already present and cannot be set again.
    """

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)


class MetObsQualityControlError(Exception):
    """
    Exception raised for errors in the dataset base.
    """

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)


class MetObsModelDataError(Exception):
    """
    Exception raised when there is an issue with model data.
    """

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)


class MetObsTimeSimplifyError(Exception):
    """
    Exception raised when there is an error with time resampling, syncing, or simplifying.
    """

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)


class MetObsStationClassError(Exception):
    """
    Exception raised for general errors related to the station class.
    """

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)


class MetObsMissingArgument(Exception):
    """
    Exception raised when an argument is required for a specific situation.
    """

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)


class MetObsGEEDatasetError(Exception):
    """
    Exception raised when there is an issue with a GEE API call.
    """

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)


class MetObsUnitsIncompatible(Exception):
    """Raised when an incompatible unit is set."""

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)


class MetObsUnitUnknown(Exception):
    """Raised when an invalid unit is set."""

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)


class MetObsTemplateError(Exception):
    """
    Exception raised for errors in the template.
    """

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)


class MetObsArgumentError(Exception):
    """Raised when an argument could not be converted to a target type."""

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)


class MetObsInconsistentStationName(Exception):
    """Special case only --> mismatch in data-metadata stationnames"""

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)


class MetObsAdditionError(Exception):
    """Raised when addition failed (often different _id())"""

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)


class MetObsNonUniqueIDs(Exception):
    """Raised when non-unique ID's are expected, but duplicates are found"""

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
