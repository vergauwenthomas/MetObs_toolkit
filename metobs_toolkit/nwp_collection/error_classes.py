class MetObsFieldNotFound(Exception):
    """
    Exception raised when requested sensor data does not exist for a station.
    """

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)