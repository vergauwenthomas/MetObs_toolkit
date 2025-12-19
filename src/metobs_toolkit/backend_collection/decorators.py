import logging
from functools import wraps

logger = logging.getLogger("<metobs_toolkit>")


def log_entry(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        logger.debug(f"Entering {func.__name__}() in {func.__code__.co_filename}")
        return func(*args, **kwargs)

    return wrapper
