import logging
from functools import wraps

logger = logging.getLogger("<metobs_toolkit>")


def log_entry(func):
    """Decorator that logs a DEBUG message each time the decorated function is entered.

    Parameters
    ----------
    func : callable
        The function to wrap.

    Returns
    -------
    callable
        Wrapped function that logs entry before delegating to the original.
    """

    @wraps(func)
    def wrapper(*args, **kwargs):
        """Call the original function after logging the entry."""
        logger.debug(f"Entering {func.__name__}() in {func.__code__.co_filename}")
        return func(*args, **kwargs)

    return wrapper
