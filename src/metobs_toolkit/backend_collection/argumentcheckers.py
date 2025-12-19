"""Collection of functions that checks/changes
arguments and input passed by the user."""

import logging
import datetime as datetimemodule
import warnings

import pandas as pd

from metobs_toolkit.backend_collection.errorclasses import MetObsArgumentError
from metobs_toolkit.backend_collection.datetime_collection import to_timedelta
from metobs_toolkit.backend_collection.decorators import log_entry


logger = logging.getLogger("<metobs_toolkit>")


@log_entry
def fmt_timedelta_arg(timedeltaarg, none_is_none=True) -> pd.Timedelta:
    if (none_is_none) & (timedeltaarg is None):
        return None
    if isinstance(timedeltaarg, datetimemodule.timedelta):
        dt = to_timedelta(timedeltaarg)
    elif isinstance(timedeltaarg, str):
        dt = to_timedelta(timedeltaarg)
    elif isinstance(timedeltaarg, pd.Timedelta):
        dt = timedeltaarg
    else:
        raise MetObsArgumentError(
            f"{timedeltaarg} could not be interpreted as a Timedelta, \
convert it to a pd.Timedelta()."
        )

    if dt == pd.Timedelta(0):
        logger.warning(f"A {dt} is given as an argument for a timedelta.")

    return dt


@log_entry
def fmt_datetime_arg(
    datetimearg, tz_if_dt_is_naive="UTC", none_is_none=True
) -> pd.Timestamp:
    """Formats a datetime given by a user in an argument.

    Conversion to a pandas.Timestamp in the tz of the records.

    """

    if (none_is_none) & (datetimearg is None):
        return None
    # check type and cast to pd.Timestamp
    if isinstance(datetimearg, datetimemodule.datetime):
        dt = pd.Timestamp(datetimearg)
    elif isinstance(datetimearg, pd.Timestamp):
        dt = datetimearg
    else:
        raise MetObsArgumentError(
            f"{datetimearg} is not in a datetime format (datetime.datetime or \
pandas.Timestamp)."
        )

    # check timezone and make tz-aware
    if dt.tz is None:
        warnings.warn(
            "The datetime argument has no timezone information. "
            f"Assuming {tz_if_dt_is_naive} timezone.",
            UserWarning,
            stacklevel=2,
        )
        # tz naive timestamp
        dt = dt.tz_localize(tz=tz_if_dt_is_naive)

    return dt
