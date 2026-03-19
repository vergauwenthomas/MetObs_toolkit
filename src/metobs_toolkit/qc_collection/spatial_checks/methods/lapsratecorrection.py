from __future__ import annotations

import logging
from typing import List, TYPE_CHECKING

import pandas as pd


logger = logging.getLogger("<metobs_toolkit>")


if TYPE_CHECKING:
    from metobs_toolkit.qc_collection.spatial_checks.buddywrapsensor import (
        BuddyWrapSensor,
    )


def correct_lapse_rate(
    widedf: pd.DataFrame,
    wrappedsensors: List[BuddyWrapSensor],
    lapserate: float | None = None,
) -> pd.DataFrame:
    """Apply a lapse-rate altitude correction to the wide observations DataFrame.

    Each station's observations are shifted by
    ``altitude * (-1) * lapserate`` so that all values are effectively
    corrected to 0 m altitude before the buddy check is performed.
    The correction term is also stored on each wrapped sensor.

    Parameters
    ----------
    widedf : pandas.DataFrame
        Wide-format DataFrame with station names as columns and timestamps
        as index.
    wrappedsensors : list of BuddyWrapSensor
        Wrapped sensors whose ``cor_term`` and
        ``flag_lapsrate_corrections`` attributes are updated in-place.
    lapserate : float or None, optional
        Lapse rate in units per metre (e.g. ``-0.0065`` K/m for
        temperature).  If None, no correction is applied.  Default is
        None.

    Returns
    -------
    pandas.DataFrame
        Updated wide observations DataFrame with corrections applied.

    Raises
    ------
    ValueError
        If ``lapserate`` is not None and at least one station has a
        missing altitude value.
    """
    if lapserate is None:
        logger.debug("No lapse rate correction applied")
        for wrapsens in wrappedsensors:
            wrapsens.flag_lapsrate_corrections = False
    else:
        logger.debug("Applying lapse rate correction with rate: %s", lapserate)
        # Test if all stations have altitude
        has_alts = [budsta.site.flag_has_altitude() for budsta in wrappedsensors]

        if not all(has_alts):
            raise ValueError(
                "At least one station has a NaN value for 'altitude', not lapse rate correction possible"
            )
        for budsta in wrappedsensors:
            budsta.flag_lapsrate_corrections = True

            # Since buddy check works with relative differences, correct all
            # stations to the 0m altitude
            correction_term = budsta.site.altitude * (-1) * lapserate
            budsta.cor_term = correction_term  # update it in the buddy station

            # apply the correction on the wide dataframe
            widedf[budsta.name] = widedf[budsta.name] + correction_term

    return widedf
