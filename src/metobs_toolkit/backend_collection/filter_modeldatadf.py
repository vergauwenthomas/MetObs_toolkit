import logging
import pandas as pd
from metobs_toolkit.backend_collection.errorclasses import (
    MetObsObstypeNotFound,
    MetObsModelDataError,
)
from metobs_toolkit.backend_collection.decorators import log_entry

logger = logging.getLogger("<metobs_toolkit>")


@log_entry
def filter_modeldatadf(
    modeldatadf: pd.DataFrame, obstype: str, modelname: str, modelvariable: str
) -> pd.DataFrame:
    """
    Filter model data DataFrame by observation type, model name, and variable.

    Parameters
    ----------
    modeldatadf : pd.DataFrame
        Multi-index DataFrame containing model data with 'obstype' level.
    obstype : str
        Target observation type to filter for.
    modelname : str or None
        Model name to filter for. If None, no filtering by model name.
    modelvariable : str or None
        Model variable to filter for. If None, no filtering by variable.

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame containing only data matching the specified criteria.

    Notes
    -----
    If multiple model names or variables remain after filtering, the function
    logs a warning and uses the first occurrence.
    """

    # filter on obstype
    if obstype not in modeldatadf.index.get_level_values("obstype"):
        raise MetObsObstypeNotFound(f"There is no modeldata present of {obstype}")
    modeldatadf = modeldatadf.xs(obstype, level="obstype", drop_level=False)

    # filter on modelname
    if modelname is not None:
        if modelname not in modeldatadf["modelname"].values:
            raise MetObsModelDataError(f"There is no modeldata present of {modelname}")
        else:
            modeldatadf = modeldatadf[modeldatadf["modelname"] == modelname]

    # filter on modelvariable
    if modelvariable is not None:
        if modelvariable not in modeldatadf["modelvariable"].values:
            raise MetObsModelDataError(
                f"There is no modeldata present of {modelvariable}"
            )
        else:
            modeldatadf = modeldatadf[modeldatadf["modelvariable"] == modelvariable]

    # If there are multiple model names or variables, warn and take first occurrence
    if len(modeldatadf["modelname"].unique()) > 1:
        unique_models = modeldatadf["modelname"].unique()
        logger.warning(
            f"Multiple model names found: {unique_models}. Using first occurrence: {unique_models[0]}"
        )
        modeldatadf = modeldatadf[modeldatadf["modelname"] == unique_models[0]]

    if len(modeldatadf["modelvariable"].unique()) > 1:
        unique_vars = modeldatadf["modelvariable"].unique()
        logger.warning(
            f"Multiple model variables found: {unique_vars}. Using first occurrence: {unique_vars[0]}"
        )
        modeldatadf = modeldatadf[modeldatadf["modelvariable"] == unique_vars[0]]

    return modeldatadf
