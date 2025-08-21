import logging
import pandas as pd

from metobs_toolkit.backend_collection.loggingmodule import log_entry

logger = logging.getLogger("<metobs_toolkit>")


@log_entry
def fill_regular_debias(df: pd.DataFrame) -> pd.DataFrame:
    """
    Fill missing values in a DataFrame by applying a regular debiasing correction.

    The function calculates the mean bias between the 'modelvalue' and 'value' columns
    for rows labeled as 'lead' or 'trail'. This bias is then used to correct the 'modelvalue'
    for all rows, and the corrected values are stored in a new column 'fillvalue'.
    Additional columns 'correction' and 'msg' are also added to the DataFrame.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame containing at least the columns 'value', 'label', and 'modelvalue'.

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns: 'value', 'label', 'modelvalue', 'fillvalue', and 'msg'.
    """
    trainset = df.loc[df["label"].isin(["lead", "trail"])]
    biasvalue = (trainset["modelvalue"] - trainset["value"]).mean()

    df["correction"] = -1.0 * biasvalue
    df["fillvalue"] = df["modelvalue"] + df["correction"]
    df["msg"] = df.apply(
        lambda x: f"bias corrected: {x['modelvalue']:.2f} + {x['correction']:.2f}",
        axis=1,
    )
    return df[["value", "label", "modelvalue", "fillvalue", "msg"]]
