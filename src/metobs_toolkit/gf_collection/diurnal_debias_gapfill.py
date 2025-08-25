import logging
import pandas as pd
import numpy as np

from metobs_toolkit.backend_collection.loggingmodule import log_entry

logger = logging.getLogger("<metobs_toolkit>")


@log_entry
def compute_diurnal_biases(df: pd.DataFrame) -> pd.DataFrame:
    """
    Compute diurnal biases for the given DataFrame.

    Calculates the mean and count of the difference between 'modelvalue' and 'value'
    for each unique combination of hour, minute, and second, considering only rows
    labeled as 'lead' or 'trail'.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame containing at least the columns 'hour', 'min', 'sec', 'diff', and 'label'.

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns 'hour', 'min', 'sec', 'diurnalbias', and 'samplesize'.
    """
    # calculate biases for diurnal records
    trainset = df.loc[df["label"].isin(["lead", "trail"])]
    diurnalbias = (
        trainset[["hour", "min", "sec", "diff"]]
        .groupby(["hour", "min", "sec"])["diff"]
        .agg(["mean", "count"])
        .reset_index()
    )

    # rename for clarity
    diurnalbias = diurnalbias.rename(
        columns={"mean": "diurnalbias", "count": "samplesize"}
    )

    return diurnalbias


@log_entry
def fill_with_diurnal_debias(df: pd.DataFrame, min_sample_size: int) -> pd.DataFrame:
    """
    Fill missing values in a DataFrame using diurnal debiasing.

    Adds time columns, computes instantaneous biases, merges diurnal biases,
    applies corrections, and fills values based on sample size.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame with a DateTimeIndex and columns 'modelvalue', 'value', and 'label'.
    min_sample_size : int
        Minimum required sample size for applying the correction.

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns: 'value', 'label', 'modelvalue', 'fillvalue', and 'msg'.
    """
    # add hour, minute, second columns
    df["hour"] = df.index.hour
    df["min"] = df.index.minute
    df["sec"] = df.index.second

    # calculate instantaneous biases
    df["diff"] = df["modelvalue"] - df["value"]

    diurnalbias = compute_diurnal_biases(df)

    # merge biases to the df
    df = df.reset_index().merge(
        diurnalbias[["diurnalbias", "samplesize", "hour", "min", "sec"]],
        how="left",
        on=["hour", "min", "sec"],
    )

    df["samplesize"] = df["samplesize"].fillna(0)
    df["correction"] = -1.0 * df["diurnalbias"]

    df["fillvalue"] = df["modelvalue"] + df["correction"]
    df["msg"] = df.apply(
        lambda x: f"diurnal bias corrected: {x['modelvalue']:.2f} + {x['correction']:.2f} based on sample of {x['samplesize']}.",
        axis=1,
    )

    # overwrite for too small sample sizes
    df.loc[df["samplesize"] < min_sample_size, "fillvalue"] = np.nan
    df.loc[df["samplesize"] < min_sample_size, "msg"] = df.loc[
        df["samplesize"] < min_sample_size, "samplesize"
    ].apply(
        lambda x: f"diurnal sample size of {x} is smaller than the minimum targeted samplesize ({min_sample_size})"
    )

    # set datetime back as index
    df = df.set_index("datetime")
    return df[["value", "label", "modelvalue", "fillvalue", "msg"]]


@log_entry
def fill_with_weighted_diurnal_debias(
    df: pd.DataFrame, min_lead_sample_size: int, min_trail_sample_size: int
) -> pd.DataFrame:
    """
    Fill missing values using weighted diurnal debiasing.

    Computes diurnal biases separately for 'lead' and 'trail' periods, applies weighted corrections,
    and fills values based on minimum sample sizes for each period.

    Parameters
    ----------
    df : pandas.DataFrame
        Input DataFrame with a DateTimeIndex and columns 'modelvalue', 'value', and 'label'.
    min_lead_sample_size : int
        Minimum required sample size for the 'lead' period.
    min_trail_sample_size : int
        Minimum required sample size for the 'trail' period.

    Returns
    -------
    pandas.DataFrame
        DataFrame with columns: 'value', 'label', 'modelvalue', 'fillvalue', and 'msg'.
    """
    # add hour, minute, second columns
    df["hour"] = df.index.hour
    df["min"] = df.index.minute
    df["sec"] = df.index.second

    # calculate instantaneous biases
    df["diff"] = df["modelvalue"] - df["value"]

    # Compute biases separately for the leading and trailing period
    for trainlabel in ["lead", "trail"]:
        trainset = df.loc[df["label"] == trainlabel]
        diurnalbias = (
            trainset[["hour", "min", "sec", "diff"]]
            .groupby(["hour", "min", "sec"])["diff"]
            .agg(["mean", "count"])
            .reset_index()
        )

        # rename for clarity
        diurnalbias = diurnalbias.rename(
            columns={
                "mean": f"diurnalbias_{trainlabel}",
                "count": f"samplesize_{trainlabel}",
            }
        )

        # merge biases to the df
        df = df.reset_index().merge(
            diurnalbias[
                [
                    f"diurnalbias_{trainlabel}",
                    f"samplesize_{trainlabel}",
                    "hour",
                    "min",
                    "sec",
                ]
            ],
            how="left",
            on=["hour", "min", "sec"],
        )

        df[f"samplesize_{trainlabel}"] = df[f"samplesize_{trainlabel}"].fillna(0)

    # compute weights for each gap record
    gapstart = df.loc[df["label"] == "gap"]["datetime"].min()
    gapend = df.loc[df["label"] == "gap"]["datetime"].max()
    gapduration = gapend - gapstart
    df["lead_weight"] = 1.0 - ((df["datetime"] - gapstart) / gapduration)
    df["trail_weight"] = 1.0 - df["lead_weight"]

    # Compute the correction
    df["correction"] = -1.0 * (
        (df["lead_weight"] * df["diurnalbias_lead"])
        + (df["trail_weight"] * df["diurnalbias_trail"])
    )

    # Compute fill values
    df["fillvalue"] = df["modelvalue"] + df["correction"]
    # overwrite fill value and msg for too small sample sizes
    df.loc[df["samplesize_lead"] < min_lead_sample_size, "fillvalue"] = np.nan
    df.loc[df["samplesize_trail"] < min_trail_sample_size, "fillvalue"] = np.nan

    # Write message
    @log_entry
    def msg_writer(row):
        """
        Write a message describing the fill status for each row.

        Parameters
        ----------
        row : pandas.Series
            Row of the DataFrame.

        Returns
        -------
        str
            Message describing the fill status.
        """
        if (row["samplesize_lead"] >= min_lead_sample_size) & (
            row["samplesize_trail"] >= min_trail_sample_size
        ):
            # Fill conditions are both met
            msg = f"Weighted diurnal bias corrected:\
{row['modelvalue']:.2f} + (-1. * ({row['lead_weight']:.2f} * {row['diurnalbias_lead']:.2f}) \
+ ({row['trail_weight']:.2f} * {row['diurnalbias_trail']:.2f})) \
based on sample of {row['samplesize_lead']} (lead) and {row['samplesize_trail']} (trail)."

        elif (row["samplesize_lead"] < min_lead_sample_size) & (
            row["samplesize_trail"] >= min_trail_sample_size
        ):
            # lead condition not met, trail met
            msg = f"Lead diurnal sample size is too small {row['samplesize_lead']} < {min_lead_sample_size}."

        elif (row["samplesize_lead"] >= min_lead_sample_size) & (
            row["samplesize_trail"] < min_trail_sample_size
        ):
            # lead condition met, trail not met
            msg = f"Trail diurnal sample size is too small {row['samplesize_trail']} < {min_trail_sample_size}."

        elif (row["samplesize_lead"] < min_lead_sample_size) & (
            row["samplesize_trail"] < min_trail_sample_size
        ):
            # both lead and trail conditions not met
            msg = f"Both lead and trail diurnal sample sizes are too small {row['samplesize_lead']} < {min_lead_sample_size} and {row['samplesize_trail']} < {min_trail_sample_size}."

        return msg

    df["msg"] = df.apply(lambda x: msg_writer(x), axis=1)  # Write messages with details

    # set datetime back as index
    df = df.set_index("datetime")
    return df[["value", "label", "modelvalue", "fillvalue", "msg"]]
