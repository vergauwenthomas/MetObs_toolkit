import pandas as pd
import numpy as np


def compute_diurnal_biases(df):

    # calculate biases for dirunal records
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


def fill_with_diurnal_debias(df, min_sample_size):

    # add hour, minute, second columns
    df["hour"] = df.index.hour
    df["min"] = df.index.minute
    df["sec"] = df.index.second

    # calculate instantanious biases
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

    # overwrite for too small samplsizes
    df.loc[df["samplesize"] < min_sample_size, "fillvalue"] = np.nan
    df.loc[df["samplesize"] < min_sample_size, "msg"] = df.loc[
        df["samplesize"] < min_sample_size, "samplesize"
    ].apply(
        lambda x: f"diurnal sample size of {x} is smaller than the minimum targeted samplesize ({min_sample_size})"
    )

    # set datetime back as index
    df = df.set_index("datetime")
    return df[["value", "label", "modelvalue", "fillvalue", "msg"]]


def fill_with_weighted_diurnal_debias(df, min_lead_sample_size, min_trail_sample_size):
    # add hour, minute, second columns
    df["hour"] = df.index.hour
    df["min"] = df.index.minute
    df["sec"] = df.index.second

    # calculate instantanious biases
    df["diff"] = df["modelvalue"] - df["value"]

    # Compute biases seperatly for the leading and trailing period
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
    # overwrite fill value and msg for too small samplsizes
    df.loc[df["samplesize_lead"] < min_lead_sample_size, "fillvalue"] = np.nan
    df.loc[df["samplesize_trail"] < min_trail_sample_size, "fillvalue"] = np.nan

    # Write message
    def msg_writer(row):
        if (
            row["samplesize_lead"]
            >= min_lead_sample_size & row["samplesize_trail"]
            >= min_trail_sample_size
        ):
            # Fill conditions are both met
            msg = f"Weighted diurnal bias corrected:\
{row['modelvalue']:.2f} + (-1. * ({row['lead_weight']:.2f} * {row['diurnalbias_lead']:.2f}) \
+ ({row['trail_weight']:.2f} * {row['diurnalbias_trail']:.2f})) \
based on sample of {row['samplesize_lead']} (lead) and {row['samplesize_trail']} (trail)."

        elif (
            row["samplesize_lead"]
            < min_lead_sample_size & row["samplesize_trail"]
            >= min_trail_sample_size
        ):
            # lead condition not met, trail met
            msg = f"Lead diurnal sample size is to small {row['samplesize_lead']} < {min_lead_sample_size}."

        elif (
            row["samplesize_lead"]
            >= min_lead_sample_size & row["samplesize_trail"]
            < min_trail_sample_size
        ):
            # lead condition met, trail not met
            msg = f"Trail diurnal sample size is to small {row['samplesize_trail']} < {min_trail_sample_size}."

        elif (
            row["samplesize_lead"]
            < min_lead_sample_size & row["samplesize_trail"]
            < min_trail_sample_size
        ):
            # both lead and trail conditions not met
            msg = f"Both lead and trail diurnal sample sizes are too small {row['samplesize_lead']} < {min_lead_sample_size} and {row['samplesize_trail']} < {min_trail_sample_size}."

        return msg

    df["msg"] = df.apply(lambda x: msg_writer(x), axis=1)  # Write messages with details

    # set datetime back as index
    df = df.set_index("datetime")
    return df[["value", "label", "modelvalue", "fillvalue", "msg"]]
