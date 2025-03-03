import logging
import pandas as pd


def fill_regular_debias(df):
    trainset = df.loc[df["label"].isin(["lead", "trail"])]
    biasvalue = (trainset["modelvalue"] - trainset["value"]).mean()

    df["correction"] = -1.0 * biasvalue
    df["fillvalue"] = df["modelvalue"] + df["correction"]
    df["msg"] = df.apply(
        lambda x: f"bias corrected: {x['modelvalue']:.2f} + {x['correction']:.2f}",
        axis=1,
    )
    return df[["value", "label", "modelvalue", "fillvalue", "msg"]]
