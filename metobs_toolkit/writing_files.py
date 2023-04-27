#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 15:30:55 2023

@author: thoverga
"""
# --
import os


def write_dataset_to_csv(
    df, metadf, filename, outputfolder, location_info, observation_types
):
    """
    Write the dataset to a file where the observations, metadata and (if available)
    the quality labels per observation type are merged together.

    A final qualty controll label for each quality-controlled-observation type
    can be added in the outputfile.

    The file will be writen to the Settings.outputfolder.

    Parameters
    ----------
    df: pandas.DataFrame
        The merged dataframe containing observations, gaps, outliers and missing timestamps.
    metadf: pandas.DataFrame
        The Dataset.metadf attribute.
    filename : string, optional
        The name of the output csv file. If none, a standard-filename is generated
        based on the period of data. The default is None.


    Returns
    -------
    None

    """

    # make column ordering 'datetime', 'name', obs, QC, Qc final metadata
    write_cols = ["datetime", "name"]
    write_cols.extend([col for col in df.columns if col in observation_types])
    write_cols.extend([col for col in df.columns if col.endswith("_label")])
    write_cols.extend(location_info)  # metadata

    df = df.reset_index()

    # merge metadata
    df = df.merge(metadf, how="left", left_on="name", right_index=True)

    # subset and order columns
    df = df[write_cols]

    # find observation type that are not present
    ignore_obstypes = [col for col in observation_types if df[col].isnull().all()]
    df = df.drop(columns=ignore_obstypes)

    # find metadata that are not present
    ignore_metadat = [col for col in location_info if df[col].isnull().all()]
    df = df.drop(columns=ignore_metadat)

    df = df.sort_values(["name", "datetime"])

    # make filename
    if isinstance(filename, type(None)):
        startstr = df["datetime"].min().strftime("%Y%m%d")
        endstr = df["datetime"].max().strftime("%Y%m%d")
        filename = "dataset_" + startstr + "_" + endstr
    else:
        if filename.endswith(".csv"):
            filename = filename[:-4]  # to avoid two times .csv.csv

    filepath = os.path.join(outputfolder, filename + ".csv")

    # write to csv in output folder
    print(f"write dataset to file: {filepath}")
    df.to_csv(path_or_buf=filepath, sep=";", na_rep="NaN", index=False)
