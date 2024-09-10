# -*- coding: utf-8 -*-

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module for computing frequency statistics of outlier labels.

@author: thoverga
"""


import pandas as pd
import logging

from metobs_toolkit.settings_files.default_formats_settings import (
    label_def,
    qc_label_group,
    gapfill_label_group,
    failed_gapfill_label_group,
)

logger = logging.getLogger(__name__)


def get_freq_statistics(comb_df, obstype, applied_qc_order):
    """Compute frequency statistics of the outliers.

    Parameters
    ----------
    comb_df : pandas.DataFrame
        The dataframe containing all observations, outliers and their labels.
    obstype : str
        The observation type to compute the frequencies of.
    applied_qc_order : pandas.DataFrame
        The _applied_qc attribute of the Dataset.

    Returns
    -------
    agg_dict : dict
        Dictionary containing occurrence frequencies for all labels.
    outl_dict : dict
        Dictionary with frequency statistics of outlier-labels.
    specific_counts : dict
        Dictionary containing the effectiveness of quality control checks
        individually.

    """
    qc_outlier_labels = [label_def[qc_check]["label"] for qc_check in qc_label_group]
    final_counts = comb_df["label"].value_counts()

    # Get all possible labels for a record:
    all_possible_labels = {
        val["label"] for val in label_def.values() if "agg_def" not in val.keys()
    }

    # all_missing = all_possible_labels-set(final_counts.index)
    all_missing_series = pd.Series(
        index=all_possible_labels, data=[0] * len(all_possible_labels)
    )

    final_counts = pd.concat([final_counts, all_missing_series])
    final_counts = final_counts[~final_counts.index.duplicated(keep="first")]

    tot_n_obs = final_counts.sum()

    # to percentages
    final_counts = (final_counts / tot_n_obs) * 100.0

    # ------- aggregate outliers ----------

    # 1 agg to ok - outlier - gap (filled/unfilled)

    try:
        agg_ok = final_counts[label_def["goodrecord"]["label"]].squeeze()
    except KeyError:
        agg_ok = 0.0

    gap_labels = [label_def["regular_gap"]["label"]]
    gap_labels.extend(
        [label_def[fill_method]["label"] for fill_method in gapfill_label_group]
    )  # all filled labels
    gap_labels.extend(
        [
            label_def[failed_fill_method]["label"]
            for failed_fill_method in failed_gapfill_label_group
        ]
    )  # all failed filled labels

    qc_specific_outliers = [
        label_def[qc_method]["label"] for qc_method in qc_label_group
    ]

    agg_dict = {
        "ok": agg_ok,
        "QC outliers": final_counts.loc[
            final_counts.index.isin(qc_specific_outliers)
        ].sum(),
        "gaps (filled/unfilled)": final_counts[
            final_counts.index.isin(gap_labels)
        ].sum(),
    }

    # 2 indevidual outliers
    outl_dict = final_counts.loc[final_counts.index.isin(qc_outlier_labels)].to_dict()

    # 3 Effectivenes per check

    specific_counts = {}
    # Note: some complexity because observations can be removed by privious executed checsk,
    # so construct the counts in the order of the applied checks

    applied_qc_order = (
        applied_qc_order.drop_duplicates()
    )  # when qc applied mulitple times on same obstype
    applied_checks = applied_qc_order.loc[applied_qc_order["obstype"] == obstype][
        "checkname"
    ].to_list()

    percent_rejected_before = 0.0

    for checkname in applied_checks:
        try:
            specific_outliers = final_counts.loc[label_def[checkname]["label"]]
        except KeyError:
            specific_outliers = 0.0

        not_checked = percent_rejected_before
        ok = 100.0 - specific_outliers - not_checked

        specific_counts[checkname] = {
            "not checked": not_checked,
            "ok": ok,
            "outlier": specific_outliers,
        }

        percent_rejected_before += specific_outliers

    # add checks that are not performed
    for qc_checkname in qc_label_group:
        if qc_checkname not in specific_counts.keys():
            specific_counts[qc_checkname] = {
                "not checked": 100.0,
                "ok": 0.0,
                "outlier": 0.0,
            }

    # add Gaps
    gap_specific_counts = {
        "not checked": 0,  # all obs are always checked
        "ok": 100.0 - final_counts[label_def["regular_gap"]["label"]],
        "outlier": final_counts[label_def["regular_gap"]["label"]],
    }
    specific_counts["Gap finder"] = gap_specific_counts

    return (agg_dict, outl_dict, specific_counts)


# =============================================================================
# Docstring test
# =============================================================================
if __name__ == "__main__":
    from metobs_toolkit.doctest_fmt import setup_and_run_doctest

    setup_and_run_doctest()
