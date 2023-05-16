# -*- coding: utf-8 -*-

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  25 13:44:54 2022

@author: thoverga
"""


import pandas as pd
import logging


logger = logging.getLogger(__name__)


def get_freq_statistics(comb_df, obstype, checks_info, gaps_info, applied_qc_order):
    outlier_labels = [qc["outlier_flag"] for qc in checks_info.values()]

    final_counts = comb_df['label'].value_counts()

    # add missing labels
    # QC labels
    non_triggered_labels_dict = {}
    # fill with zeros for non-triggered checks
    for outl_label in outlier_labels:
        if not outl_label in final_counts.index:
            non_triggered_labels_dict[outl_label] = 0

    # gaps
    if not gaps_info["gap"]["outlier_flag"] in final_counts.index:
        non_triggered_labels_dict[gaps_info["gap"]["outlier_flag"]] = 0

    # missing timestamps
    if not gaps_info["missing_timestamp"]["outlier_flag"] in final_counts.index:
        non_triggered_labels_dict[gaps_info["missing_timestamp"]["outlier_flag"]] = 0

    non_triggered_labels = pd.Series(non_triggered_labels_dict)
    final_counts = pd.concat([final_counts, non_triggered_labels])
    tot_n_obs = final_counts.sum()

    # to percentages
    final_counts = (final_counts / tot_n_obs) * 100.0

    # ------- aggregate outliers ----------

    # 1 agg to ok - outlier - gap - missing

    try:
        agg_ok = final_counts["ok"].squeeze()
    except KeyError:
        agg_ok = 0.0

    agg_dict = {
        "ok": agg_ok,
        "QC outliers": final_counts.loc[final_counts.index.isin(outlier_labels)].sum(),
        "missing (gaps)": final_counts[gaps_info["gap"]["outlier_flag"]].squeeze(),
        "missing (individual)": final_counts[
            gaps_info["missing_timestamp"]["outlier_flag"]
        ].squeeze(),
    }

    # 2 indevidual outliers
    outl_dict = final_counts.loc[final_counts.index.isin(outlier_labels)].to_dict()

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
            specific_outliers = final_counts.loc[checks_info[checkname]["outlier_flag"]]
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
    not_perf_checknames = [
        check for check in checks_info.keys() if not check in applied_checks
    ]
    for checkname in not_perf_checknames:
        specific_counts[checkname] = {"not checked": 100.0, "ok": 0.0, "outlier": 0.0}

    # add Gaps
    gap_specific_counts = {
        "not checked": 0,  # all obs are always checked
        "ok": 100.0 - final_counts[gaps_info["gap"]["outlier_flag"]],
        "outlier": final_counts[gaps_info["gap"]["outlier_flag"]],
    }
    specific_counts[gaps_info["gap"]["label_columnname"]] = gap_specific_counts

    # misssing timestamps
    missing_specific_counts = {
        "not checked": 0,  # all obs are always checked
        "ok": 100.0 - final_counts[gaps_info["missing_timestamp"]["outlier_flag"]],
        "outlier": final_counts[gaps_info["missing_timestamp"]["outlier_flag"]],
    }
    specific_counts[
        gaps_info["missing_timestamp"]["label_columnname"]
    ] = missing_specific_counts

    return (agg_dict, outl_dict, specific_counts)
