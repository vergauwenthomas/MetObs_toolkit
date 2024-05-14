#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  8 15:55:13 2023

@author: thoverga
"""


import sys, os

from pathlib import Path
import pandas as pd

import metobs_toolkit

lib_folder = Path(__file__).resolve().parents[2]
# sys.path.append(str(lib_folder))
print(str(lib_folder))

# %%


dataset = metobs_toolkit.Dataset()


dataset = dataset.import_dataset(
    folder_path=os.path.join(str(lib_folder), "tests", "test_data"),
    filename="tests_dataset3.pkl",
)


# dataset = metobs_toolkit.Dataset()
# dataset.update_settings(input_data_file=metobs_toolkit.demo_datafile,
#                         input_metadata_file=metobs_toolkit.demo_metadatafile,
#                         template_file=metobs_toolkit.demo_template,
#                         )


# dataset.import_data_from_file()

# dataset.get_lcz()
# dataset.get_landcover()

# dataset.save_dataset(outputfolder=os.path.join(str(lib_folder), "tests", "test_data"),
#                       filename="tests_dataset3.pkl")


an = dataset.get_analysis()


# =============================================================================
# test diurnal methods
# =============================================================================

teststa = ["vlinder01", "vlinder02", "vlinder03"]

from datetime import datetime

startdt = datetime(2022, 9, 4)

# Test plotting and functions
temp_diurnal = an.get_diurnal_statistics(
    colorby="lcz", stations=teststa, startdt=startdt
)


test2 = an.get_diurnal_statistics_with_reference(
    refstation="vlinder08", colorby="name", errorbands=True
)

test3 = an.get_aggregated_cycle_statistics(aggregation=["lcz"])

# =============================================================================
# test anual cycle
# =============================================================================

# test4 = an.get_anual_statistics(groupby=['name'])


# Test values

temp_diurnal_test = {
    "Low plants (LCZ D)": {
        0: 15.539583333333333,
        1: 15.297222222222224,
        2: 15.162500000000001,
        3: 15.288888888888888,
        4: 15.211111111111112,
        5: 14.987499999999999,
        6: 15.601388888888888,
        7: 16.759027777777778,
        8: 17.994444444444444,
        9: 19.257638888888888,
        10: 20.078472222222224,
        11: 20.533333333333335,
        12: 20.994444444444444,
        13: 21.1875,
        14: 20.979166666666668,
        15: 20.907638888888886,
        16: 20.69027777777778,
        17: 20.085416666666667,
        18: 18.210416666666667,
        19: 17.056944444444444,
        20: 16.257638888888888,
        21: 15.902777777777779,
        22: 15.697222222222223,
        23: 15.400694444444444,
    },
    "Open midrise": {
        0: 15.75,
        1: 15.56076388888889,
        2: 15.36840277777778,
        3: 15.242708333333333,
        4: 15.108333333333333,
        5: 15.036458333333334,
        6: 15.344791666666667,
        7: 16.242708333333333,
        8: 17.484027777777776,
        9: 19.08090277777778,
        10: 19.8375,
        11: 20.325694444444444,
        12: 20.853819444444444,
        13: 21.334375,
        14: 21.483333333333334,
        15: 21.468402777777776,
        16: 21.009722222222223,
        17: 20.249122807017542,
        18: 18.803472222222222,
        19: 17.847569444444446,
        20: 17.038888888888888,
        21: 16.404166666666665,
        22: 15.93263888888889,
        23: 15.580208333333335,
    },
}


assert (
    temp_diurnal.eq(pd.DataFrame(temp_diurnal_test)).all().all()
), f"Maybe something wrong with the verbose output, since it is not equal to hardcoded df."
# assert stats.eq(pd.DataFrame(stats_test)).all().all(), f'Maybe something wrong with the verbose output, since it is not equal to hardcoded df.'


# =============================================================================
# Test represetation
# =============================================================================

print(an)

# =============================================================================
# Test filter method
# =============================================================================

filter_an = an.apply_filter('temp < 15.5 &  hour <= 19 & lcz == "Open midrise"')

assert filter_an.df.shape == (2481, 4), "filter on analysis problem"

# =============================================================================
# aggregate method
# =============================================================================

agg_df = an.aggregate_df(agg=["lcz", "hour"])
assert agg_df.shape == (216, 4), "aggregate on analysis problem"


# =============================================================================
# Correlation check
# =============================================================================
import numpy as np

an.get_lc_correlation_matrices(
    obstype=["temp", "humidity"], groupby_labels=["lcz", "season"]
)

# plot test
an.plot_correlation_heatmap(groupby_value=("Open lowrise", "autumn"))

# value test

cor_vals = {
    "temp": {
        "temp": 1.0,
        "humidity": -0.8530694150916298,
        "water_100m": np.nan,
        "pervious_100m": -0.16467165943725465,
        "impervious_100m": 0.16467165943725534,
    },
    "humidity": {
        "temp": -0.8530694150916298,
        "humidity": 1.0,
        "water_100m": np.nan,
        "pervious_100m": 0.20280562151890283,
        "impervious_100m": -0.2028056215189032,
    },
    "water_100m": {
        "temp": np.nan,
        "humidity": np.nan,
        "water_100m": np.nan,
        "pervious_100m": np.nan,
        "impervious_100m": np.nan,
    },
    "pervious_100m": {
        "temp": -0.16467165943725465,
        "humidity": 0.20280562151890283,
        "water_100m": np.nan,
        "pervious_100m": 1.0,
        "impervious_100m": -1.0,
    },
    "impervious_100m": {
        "temp": 0.16467165943725534,
        "humidity": -0.2028056215189032,
        "water_100m": np.nan,
        "pervious_100m": -1.0,
        "impervious_100m": 1.0,
    },
}

assert (
    an.lc_cor_dict[("Open lowrise", "autumn")]["cor matrix"]
    .fillna(0)
    .eq(pd.DataFrame(cor_vals).fillna(0))
    .all()
    .all()
), "Something wrong with the lc correlations matrices"


# scatter plot test

an.plot_correlation_variation()
