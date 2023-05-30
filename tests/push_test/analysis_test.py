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


testdatafile = os.path.join(
    str(lib_folder), "tests", "test_data", "testdata_okt_small.csv"
)
static_data = os.path.join(str(lib_folder), "static_data", "vlinder_metadata.csv")

lcz_dict = {
    "vlinder01": "Low plants (LCZ D)",
    "vlinder02": "Open midrise",
    "vlinder03": "Open midrise",
    "vlinder04": "Sparsely built",
    "vlinder05": "Water (LCZ G)",
    "vlinder06": "Scattered Trees (LCZ B)",
    "vlinder07": "Compact midrise",
    "vlinder08": "Compact midrise",
    "vlinder09": "Scattered Trees (LCZ B)",
    "vlinder10": "Compact midrise",
    "vlinder11": "Open lowrise",
    "vlinder12": "Open highrise",
    "vlinder13": "Compact midrise",
    "vlinder14": "Low plants (LCZ D)",
    "vlinder15": "Sparsely built",
    "vlinder16": "Water (LCZ G)",
    "vlinder17": "Scattered Trees (LCZ B)",
    "vlinder18": "Low plants (LCZ D)",
    "vlinder19": "Compact midrise",
    "vlinder20": "Compact midrise",
    "vlinder21": "Sparsely built",
    "vlinder22": "Low plants (LCZ D)",
    "vlinder23": "Low plants (LCZ D)",
    "vlinder24": "Dense Trees (LCZ A)",
    "vlinder25": "Water (LCZ G)",
    "vlinder26": "Open midrise",
    "vlinder27": "Compact midrise",
    "vlinder28": "Open lowrise",
}


dataset = metobs_toolkit.Dataset()
dataset.update_settings(
    input_data_file=testdatafile,
    input_metadata_file=static_data,
    # data_template_file= template,
    output_folder="/home/thoverga/Documents",
)


dataset.import_data_from_file()

dataset.metadf["lcz"] = pd.Series(lcz_dict)  # to avoid gee interaction
# %%


an = dataset.get_analysis()

# %%
teststa = ["vlinder01", "vlinder02", "vlinder03"]

from datetime import datetime

startdt = datetime(2022, 10, 6)

# %% Test plotting and functions
(temp_diurnal, stats) = an.get_diurnal_statistics(
    colorby="lcz", stations=teststa, startdt=startdt, verbose=True
)


test2 = an.get_diurnal_statistics_with_reference(
    refstation="vlinder08", colorby="name", errorbands=True, verbose=True
)

test3 = an.get_aggregated_diurnal_statistics(aggregation=["lcz"], verbose=True)


# %% Test values

temp_diurnal_test = {
    "vlinder01": {
        0: 8.927777777777777,
        1: 8.706944444444444,
        2: 8.397222222222222,
        3: 8.05,
        4: 7.697222222222223,
        5: 7.591666666666667,
        6: 7.668055555555556,
        7: 9.37361111111111,
        8: 11.205555555555556,
        9: 12.833333333333334,
        10: 15.100000000000001,
        11: 15.808333333333334,
        12: 16.266666666666666,
        13: 16.59027777777778,
        14: 16.12222222222222,
        15: 15.865277777777777,
        16: 15.155714285714287,
        17: 12.578333333333335,
        18: 10.593333333333334,
        19: 9.398333333333333,
        20: 8.841666666666667,
        21: 8.708571428571428,
        22: 8.504166666666666,
        23: 7.872222222222222,
    },
    "vlinder02": {
        0: 9.422222222222222,
        1: 8.9875,
        2: 8.431944444444445,
        3: 8.013043478260869,
        4: 8.129577464788733,
        5: 8.005555555555555,
        6: 7.86388888888889,
        7: 9.452777777777778,
        8: 11.801388888888889,
        9: 14.016666666666667,
        10: 15.751388888888888,
        11: 16.479166666666668,
        12: 17.047222222222224,
        13: 17.109722222222224,
        14: 16.47222222222222,
        15: 16.140277777777776,
        16: 15.333333333333334,
        17: 13.320833333333333,
        18: 11.219444444444443,
        19: 10.155555555555557,
        20: 9.6125,
        21: 9.037500000000001,
        22: 8.591666666666667,
        23: 8.308333333333334,
    },
    "vlinder03": {
        0: 9.322222222222223,
        1: 8.880555555555555,
        2: 8.38888888888889,
        3: 7.933333333333334,
        4: 7.640277777777778,
        5: 7.38611111111111,
        6: 7.293055555555556,
        7: 8.093055555555557,
        8: 9.912500000000001,
        9: 12.093055555555557,
        10: 13.814492753623188,
        11: 14.773333333333333,
        12: 15.918181818181816,
        13: 16.544444444444444,
        14: 16.206944444444446,
        15: 15.593055555555557,
        16: 14.70138888888889,
        17: 13.274999999999999,
        18: 11.997222222222222,
        19: 11.069444444444445,
        20: 10.274999999999999,
        21: 9.754166666666666,
        22: 9.177777777777777,
        23: 8.684722222222222,
    },
}

stats_test = {
    "mean": {
        ("vlinder01", 0): 8.927777777777777,
        ("vlinder01", 1): 8.706944444444444,
        ("vlinder01", 2): 8.397222222222222,
        ("vlinder01", 3): 8.05,
        ("vlinder01", 4): 7.697222222222223,
        ("vlinder01", 5): 7.591666666666667,
        ("vlinder01", 6): 7.668055555555556,
        ("vlinder01", 7): 9.37361111111111,
        ("vlinder01", 8): 11.205555555555556,
        ("vlinder01", 9): 12.833333333333334,
        ("vlinder01", 10): 15.100000000000001,
        ("vlinder01", 11): 15.808333333333334,
        ("vlinder01", 12): 16.266666666666666,
        ("vlinder01", 13): 16.59027777777778,
        ("vlinder01", 14): 16.12222222222222,
        ("vlinder01", 15): 15.865277777777777,
        ("vlinder01", 16): 15.155714285714287,
        ("vlinder01", 17): 12.578333333333335,
        ("vlinder01", 18): 10.593333333333334,
        ("vlinder01", 19): 9.398333333333333,
        ("vlinder01", 20): 8.841666666666667,
        ("vlinder01", 21): 8.708571428571428,
        ("vlinder01", 22): 8.504166666666666,
        ("vlinder01", 23): 7.872222222222222,
        ("vlinder02", 0): 9.422222222222222,
        ("vlinder02", 1): 8.9875,
        ("vlinder02", 2): 8.431944444444445,
        ("vlinder02", 3): 8.013043478260869,
        ("vlinder02", 4): 8.129577464788733,
        ("vlinder02", 5): 8.005555555555555,
        ("vlinder02", 6): 7.86388888888889,
        ("vlinder02", 7): 9.452777777777778,
        ("vlinder02", 8): 11.801388888888889,
        ("vlinder02", 9): 14.016666666666667,
        ("vlinder02", 10): 15.751388888888888,
        ("vlinder02", 11): 16.479166666666668,
        ("vlinder02", 12): 17.047222222222224,
        ("vlinder02", 13): 17.109722222222224,
        ("vlinder02", 14): 16.47222222222222,
        ("vlinder02", 15): 16.140277777777776,
        ("vlinder02", 16): 15.333333333333334,
        ("vlinder02", 17): 13.320833333333333,
        ("vlinder02", 18): 11.219444444444443,
        ("vlinder02", 19): 10.155555555555557,
        ("vlinder02", 20): 9.6125,
        ("vlinder02", 21): 9.037500000000001,
        ("vlinder02", 22): 8.591666666666667,
        ("vlinder02", 23): 8.308333333333334,
        ("vlinder03", 0): 9.322222222222223,
        ("vlinder03", 1): 8.880555555555555,
        ("vlinder03", 2): 8.38888888888889,
        ("vlinder03", 3): 7.933333333333334,
        ("vlinder03", 4): 7.640277777777778,
        ("vlinder03", 5): 7.38611111111111,
        ("vlinder03", 6): 7.293055555555556,
        ("vlinder03", 7): 8.093055555555557,
        ("vlinder03", 8): 9.912500000000001,
        ("vlinder03", 9): 12.093055555555557,
        ("vlinder03", 10): 13.814492753623188,
        ("vlinder03", 11): 14.773333333333333,
        ("vlinder03", 12): 15.918181818181816,
        ("vlinder03", 13): 16.544444444444444,
        ("vlinder03", 14): 16.206944444444446,
        ("vlinder03", 15): 15.593055555555557,
        ("vlinder03", 16): 14.70138888888889,
        ("vlinder03", 17): 13.274999999999999,
        ("vlinder03", 18): 11.997222222222222,
        ("vlinder03", 19): 11.069444444444445,
        ("vlinder03", 20): 10.274999999999999,
        ("vlinder03", 21): 9.754166666666666,
        ("vlinder03", 22): 9.177777777777777,
        ("vlinder03", 23): 8.684722222222222,
    },
    "std": {
        ("vlinder01", 0): 3.463094042079706,
        ("vlinder01", 1): 3.408329746201494,
        ("vlinder01", 2): 3.4654419559302614,
        ("vlinder01", 3): 3.3482810389053306,
        ("vlinder01", 4): 3.2662582383362198,
        ("vlinder01", 5): 3.2377069785379358,
        ("vlinder01", 6): 3.1414460825912403,
        ("vlinder01", 7): 2.3032734241408854,
        ("vlinder01", 8): 2.1446880917640714,
        ("vlinder01", 9): 2.619375261500992,
        ("vlinder01", 10): 0.7492719470954619,
        ("vlinder01", 11): 0.6806686853050685,
        ("vlinder01", 12): 0.8722643351774357,
        ("vlinder01", 13): 0.8680818618862144,
        ("vlinder01", 14): 1.2548746890160087,
        ("vlinder01", 15): 1.5283388498832982,
        ("vlinder01", 16): 1.2434915446162427,
        ("vlinder01", 17): 1.5662858654886016,
        ("vlinder01", 18): 1.951085457707282,
        ("vlinder01", 19): 2.3256485639360593,
        ("vlinder01", 20): 2.668414857758692,
        ("vlinder01", 21): 2.9161093226730266,
        ("vlinder01", 22): 3.188113376498465,
        ("vlinder01", 23): 3.3695497347138965,
        ("vlinder02", 0): 2.9830139887308635,
        ("vlinder02", 1): 2.936172534938908,
        ("vlinder02", 2): 2.7991863198451328,
        ("vlinder02", 3): 2.637451591051608,
        ("vlinder02", 4): 2.3485189476529467,
        ("vlinder02", 5): 2.3652941864426524,
        ("vlinder02", 6): 2.341078694164118,
        ("vlinder02", 7): 1.9347796921003715,
        ("vlinder02", 8): 1.6911112319265236,
        ("vlinder02", 9): 0.953422301242761,
        ("vlinder02", 10): 0.7317725067020114,
        ("vlinder02", 11): 0.6844217256864004,
        ("vlinder02", 12): 0.6346724168280178,
        ("vlinder02", 13): 0.9678143987796058,
        ("vlinder02", 14): 1.5775022351915797,
        ("vlinder02", 15): 1.732510977728548,
        ("vlinder02", 16): 1.5850534269220025,
        ("vlinder02", 17): 1.798586103535493,
        ("vlinder02", 18): 1.9580558688664906,
        ("vlinder02", 19): 2.3030835966608194,
        ("vlinder02", 20): 2.4443633060228565,
        ("vlinder02", 21): 2.4320881597716406,
        ("vlinder02", 22): 2.5443530404340904,
        ("vlinder02", 23): 2.8625680384312053,
        ("vlinder03", 0): 2.497955345411165,
        ("vlinder03", 1): 2.6470626233773498,
        ("vlinder03", 2): 2.5714880634993955,
        ("vlinder03", 3): 2.6664319145497757,
        ("vlinder03", 4): 2.737515050944396,
        ("vlinder03", 5): 2.709987035754841,
        ("vlinder03", 6): 2.575962907007507,
        ("vlinder03", 7): 2.3697819706503402,
        ("vlinder03", 8): 1.9990094377952237,
        ("vlinder03", 9): 1.676038789837253,
        ("vlinder03", 10): 1.3043226056101083,
        ("vlinder03", 11): 1.1559524162786996,
        ("vlinder03", 12): 1.1226218064333688,
        ("vlinder03", 13): 1.2640075869901661,
        ("vlinder03", 14): 1.1605397468042047,
        ("vlinder03", 15): 1.3464953349384277,
        ("vlinder03", 16): 1.458872133866912,
        ("vlinder03", 17): 1.2543726337198846,
        ("vlinder03", 18): 1.3714361981865564,
        ("vlinder03", 19): 1.5316969222593748,
        ("vlinder03", 20): 1.740952467173569,
        ("vlinder03", 21): 1.8825318208035762,
        ("vlinder03", 22): 2.0071156673282258,
        ("vlinder03", 23): 2.1188396852743696,
    },
    "median": {
        ("vlinder01", 0): 9.5,
        ("vlinder01", 1): 9.9,
        ("vlinder01", 2): 9.8,
        ("vlinder01", 3): 9.5,
        ("vlinder01", 4): 9.15,
        ("vlinder01", 5): 9.15,
        ("vlinder01", 6): 9.4,
        ("vlinder01", 7): 10.55,
        ("vlinder01", 8): 12.2,
        ("vlinder01", 9): 13.7,
        ("vlinder01", 10): 15.0,
        ("vlinder01", 11): 15.8,
        ("vlinder01", 12): 16.6,
        ("vlinder01", 13): 16.75,
        ("vlinder01", 14): 16.35,
        ("vlinder01", 15): 16.3,
        ("vlinder01", 16): 15.4,
        ("vlinder01", 17): 12.649999999999999,
        ("vlinder01", 18): 10.9,
        ("vlinder01", 19): 8.3,
        ("vlinder01", 20): 7.65,
        ("vlinder01", 21): 7.05,
        ("vlinder01", 22): 8.3,
        ("vlinder01", 23): 6.75,
        ("vlinder02", 0): 10.1,
        ("vlinder02", 1): 9.8,
        ("vlinder02", 2): 8.75,
        ("vlinder02", 3): 8.4,
        ("vlinder02", 4): 9.0,
        ("vlinder02", 5): 9.0,
        ("vlinder02", 6): 8.8,
        ("vlinder02", 7): 10.2,
        ("vlinder02", 8): 12.3,
        ("vlinder02", 9): 14.25,
        ("vlinder02", 10): 15.6,
        ("vlinder02", 11): 16.6,
        ("vlinder02", 12): 17.1,
        ("vlinder02", 13): 17.25,
        ("vlinder02", 14): 16.8,
        ("vlinder02", 15): 16.55,
        ("vlinder02", 16): 15.6,
        ("vlinder02", 17): 13.2,
        ("vlinder02", 18): 11.149999999999999,
        ("vlinder02", 19): 9.7,
        ("vlinder02", 20): 9.3,
        ("vlinder02", 21): 8.65,
        ("vlinder02", 22): 8.15,
        ("vlinder02", 23): 7.949999999999999,
        ("vlinder03", 0): 8.7,
        ("vlinder03", 1): 8.25,
        ("vlinder03", 2): 8.1,
        ("vlinder03", 3): 7.5,
        ("vlinder03", 4): 7.4,
        ("vlinder03", 5): 7.3,
        ("vlinder03", 6): 7.4,
        ("vlinder03", 7): 8.4,
        ("vlinder03", 8): 10.850000000000001,
        ("vlinder03", 9): 12.649999999999999,
        ("vlinder03", 10): 14.2,
        ("vlinder03", 11): 14.5,
        ("vlinder03", 12): 15.850000000000001,
        ("vlinder03", 13): 16.7,
        ("vlinder03", 14): 16.5,
        ("vlinder03", 15): 15.65,
        ("vlinder03", 16): 14.75,
        ("vlinder03", 17): 13.2,
        ("vlinder03", 18): 11.8,
        ("vlinder03", 19): 10.75,
        ("vlinder03", 20): 9.8,
        ("vlinder03", 21): 9.05,
        ("vlinder03", 22): 8.55,
        ("vlinder03", 23): 7.949999999999999,
    },
}

assert (
    temp_diurnal.eq(pd.DataFrame(temp_diurnal_test)).all().all()
), f"Maybe something wrong with the verbose output, since it is not equal to hardcoded df."
assert (
    stats.eq(pd.DataFrame(stats_test)).all().all()
), f"Maybe something wrong with the verbose output, since it is not equal to hardcoded df."
