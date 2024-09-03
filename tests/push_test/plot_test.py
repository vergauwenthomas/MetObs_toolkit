# -*- coding: utf-8 -*-

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 08:29:37 2022

@author: thoverga
"""

import sys, os

from pathlib import Path


lib_folder = Path(__file__).resolve().parents[2]
# sys.path.append(str(lib_folder))


import metobs_toolkit

# %% IO testdata

testdatafile = os.path.join(
    str(lib_folder), "tests", "test_data", "vlinderdata_small.csv"
)
metadatafile = os.path.join(str(lib_folder), "static_data", "vlinder_metadata.csv")


dataset = metobs_toolkit.Dataset()
dataset.update_file_paths(
    input_data_file=testdatafile,
    template_file=metobs_toolkit.demo_template,
    input_metadata_file=metadatafile,
)


dataset.import_data_from_file()
dataset.coarsen_time_resolution()

# %% timeseries plots of dataset
dataset.make_plot()
dataset.make_plot(
    obstype="humidity", stationnames=["vlinder02", "vlinder17"], legend=False
)

from datetime import datetime

dataset.make_plot(
    starttime=datetime(2022, 9, 4), endtime=datetime(2022, 9, 6), title="test"
)

# %% timeseries plot of station

dataset.get_station("vlinder05").make_plot()

# %% Make spatial plot
dataset.make_geo_plot()
dataset.make_geo_plot(
    variable="wind_direction", timeinstance=datetime(2022, 9, 5, 12, 0)
)

dataset.make_geo_plot(variable="temp", timeinstance=datetime(2022, 9, 5, 12, 0))

# %% Interactive spatial plot
outfolder = os.path.join(str(lib_folder), "development")
outfile = "deletame"


outfile = os.path.join(outfolder, outfile)
dataset.make_interactive_plot(outputfile=outfile, obstype="humidity", radius=11)

assert os.path.exists(outfile + ".html"), "interactive html is not saved!"
os.remove(outfile + ".html")

# %% GEE plots


# interactive gee static model data

lcz_map = dataset.gee_datasets["lcz"]
lcz_map.set_metadf(dataset.metadf)


dataset.make_gee_static_spatialplot(
    Model="lcz", outputfolder=outfolder, filename=outfile, overwrite=True
)

trgpath = os.path.join(outfolder, outfile + ".html")

assert os.path.exists(trgpath), "interactive static geeplot is not saved!"
# os.remove(trgpath)

# interactive geedynamic model data

# plot ERA5 'temperature' on a specific date

import datetime

inst = datetime.datetime(2004, 9, 16, 22, 18)

dataset.make_gee_dynamic_spatialplot(
    timeinstance=inst, outputfolder=outfolder, filename=outfile, overwrite=True
)
