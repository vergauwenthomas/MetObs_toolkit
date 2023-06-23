#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Idea is to test all testfiles with all kind of datasets

@author: thoverga
"""


import sys, os

from pathlib import Path
lib_folder = Path(__file__).resolve().parents[1]

# print(str(lib_folder))

sys.path.insert(0, str(lib_folder))
import metobs_toolkit

# %%
testfolder=os.path.join(str(lib_folder), 'tests', 'push_test')

from tests.push_test.test_data_paths import testdata




# %%
def read_in_the_dataset(dataname, testdatadict):
    print(f'\n ------ read dataset ({dataname}) ---------\n')
    datafile = testdatadict[dataname]['datafile']
    metafile =testdatadict[dataname]['metadatafile']
    template = testdatadict[dataname]['template']
    kwargsdict = testdatadict[dataname]['kwargs']




    dataset = metobs_toolkit.Dataset()
    dataset.update_settings(input_data_file=datafile,
                            input_metadata_file=metafile,
                            data_template_file=template,
                            metadata_template_file=template)

    dataset.import_data_from_file(**kwargsdict)
    return dataset


def IO_test(dataset, name):
    print(f'\n ------ IO tests ({name}) ---------\n')
    def del_file(file_path):
        if os.path.isfile(file_path):
            os.remove(file_path)
            print(f"{file_path} deleted.")
        else:
            print(f"{file_path} not found.")


    # Sycnronize data
    test = dataset.sync_observations(tollerance='5T', verbose=True)

    # pickel test
    outfolder =os.path.join(str(lib_folder), 'tests', 'test_data')
    file='dataset_IO_test'


    del_file(os.path.join(outfolder, file+'.pkl'))


    # save dataset as pickle
    dataset.update_default_name('this_is_a_test_name')

    dataset.save_dataset(outputfolder=outfolder,
                         filename=file)


    del dataset #remove from kernel


    # read dataset
    new_dataset = metobs_toolkit.Dataset()
    new_dataset = new_dataset.import_dataset(folder_path=outfolder,
                               filename=file +'.pkl')

    del_file(os.path.join(outfolder, file+'.pkl'))

    assert new_dataset.settings.app["default_name"] == 'this_is_a_test_name', 'some attributes are not correctly saved when pickled.'




def qc_testing(dataset, name):
    print(f'\n ------ QC tests ({name}) ---------\n')


    # on get station
    stationname = dataset.metadf.index[0]
    station = dataset.get_station(stationname)
    station.apply_quality_control(obstype='temp')
    station.get_qc_stats(make_plot=False)


    #on dataset
    dataset.get_qc_stats(make_plot=False)
    dataset.apply_quality_control(obstype="temp")
    dataset.get_qc_stats(make_plot=True)


    # titan test
    dataset.update_titan_qc_settings(obstype='temp',
                                  buddy_radius=50000,
                                  buddy_num_min=3,
                                  buddy_max_elev_diff=200,
                                  buddy_threshold=3)


def gapfill_testing(dataset, name):
    print(f'\n ------ gaps missing tests ({name})---------\n')

    # testing conversion to df + update from outliers
    _ = dataset.get_gaps_df()

    init_outl_shape = dataset.outliersdf.shape
    dataset.update_gaps_and_missing_from_outliers(n_gapsize=3)
    assert init_outl_shape!=dataset.outliersdf.shape, 'outliers still the same as before updateing to gaps'


    _ = dataset.get_gaps_df()

    dataset.get_gaps_info()
    dataset.get_missing_obs_info()
    # filling
    datset.fill_missing_obs_linear()
    datset.fill_gaps_linear()

    dataset.get_gaps_info()
    dataset.get_missing_obs_info()

    # plot test
    dataset.make_plot(colorby='label', title='AFTER GAP AND MISSING FILL')




def plot_testing(dataset, name):
    print(f'\n ------ plot tests ({name})---------\n')

    dataset.make_plot(colorby='name')
    dataset.make_plot(colorby='label')

    if not dataset.metadf.empty:
        dataset.make_geo_plot(variable='temp')


# %%



for name in testdata:
    print(f'\n ************ {name} *************\n')
    dataset = read_in_the_dataset(name, testdata)
    datset = dataset.coarsen_time_resolution(freq=testdata[name]['coarsen'])
    IO_test(dataset, name)
    qc_testing(dataset, name)
    plot_testing(dataset, name)
    gapfill_testing(dataset, name)





