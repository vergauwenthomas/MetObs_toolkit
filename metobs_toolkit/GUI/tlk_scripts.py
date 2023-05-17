#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 31 09:40:08 2023

@author: thoverga
"""

#%% Load vlinder toolkit (not shure what is best)
import sys
from io import StringIO

from PyQt5.QtWidgets import QSpinBox, QDoubleSpinBox


import metobs_toolkit.GUI.path_handler as path_handler
from metobs_toolkit.GUI.errors import Error, Notification
# method1: add the toolkit to path,
# sys.path.append(path_handler.TLK_dir)
# from vlinder_toolkit import Dataset

# method2: loead as package
import metobs_toolkit

#%% Mapping
# when specific manipulation has to be done on values

tlk_to_gui = {
    'window_variation__max_increase_per_second': {'multiply': 3600.0},
    'window_variation__max_decrease_per_second' : {'multiply': 3600.0},
    'step__max_increase_per_second' : {'multiply': 3600.0},
    'step__max_decrease_per_second' : {'multiply': 3600.0},

    'persistance__time_window_to_check': {'remove': 'h'},
    'window_variation__time_window_to_check': {'remove': 'h'},
}


gui_to_tlk = {
    'window_variation__max_increase_per_second': {'devide': 3600.0},
    'window_variation__max_decrease_per_second' : {'devide': 3600.0},
    'step__max_increase_per_second' : {'devide': 3600.0},
    'step__max_decrease_per_second' : {'devide': 3600.0},

    'persistance__time_window_to_check': {'append': 'h'},
    'window_variation__time_window_to_check': {'append': 'h'},
}

def _append(x, val):
    return str(x)+str(val)
def _div(x, val):
    return float(x)/float(val)
def _multiply(x, val):
    return float(x) * float(val)
def _remove(x, val):
    return str(x).replace(str(val), '')


qc_not_in_gui = ['duplicated_timestamp', 'internal_consistency'] #not present in gui

#%% Helpers

class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self
    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        del self._stringio    # free up some memory
        sys.stdout = self._stdout

#%% Initialisation of widgets


def get_default_settings():
    "Get default settings for initiating the widgets"
    _dummy = metobs_toolkit.Dataset()
    return _dummy.settings



def set_qc_default_settings(main):
    # set default settings

    main.obstype_selector.addItems(['temp']) #TODO: for now only QC on temp
    obstype = main.obstype_selector.currentText()

    # Get standard tlk settings
    qc_set = main.default_settings.qc['qc_check_settings']

    for checkname in qc_set:
        if checkname in qc_not_in_gui:
            continue
        check_set = qc_set[checkname][obstype]
        for set_name, set_val in check_set.items():
            widgetname = checkname+'__'+set_name
            if widgetname in tlk_to_gui:
                apl=list(tlk_to_gui[widgetname].keys())[0]
                val =list(tlk_to_gui[widgetname].values())[0]

                if apl == 'multiply': set_val = _multiply(set_val, val)
                elif apl == 'devide': set_val = _div(set_val, val)
                elif apl == 'append': set_val = _append(set_val, val)
                elif apl == 'remove': set_val = _remove(set_val, val)



            print(f' {widgetname} : {set_val}')

            widg = getattr(main, widgetname)
            if isinstance(widg, type(QSpinBox())):
                print('spinbox')
                widg.setValue(int(set_val))
            elif isinstance(widg, type(QDoubleSpinBox())):
                widg.setValue(float(set_val))
                print('doublspinbox')


def apply_qualitycontrol(main):
    # check if dataset exist
    if isinstance(main.dataset, type(None)):
        Error('Cannot apply quality control, first make a Dataset.')

    with Capturing() as terminaloutput:
        main.dataset.apply_quality_control()
        print('QC IS DONE')

    # make mergedf for visualising
    with Capturing() as terminaloutput:
        comb_df = main.dataset.combine_all_to_obsspace()

    # write terminal output to
    main.prompt.append('\n \n ------ Qualtiy control ----- \n \n')
    for line in terminaloutput:
        main.prompt.append(line + '\n')

    return comb_df


def update_qc_settings(main, dataset, obstype):
    qc_set = dataset.settings.qc['qc_check_settings']

    for checkname in qc_set:
        if checkname in qc_not_in_gui:
            continue
        check_set = qc_set[checkname][obstype]
        for set_name, _ in check_set.items():
            widgetname = checkname+'__'+set_name

            widget = getattr(main, widgetname)
            widget_val = widget.value()

            if widgetname in gui_to_tlk:
                apl=list(gui_to_tlk[widgetname].keys())[0]
                val =list(gui_to_tlk[widgetname].values())[0]

                if apl == 'multiply': set_val = _multiply(widget_val, val)
                elif apl == 'devide': set_val = _div(widget_val, val)
                elif apl == 'append': set_val = _append(widget_val, val)
                elif apl == 'remove': set_val = _remove(widget_val, val)
            else:
                set_val = widget_val
            #set value
            dataset.settings.qc['qc_check_settings'][checkname][obstype][set_name] = set_val








#%% Toolkit functions

def load_dataset(main):


    print('in the tlk module to load dataset')
    # Get IO information

    data_file = main.data_file_T_2.text()
    metadata_file = main.metadata_file_T_2.text()
    template_name =str(main.select_temp.currentText())+'.csv'

    template_file = main.template_dict[template_name]

    # Basic checks
    for file in [data_file, metadata_file, template_file]:
        if not path_handler.file_exist(file):
            Error(f'{file} does not exist!')
            return


    test=Capturing()

    with Capturing() as terminaloutput:

        # init dataset
        dataset = metobs_toolkit.Dataset()

        #use metadata
        use_metadata = main.metadata_box.isChecked()
        if use_metadata:
            dataset.update_settings(input_data_file=data_file,
                                    input_metadata_file=metadata_file,
                                    data_template_file=template_file,
                                    metadata_template_file=template_file,
                                    )
        else:
            dataset.update_settings(input_data_file=data_file,
                                    data_template_file=template_file,
                                    )

        # set timezone
        tz_select = str(main.tz_selector.currentText())
        dataset.update_timezone(tz_select)

        # create dataset

        dataset.import_data_from_file()

        # coarsen settings
        use_coarsening = main.resample_box.isChecked()
        resolution = str(main.resample_spinbox.value())+'T'
        if use_coarsening:
            dataset.coarsen_time_resolution(freq=resolution)

        # Update QC settings
        update_qc_settings(main, dataset, main.obstype_selector.currentText())

        # make mergedf (to visualise)
        comb_df= dataset.combine_all_to_obsspace()

        print('DATASET IS LOADED')




    # write terminal output to
    main.prompt.append('------ Make Dataset ----- \n \n')
    for line in terminaloutput:
        main.prompt.append(line + '\n')

    return dataset, comb_df


