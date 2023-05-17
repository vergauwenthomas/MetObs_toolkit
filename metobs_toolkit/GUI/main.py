#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 09:47:48 2023

@author: thoverga
"""

# command: designer in terminal to open the desinger tool

# DEBUG
import sys
sys.path.insert(0, '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/metobs_toolkit')
import metobs_toolkit


import os, sys
from pathlib import Path
import matplotlib
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QDialog, QApplication, QFileDialog, QMainWindow
from PyQt5.uic import loadUi

import pandas as pd
import pytz

import metobs_toolkit.GUI.data_func as data_func
import metobs_toolkit.GUI.template_func as template_func
from metobs_toolkit.GUI.json_save_func import get_saved_vals, update_json_file
from metobs_toolkit.GUI.pandasmodel import DataFrameModel
import metobs_toolkit.GUI.path_handler as path_handler
import metobs_toolkit.GUI.tlk_scripts as tlk_scripts
from metobs_toolkit.GUI.errors import Error, Notification

from metobs_toolkit.GUI.extra_windows import MergeWindow, TimeSeriesWindow


class MainWindow(QDialog):
    def __init__(self):
        super(MainWindow, self).__init__()
        # Load gui widgets
        loadUi(os.path.join(path_handler.GUI_dir,'templ_build.ui'),
               self) #open the ui file

        # -------- Information to pass beween different triggers ----------
        self.dataset = None #the vlindertoolkit dataset instance
        self.merge_window = MergeWindow() #New  ui window

        self.tswindow = TimeSeriesWindow(dataset=None) #for plots

        # P1 ------------------------



        # P2 ------------------------
        self.template_dict = None #dict; all available templname : templpath
        self.default_settings = tlk_scripts.get_default_settings()





        # ------- Setup (widgets and default values) ---------------

        # Setup error message dialog
        self.error_dialog = QtWidgets.QErrorMessage(self)
        self.set_datapaths_init()

        # link dfmodels to tables
        self.templmodel = DataFrameModel()
        self.table.setModel(self.templmodel)



        # ------- Callbacks (Getters) ---------
        self.Browse_data_B.clicked.connect(lambda: self.browsefiles_data()) #browse datafile
        self.Browse_data_B_2.clicked.connect(lambda: self.browsefiles_data_p2()) #browse datafile
        self.Browse_metadata_B.clicked.connect(lambda: self.browsefiles_metadata()) #browse metadatafile
        self.Browse_metadata_B_2.clicked.connect(lambda: self.browsefiles_metadata_p2()) #browse metadatafile

        # ------- Callbacks (triggers) ------------
        # P1 ------
        # save paths when selected
        self.save_data_path.clicked.connect(lambda: self.save_path(savebool=self.save_data_path.isChecked(),
                                                                   savekey='data_file_path',
                                                                   saveval=self.data_file_T.text()))
        self.save_metadata_path.clicked.connect(lambda: self.save_path(savebool=self.save_metadata_path.isChecked(),
                                                                   savekey='metadata_file_path',
                                                                   saveval=self.metadata_file_T.text()))

        self.data_file_T.editingFinished.connect(lambda: self.link_datapath())
        self.metadata_file_T.editingFinished.connect(lambda: self.link_datapath())
        # initiate the start mapping module
        self.start_mapping_B.clicked.connect(lambda: self.set_templ_map_val())

        # construnct the mappindict
        self.build_B.clicked.connect(lambda: self.build_template())

        # save template
        self.save_template.clicked.connect(lambda: self.save_template_call())

        # P2 -------
        self.make_dataset.clicked.connect(lambda: self.make_tlk_dataset())
        self.apply_qc.clicked.connect(lambda: self.apply_qc_on_dataset())

        self.show_dataset.clicked.connect(lambda: self.create_dataset_window())


         # P3 -------
        self.plot_dataset.clicked.connect(lambda: self.make_figure())




        # ------- Initialize --------------



        # ----------P2 -----------------
        self.set_possible_templates() #load available dataset
        self.set_timezone_spinners()




        tlk_scripts.set_qc_default_settings(self)


        #----- Cleanup files ----------------
        path_handler.clear_dir(path_handler.TMP_dir) #cleanup tmp folder

#%%
# =============================================================================
# Helpers
# =============================================================================
    def get_val(self, data_func_return):

        if not data_func_return[1]:
            self.error_dialog.showMessage(data_func_return[2])
        return data_func_return[0]


# =============================================================================
# Init values
# =============================================================================
    def set_templ_map_val(self):
        # First try reading the datafile
        data_columns = self.read_datafiles(self.data_file_T.text())

        # Check if meta data is given, if so read the metadata columns
        if len(self.metadata_file_T.text()) > 2:
            metadata_columns = self.read_datafiles(self.metadata_file_T.text())
        else:
            metadata_columns = []

        # Set defaults appropriate
        template_func.set_templ_vals(self, data_columns, metadata_columns)



    def set_datapaths_init(self):
        saved_vals = get_saved_vals()

        # set datafile path
        if 'data_file_path' in saved_vals:
            self.data_file_T.setText(str(saved_vals['data_file_path']))
            self.data_file_T_2.setText(str(saved_vals['data_file_path']))

        # set metadata file path
        if 'metadata_file_path' in saved_vals:
            self.metadata_file_T.setText(str(saved_vals['metadata_file_path']))
            self.metadata_file_T_2.setText(str(saved_vals['metadata_file_path']))


    def set_possible_templates(self):
        templ_dict = template_func.get_all_templates()

        # remove .csv for presenting
        templ_names = [name.replace('.csv', '') for name in templ_dict.keys()]

        #update spinner
        self.select_temp.clear()
        self.select_temp.addItems(templ_names)

        # Store information to pass between triggers
        self.template_dict = templ_dict #dict templname : templpath

    def set_timezone_spinners(self):
        # add all common tz to options
        tzlist = pytz.common_timezones
        self.tz_selector.addItems(tzlist)

        # set default
        default = self.default_settings.time_settings['timezone']
        self.tz_selector.setCurrentText(default)


# =============================================================================
# Save values
# =============================================================================
    def save_path(self, savebool, savekey, saveval):
        if savebool:
            savedict = {str(savekey): str(saveval)}
            update_json_file(savedict)

# =============================================================================
# Triggers
# =============================================================================
    def browsefiles_data(self):
        fname=QFileDialog.getOpenFileName(self, 'Select data file', str(Path.home()))
        self.data_file_T.setText(fname[0]) #update text


    def browsefiles_data_p2(self):
        fname=QFileDialog.getOpenFileName(self, 'Select data file', str(Path.home()))
        self.data_file_T_2.setText(fname[0])


    def browsefiles_metadata(self):
        fname=QFileDialog.getOpenFileName(self, 'Select metadata file', str(Path.home()))
        self.metadata_file_T.setText(fname[0]) #update text

    def browsefiles_metadata_p2(self):
        fname=QFileDialog.getOpenFileName(self, 'Select metadata file', str(Path.home()))
        self.metadata_file_T_2.setText(fname[0]) #update text

    def link_datapath(self):
        self.data_file_T_2.setText(str(self.data_file_T.text()))
        self.metadata_file_T_2.setText(str(self.metadata_file_T.text()))


    def read_datafiles(self, filepath):
        if not path_handler.file_exist(filepath):
            Error(f'{filepath} is not a file.')
        _return = data_func.get_columns(filepath=filepath)
        columns = self.get_val(_return)

        return columns

    def build_template(self):
        df = template_func.make_template_build(self)
        if df.empty:
            Error('There are no mapped values.')
        self.templmodel.setDataFrame(df)

        self.save_template.setEnabled(True) #enable the save button


    def save_template_call(self):
        # copy the template to the templates dir of the toolkit
        templ_loc = os.path.join(path_handler.TMP_dir, 'template.csv')

        filename = str(self.templatename.text())
        if not filename.endswith('.csv'):
            filename = filename + '.csv'

        target_loc = os.path.join(path_handler.CACHE_dir, filename)

        # check if templatefile already exists.
        if path_handler.file_exist(target_loc):
            Error(f'{target_loc} already exists! Change name of the template file.')
            return
        path_handler.copy_file(templ_loc, target_loc)

        self.set_possible_templates() #update widget

        Notification(f'Template ({filename}) is saved!')



    def make_tlk_dataset(self):
        self.dataset, self.merge_window.comb_df = tlk_scripts.load_dataset(self)

        # trigger update in seperate window
        self.merge_window.trigger_update()


    def apply_qc_on_dataset(self):
        self.merge_window.comb_df = tlk_scripts.apply_qualitycontrol(self)

        # trigger update in seperate window
        self.merge_window.trigger_update()


    def create_dataset_window(self):

        self.merge_window.show()

# =============================================================================
# Testing
# =============================================================================


    def make_figure(self):
        self.tswindow.set_dataset(self.dataset)
        print(self.dataset)
        print('tswindow set')
        self.tswindow.make_plot()
        print('plot made')
        self.tswindow.show()


#%%

# =============================================================================
# Main and protector
# =============================================================================

def main():
    try:
        app=QApplication(sys.argv)

        mainwindow = MainWindow()
        widget = QtWidgets.QStackedWidget()
        widget.addWidget(mainwindow)
        widget.show()

        succesfull=True
        print('einde succes')
    except:
        print('begin fail')
        # sys.exit('Something went wrong in the GUI')
        succesfull=False
        pass
    sys.exit(app.exec_())

    return succesfull


if __name__ == '__main__':
    matplotlib.use('Qt5Agg') #in protector because server runners do not support this, when this module is imported from the __init__
    app=QApplication(sys.argv)

    widget = main()
    widget.show()
    sys.exit(app.exec_())



