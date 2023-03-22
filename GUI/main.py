#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 09:47:48 2023

@author: thoverga
"""

# command: designer in terminal to open the desinger tool



import sys
from pathlib import Path
from PyQt5 import QtWidgets
from PyQt5.QtWidgets import QDialog, QApplication, QFileDialog, QMainWindow
from PyQt5.uic import loadUi

import data_func
import template_func






class MainWindow(QDialog):
    def __init__(self):
        super(MainWindow, self).__init__()
        # Load gui widgets
        loadUi('templ_build.ui', self) #open the ui file
        
        # Setup error message dialog
        self.error_dialog = QtWidgets.QErrorMessage(self)
        
        
        # Callbacks (Getters)
        self.Browse_data_B.clicked.connect(lambda: self.browsefiles_data()) #browse datafile
        self.Browse_metadata_B.clicked.connect(lambda: self.browsefiles_metadata()) #browse metadatafile

        # self.start_mapping_B.clicked.connect(lambda: self.read_datafiles())
        # debug
        self.start_mapping_B.clicked.connect(lambda: self.set_templ_map_val())
        default_path = '/home/thoverga/Documents/VLINDER_github/vlinder_toolkit/tests/test_data/vlinderdata_small.csv'
        self.data_file_T.setText(default_path)
        
        self.build_B.clicked.connect(lambda: self.build_template())


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
        print(self.data_file_T.text())
        # First try reading the datafile
        data_columns = self.read_datafiles(self.data_file_T.text())
        
        
        # Check if meta data is given, if so read the metadata columns
        if len(self.metadata_file_T.text()) > 2:
            metadata_columns = self.read_datafiles(self.data_file_T.text())
        else:
            metadata_columns = []
        
        # Set defaults appropriate
        template_func.set_templ_vals(self, data_columns, metadata_columns)
        

# =============================================================================
# Triggers
# =============================================================================
    def browsefiles_data(self):
        fname=QFileDialog.getOpenFileName(self, 'Select data file', str(Path.home()))
        self.data_file_T.setText(fname[0]) #update text
        

    def browsefiles_metadata(self):
        fname=QFileDialog.getOpenFileName(self, 'Select metadata file', str(Path.home()))
        self.metadata_file_T.setText(fname[0]) #update text

    def read_datafiles(self, filepath):
        _return = data_func.get_columns(filepath=filepath)
        print(_return)
        columns = self.get_val(_return)
        
        return columns

    def build_template(self):
        template_func.make_template_build(self)




# =============================================================================
# Main and protector
# =============================================================================

def main():
    mainwindow = MainWindow()

    widget = QtWidgets.QStackedWidget()
    widget.addWidget(mainwindow)
    # widget.show()
    
    return widget
    
    
    
    
    
if __name__ == '__main__':
    app=QApplication(sys.argv)
    
    widget = main()
    widget.show()
    sys.exit(app.exec_())



