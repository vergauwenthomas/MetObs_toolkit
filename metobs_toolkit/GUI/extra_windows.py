#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 16:17:53 2023

@author: thoverga
"""

# DEBUG
import sys
sys.path.insert(0, '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/metobs_toolkit')
import metobs_toolkit
# END DEBUG

import sys, os
import pandas as pd
from PyQt5.QtWidgets import QApplication, QVBoxLayout, QWidget, QPushButton, QMainWindow
from PyQt5.uic import loadUi

import metobs_toolkit.GUI.path_handler as path_handler
from metobs_toolkit.GUI.pandasmodel import DataFrameModel
import metobs_toolkit.GUI.template_func as template_func






# =============================================================================
# Figure windows
# =============================================================================

from metobs_toolkit.GUI.figuremodel import MplCanvas

class TimeSeriesWindow(QMainWindow):
    """ Creates new window """

    def __init__(self, dataset):
        super().__init__()
        loadUi(os.path.join(path_handler.GUI_dir,'fig_window.ui'), self)


        self.canvas=MplCanvas(dataset)

        # triggers
        self.update_plot_box.clicked.connect(lambda: self.update_plot())



    def set_dataset(self, dataset):
        self.canvas.set_dataset(dataset)
        self.init_widgets(dataset)

    def init_widgets(self, dataset):
        self.select_colorby.clear()
        self.select_colorby.addItems(['name','label'])

        stationnames = dataset.df.index.get_level_values('name').unique().to_list()
        stationnames.insert(0, 'ALL')
        self.select_subset.clear()
        self.select_subset.addItems(stationnames)

        self.select_obstype.clear()
        self.select_obstype.addItems(dataset.df.columns.to_list())
    def _make_axes(self):
        # self.canvas.testplot() #create mpl axes plot
        self.canvas.timeseriesplot() #create mpl axes plot



        # self.toolbar = self.canvas.create_toolbar()


    def _pass_to_layout(self):
        self.vert_layout.addWidget(self.canvas.create_toolbar())
        self.vert_layout.addWidget(self.canvas)

    def update_plot(self):
        obstype =self.select_obstype.currentText()
        subset=self.select_subset.currentText()
        colorby = self.select_colorby.currentText()
        show_outliers = self.select_show_outliers.isChecked()

        if subset=='ALL':
            stationnames=None
        else:
            stationnames = [subset]

        self.canvas._clear_axis()
        self.canvas.timeseriesplot(obstype=obstype,
                                   colorby=colorby,
                                   stationnames=stationnames,
                                   show_outliers=show_outliers)
        self.canvas.draw()


    def make_plot(self):
        self._make_axes()
        self._pass_to_layout()






# =============================================================================
# Table windows
# =============================================================================

class MergeWindow(QMainWindow):
    """ Creates new window """

    def __init__(self):
        super().__init__()
        loadUi(os.path.join(path_handler.GUI_dir,'merge_overview.ui'), self)

        # Define data attributes
        self.comb_df = pd.DataFrame()

        self.combmodel = DataFrameModel()
        self.merge_table.setModel(self.combmodel)


        # initialise values
        self.set_obstype_subsetting_spinner()


        # Triggers
        self.subset_merged_obstype.activated.connect(lambda: self.subset_comb_table())

    def trigger_update(self):
        self.subset_comb_table()

    def set_obstype_subsetting_spinner(self):
        # check if dataset is available, if not use all possible obstypes
        all_obstypes = list(template_func.Obs_map_values.keys())
        if self.comb_df.empty:
            obstypes = all_obstypes
        else:
            obstypes = [obstype for obstype in all_obstypes if obstype in self.comb_df.columns]

        # insert 'NO SELECTION'
        if not 'NO SELECTION' in obstypes:
            obstypes.insert(0,'NO SELECTION')

        # Update spinner
        self.subset_merged_obstype.clear()
        self.subset_merged_obstype.addItems(obstypes)
        self.subset_merged_obstype.setCurrentText('NO SELECTION')

    def subset_comb_table(self):
        # only if the combdf dataframe is not empty, and
        #  an obstype is selected, and the obstype is in the dataframe
        obstype = self.subset_merged_obstype.currentText()
        comb_df = self.comb_df.reset_index()

        if ((not comb_df.empty) &
            (obstype != 'NO SELECTION') &
            (obstype in comb_df.columns)):

            subsetcols = ['name', 'datetime', obstype, obstype+'_final_label']
            # update model
            self.combmodel.setDataFrame(comb_df[subsetcols])

        else:
            self.combmodel.setDataFrame(comb_df)