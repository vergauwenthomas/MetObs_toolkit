#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 12:41:38 2023

@author: thoverga
"""

import sys, os
import matplotlib
# matplotlib.use('Qt5Agg')

from PyQt5 import QtCore, QtGui, QtWidgets

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure


from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QLabel, QVBoxLayout, QWidget

import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from PyQt5.QtWidgets import QWidget, QVBoxLayout


import sys
from PyQt5.QtWidgets import QApplication, QVBoxLayout, QWidget, QPushButton, QMainWindow
from PyQt5.uic import loadUi





# MatplotlibCanvas creator for dataset plots
class MplCanvas(FigureCanvasQTAgg):

    def __init__(self, dataset=None, width=5, height=4, dpi=100):
        # Data object
        self.dataset = dataset

        # Figure object
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)


    def set_dataset(self, dataset):
        self.dataset = dataset

    def create_toolbar(self):
        # create toolbar connectd to the canvas
        return NavigationToolbar(canvas=self)
    def _clear_axis(self):
        self.axes.cla()
    # =============================================================================
    #     fill axes methods
    # =============================================================================
    def timeseriesplot(self, obstype='temp', colorby='name',
                       stationnames=None, show_outliers=True):


        self.axes = self.dataset.make_plot(stationnames=stationnames,
                                           obstype=obstype,
                                           colorby=colorby,
                                           starttime=None,
                                           endtime=None,
                                           _ax=self.axes,
                                           title=None,
                                           legend=False,
                                           show_outliers=show_outliers)
        print('figure is made')




