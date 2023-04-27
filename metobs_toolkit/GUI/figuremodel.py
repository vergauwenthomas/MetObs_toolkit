#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 12:41:38 2023

@author: thoverga
"""

import sys, os
import matplotlib
matplotlib.use('Qt5Agg')

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

class NewWindow(QMainWindow):
    """ Creates new window """

    def __init__(self, ui_templ_file=None):
        super().__init__()

        if not isinstance(ui_templ_file, type(None)):
            loadUi(ui_templ_file,self)






# def create_window():
#     # create a new QWindow
#     window = QMainWindow()
#     window.setGeometry(100, 100, 300, 200)
#     window.setWindowTitle("Vertical Layout Example")

#     # create a QVBoxLayout widget to hold the two child widgets
#     layout = QVBoxLayout()

#     # create the two child widgets to be added to the layout
#     button1 = QPushButton("Button 1")
#     button2 = QPushButton("Button 2")

#     # add the child widgets to the layout in vertical alignment
#     layout.addWidget(button1)
#     layout.addWidget(button2)

#     # create a QWidget to add the layout to
#     widget = QWidget()
#     widget.setLayout(layout)

#     # add the QWidget to the main window as the central widget
#     window.setCentralWidget(widget)

#     # show the window
#     window.show()
#     return window




# class test(QWidget):
#     def __init__(self, parent=None):
#         super().__init__(parent)

#         # create a Figure object
#         self.figure = plt.Figure()

#         # create a Canvas widget to display the Figure
#         self.canvas = FigureCanvas(self.figure)

#         # create a Vertical Box layout to hold the Canvas widget
#         self.vertical_layout = QVBoxLayout()
#         self.vertical_layout.addWidget(self.canvas)

#         # set the layout to the Matplotlib widget
#         self.setLayout(self.vertical_layout)

#         # create a subplot on the Figure
#         ax = self.figure.add_subplot(111)

#         # plot some data on the subplot
#         x = [1, 2, 3, 4, 5]
#         y = [2, 4, 6, 8, 10]
#         ax.plot(x, y)

#         # set the title of the subplot
#         ax.set_title('Simple Matplotlib Plot')
#         print('ax gemaakt')
#         # refresh the Canvas
#         self.canvas.draw()
#         print('ax drawn')
# # MatplotlibCanvas
# class MplCanvas(FigureCanvasQTAgg):

#     def __init__(self, parent=None, width=5, height=4, dpi=100):
#         self.parent = parent
#         self.widget = QtWidgets.QWidget() #to call with show
#         self.layouts = QtWidgets.QVBoxLayout() # add visual children here

#         fig = Figure(figsize=(width, height), dpi=dpi)
#         self.axes = fig.add_subplot(111)
#         super(MplCanvas, self).__init__(fig)

#     def create_toolbar(self):
#         # create toolbar
#         toolbar = NavigationToolbar(canvas=self,
#                                     parent=self.parent)

#         # add it to the layouts
#         self.layouts.addWidget(toolbar)

#     def testplot(self):
#         self.axes.plot([0,1,2,3,4], [10,1,20,3,40])

#         # add it to the layouts
#         self.layouts.addWidget(toolbar)


# class AnotherWindow(QWidget):
#     """
#     This "window" is a QWidget. If it has no parent, it
#     will appear as a free-floating window as we want.
#     """
#     def __init__(self):
#         super().__init__()
#         self.layout = QVBoxLayout()


#         # self.label = QLabel("Another Window")
#         # self.layout.addWidget(self.label)



#         # self.setLayout(self.layout)
#     def get_plot(self):
#         plot=MplCanvas()
#         plot.axes.plot([0,1,2,3,4], [10,1,20,3,40])
#         self.layout.addWidget(plot)

#     def set_layout(self):
#         self.setLayout(self.layout)


# class NewWindow(QtWidgets.QMainWindow):

#     def __init__(self, *args, **kwargs):
#         super(MainWindow, self).__init__(*args, **kwargs)

#         sc = MplCanvas(self, width=5, height=4, dpi=100)
#         sc.axes.plot([0,1,2,3,4], [10,1,20,3,40])

#         # Create toolbar, passing canvas as first parament, parent (self, the MainWindow) as second.
#         toolbar = NavigationToolbar(sc, self)
#         print('toolbar: ',type(toolbar))
#         print('sc: ', type(sc))

#         layout = QtWidgets.QVBoxLayout()
#         layout.addWidget(toolbar)
#         layout.addWidget(sc)

#         # Create a placeholder widget to hold our toolbar and canvas.
#         widget = QtWidgets.QWidget()
#         widget.setLayout(layout)
#         self.setCentralWidget(widget)

#         self.show()


# app = QtWidgets.QApplication(sys.argv)
# w = MainWindow()
# app.exec_()