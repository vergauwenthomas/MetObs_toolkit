#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:25:02 2022
@author: thoverga
"""

#%%
import metobs_toolkit
import os
import sys
from pathlib import Path


lib_folder = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(lib_folder))


import metobs_toolkit


#%%
metobs_toolkit.launch_gui()

#%%

# from PyQt5.QtWidgets import QApplication, QMainWindow, QWidget, QVBoxLayout, QPushButton


# class MainWindow(QMainWindow):
#     def __init__(self):
#         super().__init__()
#         self.setGeometry(100, 100, 300, 200)
#         self.button = QPushButton("Show Layout", self)
#         self.button.clicked.connect(self.show_layout)

#     def show_layout(self):
#         layout_window = QWidget(self)
#         layout_window.setGeometry(200, 200, 300, 200)
#         layout = QVBoxLayout(layout_window)
#         layout.addWidget(QPushButton("Button 1", layout_window))
#         layout.addWidget(QPushButton("Button 2", layout_window))
#         layout.addWidget(QPushButton("Button 3", layout_window))
#         layout_window.show()


# if __name__ == '__main__':
#     app = QApplication([])
#     window = MainWindow()
#     window.show()
#     app.exec_()





#%% % Import


# # testdatafile = os.path.join(str(lib_folder), 'tests', 'test_data',  'testdata_okt_small.csv')
# testdatafile = os.path.join(str(lib_folder), 'tests', 'test_data',  'testdata_breaking.csv')


# template = os.path.join(str(lib_folder), 'tests', 'test_data',  'template_breaking.csv')

# static_data = os.path.join(
#     str(lib_folder), 'static_data', 'vlinder_metadata.csv')



# # #% Setup dataset

# dataset = metobs_toolkit.Dataset()
# dataset.update_settings(input_data_file=testdatafile,
#                         # input_metadata_file=static_data,
#                         data_template_file= template,
#                         output_folder='/home/thoverga/Documents'
#                         )


# dataset.import_data_from_file(coarsen_timeres=True)
# # dataset.apply_quality_control()



# dataset.import_data_from_file()
# dataset.apply_quality_control(gross_value=True,
#                               persistance=False)


# dataset.get_qc_stats('temp', make_plot=True)

