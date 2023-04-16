#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 16 20:20:47 2023

@author: thoverga
"""
import sys
import metobs_toolkit

is_succesfull = metobs_toolkit.launch_gui()
if is_succesfull:
    sys.exit() #normal exit
else:
    sys.exit('Problem launching GUI!')



