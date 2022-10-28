#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 08:18:09 2022

@author: thoverga
"""
import sys

from pathlib import Path


lib_folder = Path(__file__).resolve().parents[2]
sys.path.append(str(lib_folder))
print(str(lib_folder))

from src import vlinder_toolkit



