#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 09:48:45 2023

@author: thoverga
"""

import sys, os

from pathlib import Path


lib_folder = Path(__file__).resolve().parents[2]

import metobs_toolkit


print(f'Succesfull import of the metobs_toolkit version: {metobs_toolkit.__version__}')