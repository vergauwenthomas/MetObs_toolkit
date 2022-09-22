#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 13:53:08 2022

@author: thoverga
"""
from pathlib import Path
import os, sys

lib_folder = Path(__file__).resolve().parents[1]
sys.path.append(str(lib_folder))
from src import vlinder_toolkit


#%%
settings = vlinder_toolkit.Settings()


settings.check_settings()

testvariable='Testtest'
settings.update_settings(output_data_folder=testvariable)
testvariable='changed test'

if 'Testtest'!=settings.output_data_folder:
    sys.exit()

print('Settings test passed!')



