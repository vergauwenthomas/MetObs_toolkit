#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 14:05:22 2023

@author: thoverga
"""

import os
import shutil
from pathlib import Path



# =============================================================================
# Main folders
# =============================================================================

GUI_dir = str( Path(__file__).resolve().parents[0])

TLK_dir = str( Path(__file__).resolve().parents[1])


# =============================================================================
# Derived locations
# =============================================================================

TMP_dir = os.path.join(GUI_dir, 'tmp')

CACHE_dir = os.path.join(GUI_dir, 'cache')
# toolkit location of templates
tlk_default_template = os.path.join(TLK_dir, 'data_templates',
                                    'template_defaults', 'default_template.csv')


# =============================================================================
# Helper functions
# =============================================================================

def make_dir(dir_path):
    if os.path.isdir(dir_path):
        return
    else:
        os.mkdir(dir_path)
        return


def file_exist(filepath):
    if os.path.isfile(filepath):
        return True
    elif os.path.islink(filepath):
        return True
    else:
        return False


def copy_file(filepath, targetpath):
    shutil.copyfile(filepath, targetpath)

def list_csv_filenames(folderpath):
    all_stuff = os.listdir(folderpath)
    all_files = [file for file in all_stuff if os.path.isfile(os.path.join(
                                                folderpath, file))]

    filenames =[file for file in all_files if file.endswith('.csv')]
    filepaths =[os.path.join(folderpath, file) for file in filenames]

    return filenames, filepaths


def clear_dir(directory):
    make_dir(directory)
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                shutil.rmtree(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))

