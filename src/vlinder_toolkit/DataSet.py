#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 08:32:35 2022

@author: thoverga
"""
import pandas as pd
import Station




class Dataset:
    def __init__(self, observations):
        self.df = observations
        
    
    
    def print_dataset(self):
        print(self.df.head())
    