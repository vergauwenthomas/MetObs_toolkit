#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 27 10:30:10 2023

@author: thoverga
"""


from bokeh.models import Button, Div, TextInput, Select
from bokeh.layouts import layout, row, column



def make_data_template_spinners(data_file_columns, tlk_meta_names, tlk_dt_names, tlk_obs_names):
    
    data_file_columns.insert(0, 'No Mapping')
    
    #get all toolkit variables
    widget_grid = [] #will be a nested list if widgets
    spinners_dict = {} #save all selected mappings
    
    
    #metadata widgest
    if not isinstance(tlk_meta_names,type(None)):
        for key, val in tlk_meta_names.items():
            spinners_dict[key] = {}
            row_widgets = []
            row_widgets.append(Div(text=text_increase(key,1), sizing_mode='stretch_width'))
            for topic, topictext in val.items():
                spinners_dict[key][topic] = TextInput(title=topic,
                                                      value=topictext,
                                                      sizing_mode='stretch_width',
                                                                    )
                row_widgets.append(spinners_dict[key][topic])
            #add selector
            spinners_dict[key]['template_column_name'] = Select(options=data_file_columns, value=data_file_columns[0])
            row_widgets.append(spinners_dict[key]['template_column_name'])
            
            widget_grid.append(row(row_widgets))
    
    #datetime widgest
    if not isinstance(tlk_dt_names,type(None)):
        for key, val in tlk_dt_names.items():
            spinners_dict[key] = {}
            row_widgets = []
            row_widgets.append(Div(text=text_increase(key,1), sizing_mode='stretch_width'))
            for topic, topictext in val.items():
                spinners_dict[key][topic] = TextInput(title=topic,
                                                                     value=topictext,
                                                                     sizing_mode='stretch_width',
                                                                    )
                row_widgets.append(spinners_dict[key][topic])
            #add selector
            spinners_dict[key]['template_column_name'] = Select(options=data_file_columns, value=data_file_columns[0])
            row_widgets.append(spinners_dict[key]['template_column_name'])
            
            widget_grid.append(row(row_widgets))
    
    
    
    #obstype widgest
    if not isinstance(tlk_obs_names,type(None)):
        for key, val in tlk_obs_names.items():
            spinners_dict[key] = {}
            row_widgets = []
            row_widgets.append(Div(text=text_increase(key,1), sizing_mode='stretch_width'))
            for topic, topictext in val.items():
                if isinstance(topictext, str):
                    spinners_dict[key][topic] = TextInput(title=topic,
                                                                         value=topictext,
                                                                         sizing_mode='stretch_width',
                                                                         )
                elif isinstance(topictext, list):
                    spinners_dict[key][topic] = Select(title=topic,
                                                                         options=topictext,
                                                                         sizing_mode='stretch_width',
                                                                         value=topictext[0]
                                                                         )
                    
                row_widgets.append(spinners_dict[key][topic])
            #add selector
            spinners_dict[key]['template_column_name'] = Select(options=data_file_columns, value=data_file_columns[0])
            row_widgets.append(spinners_dict[key]['template_column_name'])
            widget_grid.append(row(row_widgets))
        
        

    return widget_grid, spinners_dict

        
        
        
        
        
        
        
        
def text_increase(string, num):
    formatstring ='<font size="+'+str(int(num))+'">' + string + '</font>'
    return formatstring