#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 11:53:45 2023

@author: thoverga
"""

import pandas as pd
import numpy as np
import bokeh


from bokeh.models import Button, Div, TextInput, Select
from bokeh.layouts import layout, row, column
from bokeh.models import ColumnDataSource, DataTable, DateFormatter, TableColumn

#%%
import gui_settings

from gui_modules.set_widget_values import *
from gui_modules.data_handlers import read_csv, init_datatable_with_default

#%%
import os, sys
from pathlib import Path
toolkit_folder = Path(__file__).resolve().parents[2]
sys.path.append(str(toolkit_folder)) #turn of in operation mode


import vlinder_toolkit
#%%


class App_layout:
    def set_root(self): #pass app
        """ Set the root of the current document """
        
        
        bokeh.io.curdoc().clear()
        bokeh.io.curdoc().add_root(self.IO_page) 
    
    
    def build_io_page(self, widgets):
        
        #Define rows and scaling
        datarow = row(widgets['datapath'], widgets['datatemplatepath'], widgets['make_data_template'])
        metadatarow = row(widgets['metadatapath'], widgets['metadatatemplatepath'], widgets['make_metadata_template'])
        
        widget_colmn = column([widgets['text'], datarow, metadatarow], sizing_mode='inherit')
        
        # init datatable

        source = init_datatable_with_default(settingslist=[gui_settings.dt_names,
                                                           gui_settings.meta_names,
                                                           gui_settings.obs_names])
        columns = [TableColumn(field=col, title=col) for col in gui_settings.data_table_column_order]
   
        self.data_table = DataTable(source=source, columns=columns, width=1400, height=400)
        
        table_column = column(self.data_table)
        
        
        
        
        #Define glyphs
        
        
        #create layout + SIZING MODE ON PARENT LEVEL (otherwise too complex)
        layout_io_page = layout([
            row([widget_colmn, table_column])
            ], 
            sizing_mode='stretch_width')
        
        return layout_io_page
    def cleanup_io_page(self):
    
        self.IO_page.children = self.IO_page.children[:3]
    
    
    def add_extra_row_io_page(self, rows, page='IO_page', column=None, sizemode='inherit', blank_avobe=True):
        """
        Add extra widgets (structured in rows or columns) to the IO page as a row at the bottom. 
        

        Parameters
        ----------
        rows : bokeh.row
            row containing widgets or glyphs.
        sizemode : TYPE, optional
            DESCRIPTION. The default is 'inherit'.
        blank_avobe : TYPE, optional
            DESCRIPTION. The default is True.

        Returns
        -------
        None.

        """
        blk_space=20
        new_column = []
        if isinstance(rows, list):
            i = 0
            for row_in_rows in rows:
                row_in_rows.sizing_mode = sizemode
                if (i == 0) & (blank_avobe):
                    row_in_rows.margin=(blk_space,0,0,0)
                new_column.append(row_in_rows)
                # children = getattr(self, page).children
                # children.append(row_in_rows)
                # setattr(getattr(self, page), 'children', children)
                
                
                i += 1
        else:
            rows.sizing_mode=sizemode
            if (blank_avobe):
                rows.margin=(blk_space,0,0,0)
            new_column.append(rows)
            # children = getattr(self, page).children
            # children.append(rows)
            # setattr(getattr(self, page), 'children', children)
          
            
        # # get columns current children
        # if not isinstance(column, type(None)):
        #     children = self.IO_page.children[0].children[column].children
        # else: 
        children = getattr(self, page).children

        children.extend(new_column)
        setattr(getattr(self, page), 'children', children)
    
    
    
    def __init__(self, app):
        
        self.IO_page = self.build_io_page(app.widg['io']) #build page IO
        
        
        # self.set_root()

    

class App:
    def __init__(self):
        
        self.widg={} #nested dict
        
        
        #init data objects to update
        self.mapperdf = pd.DataFrame()
        
        
        
        #init widgetsn
        
        self.define_widgets_io()
        self.layout = App_layout(self)
        

    def mapper_to_df(self, attr):
        print('mapper triggered')
        
        
        df = pd.DataFrame()
        
        categories = [gui_settings.dt_names, gui_settings.meta_names, gui_settings.obs_names]
        
        
        
        # dt spinners to dataframe
        for category in categories:
            mapper = {key : {colname: self.widg['io']['selectors'][key][colname].value for colname in val.keys()} for key, val in category.items()}
         
            for key in mapper:
                mapper[key]['template_column_name'] = self.widg['io']['selectors'][key]['template_column_name'].value
    
            # convert to pandas dataframe
            cat_df = pd.DataFrame(mapper).transpose()
            df = pd.concat([df, cat_df])
       
        df = df.reset_index()
        df = df.rename(columns={'index': 'toolkit_name'})
        print('df: ', df)
        # update table
        self.layout.data_table.source.data = df
        print(self.layout.data_table.source.data)
        print('updated ???')        
        
        
    def create_data_table(self, dataframe):
        source = ColumnDataSource(data)

        columns = [
                TableColumn(field="dates", title="Date", formatter=DateFormatter()),
                TableColumn(field="downloads", title="Downloads"),
            ]
        data_table = DataTable(source=source, columns=columns, width=400, height=280)
        
        show(data_table)
        
    def trig_make_data_template(self, attr):
        
        self.layout.cleanup_io_page()
        
        data_file_path = self.widg['io']['datapath'].value
        # check if path is empty
        if data_file_path=='':
            self.widg['io']['text'].text = text_increase('Path to datafile is empty!', 2)
            return
        # add csv extension if needed
        if data_file_path[-4:]!='.csv':
            data_file_path = data_file_path+'.csv'
        
        #check if file exists
        if not os.path.isfile(data_file_path):
            self.widg['io']['text'].text = text_increase(f'{data_file_path} is not found!', 2)
            return
        
        #read file, and get columnnames
        df = read_csv(data_file_path)
        data_columns  = list(df.columns)
        
        
        #Create grid and store values in app class
        #TODO: een lijst meegeven in plaats van argumenten, zo kan ik dezlfde functie gebruiken voor de metadata
        widget_grid, self.widg['io']['selectors'] = make_data_template_spinners(data_file_columns = data_columns,
                                                                                        tlk_meta_names = gui_settings.meta_names,
                                                                                        tlk_dt_names = gui_settings.dt_names,
                                                                                        tlk_obs_names = gui_settings.obs_names)
        
        # self.layout.add_extra_row_io_page(widget_grid)
        
        
        #add button to bottom
        test = Button(button_type='success', label="Create template") #define button
        test.on_click(self.mapper_to_df) #make callback
        self.layout.add_extra_row_io_page(row(test), column=0) #add to layout
        
    def trig_make_metadata_template(self, attr):
        
        
      
        self.layout.cleanup_io_page()
        
        data_file_path = self.widg['io']['metadatapath'].value
        # check if path is empty
        if data_file_path=='':
            self.widg['io']['text'].text = text_increase('Path to Metadatafile is empty!', 2)
            return
        # add csv extension if needed
        if data_file_path[-4:]!='.csv':
            data_file_path = data_file_path+'.csv'
        
        #check if file exists
        if not os.path.isfile(data_file_path):
            self.widg['io']['text'].text = text_increase(f'{data_file_path} is not found!', 2)
            return
        
        #read file, and get columnnames
        df = read_csv(data_file_path)
        data_columns  = list(df.columns)
        
        
        #Create grid and store values in app class
        #TODO: een lijst meegeven in plaats van argumenten, zo kan ik dezlfde functie gebruiken voor de metadata
        widget_grid, self.widg['io']['selectors'] = make_data_template_spinners(data_file_columns = data_columns,
                                                                                        tlk_meta_names = gui_settings.meta_names,
                                                                                        tlk_dt_names = gui_settings.dt_names,
                                                                                        tlk_obs_names = None)
        
        self.layout.add_extra_row_io_page(widget_grid)
        
        
        #add button to bottom
        print("add button")
        test = Button(button_type='success', label="Create template") #define button
        test.on_click(self.mapper_to_df) #make callback
        self.layout.add_extra_row_io_page(row(test)) #add to layout
    
    
    
        
        
    def define_widgets_io(self):
        self.widg['io'] = {} #init io dict
        
        #Create widgets
        self.widg['io']['text'] = Div(text=text_increase('init text blabla',2),
                                      sizing_mode='stretch_width')
        self.widg['io']['datapath'] = TextInput(title='Path to data file',
                                                value='/home/thoverga/Documents/VLINDER_github/vlinder_toolkit/tests/test_data/vlinderdata_small',
                                                sizing_mode='stretch_width',
                                                suffix='.csv')
        self.widg['io']['metadatapath'] = TextInput(title='Path to metadata file (if needed)',
                                                    value='/home/thoverga/Documents/VLINDER_github/vlinder_toolkit/vlinder_toolkit/data_templates/template_defaults/vlinder_metadata_template.csv',
                                                sizing_mode='stretch_width',
                                                suffix='.csv')
        self.widg['io']['datatemplatepath'] = TextInput(title='Path to template for data file',
                                                sizing_mode='stretch_width',
                                                suffix='.csv')
        self.widg['io']['metadatatemplatepath'] = TextInput(title='Path to template for metadata file (if needed)',
                                                sizing_mode='stretch_width',
                                                suffix='.csv')
        self.widg['io']['make_data_template'] = Button(button_type='warning', label="Make data template")
        self.widg['io']['make_metadata_template'] = Button(button_type='warning', label="Make metadata template")
        
        
        
        #Setup callbacks
        self.widg['io']['make_data_template'].on_click(self.trig_make_data_template)
        self.widg['io']['make_metadata_template'].on_click(self.trig_make_metadata_template)
        

    

        
dummy = App()


        