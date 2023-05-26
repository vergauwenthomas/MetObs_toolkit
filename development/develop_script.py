#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:25:02 2022
@author: thoverga
"""

#%%

# import metobs_toolkit

import os
import sys
from pathlib import Path


lib_folder = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(lib_folder))

import metobs_toolkit


#%%



# Make an empty dataset
dataset = metobs_toolkit.Dataset()

# Add the demo data files to the dataset settings
dataset.update_settings(input_data_file = metobs_toolkit.demo_datafile,
                        input_metadata_file = metobs_toolkit.demo_metadatafile,
                        data_template_file = metobs_toolkit.demo_template,
                        metadata_template_file = metobs_toolkit.demo_template, # Contains also the metadata mapping
                        output_folder = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/development'
                        )

dataset.import_data_from_file()


dataset.coarsen_time_resolution(freq='30T')
dataset.apply_quality_control()


# deze functie zal outliers omzetten naar missing observations, en indien ze voldoen aan de defenitie ook gaps.
# Ik zet hier de gapsize laag (defenitie van een gaps is hier minimum 1 uur)
dataset.update_gaps_and_missing_from_outliers(n_gapsize=2)

# Gapfill settings kan je hier nog aanpassen als je wil:
dataset.update_gap_and_missing_fill_settings()





# Ofwel interpoleer je de gaps:

dataset.fill_gaps_linear(obstype='temp', overwrite=True)

# Ofwel gebruik je model data. Hiervoor moet je eerst een modelobject aanmaken:

# Er zijn twee manieren om dat te doen, afhankelijk van hoeveel data je nodig hebt.

era_model = dataset.get_modeldata(
                    modelname='ERA5_hourly',
                    stations=None, #None is allemaal
                    startdt=None, #neem de startdt van de observations
                    enddt=None) #neem de enddt van de observaties

# Als er weinig data nodg is dan wordt dit Modeldata object rechtstreeks opgevuld,
# Als er meer data nodig is dan zal er een boodschap verschijnen in de terminal
# dat de modeldata wordt weggeschreven in een file in de google drive  (folder: era5_timeseries/era5_data )

#Dan is het de bedoeling dat de gebruiker deze file teruggeeft aan het Modeldata object
era_model = metobs_toolkit.Modeldata('era5')
# donwload de csv file lokaal, of link naar het path in Colab
era_model.set_model_from_csv(csvpath='/home/thoverga/Downloads/era_5_data.csv')



# En dan geef je de modeldata mee aan de gapfill functie:

dataset.fill_gaps_era5(modeldata=era_model,
                       method='debias',
                       obstype='temp',
                       overwrite=True)

# (vul ook de de missing observations op (enkel interpolation))
dataset.fill_missing_obs_linear(obstype='temp')


# Als je nu plot (met colorby ='label') dan zal je zien dat de gaten opgevuld zijn
dataset.make_plot(colorby='label')


# Er is ook de automatische gapfill functie die een methode selecteerd aan de hand
# van de gapsize zoals je hebt voorgesteld:

# dataset.fill_gaps_automatic(modeldata=era_model,
#                             obstype="temp",
#                             max_interpolate_duration_str='6H',
#                             overwrite=True)




# om een gap en de fill techniek in detail te bekijken, dan gebruik je de .get_info() methode:

for gap in dataset.gaps:
    gap.get_info() #voor alle gaps

# als je de gafill details wilt, dan worden die opgeslagen in het .gapfill_info attribute:


dataset.gaps[2].gapfill_info


# PS: Je kan ook 'dataset' overal vervangen door een Station. Of eerst alles opvullen in de hele dataset,
# en dan een station eruit halen om te analyseren:

sta = dataset.get_station('vlinder05')

for gap in sta.gaps:
    gap.get_info()





#%%

# sta = dataset.get_station('vlinder05')
# sta.make_plot(title = 'before qc')

#%%


# #%%
# dataset.apply_quality_control()
# sta.make_plot(colorby='label', title='after qc')

#%%
# dataset.update_gaps_and_missing_from_outliers()
# sta.make_plot(colorby='label', title='after updating gaps')


#%%
# sta.fill_gaps_linear()
# sta.make_plot(colorby='label', title='after gaps fill')

# print(sta.gaps[0].get_info())

