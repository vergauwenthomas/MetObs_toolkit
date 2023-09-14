#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 13:43:15 2023

@author: thoverga
"""
import sys
from metobs_toolkit.obstypes import Obstype

from metobs_toolkit.obstypes import (temperature,
                                     pressure)




#%%
# =============================================================================
# Standard modeldata equivalents
# =============================================================================
tlk_std_modeldata_obstypes = {
    'temp': { "ERA5_hourly": {'name': 'temperature_2m', 'units': 'Kelvin',
                             'band_desc': "Temperature of air at 2m above the surface of land, sea or in-land waters. 2m temperature is calculated by interpolating between the lowest model level and the Earth's surface, taking account of the atmospheric conditions."}},

    'pressure': {'ERA5_hourly': {'name': 'surface_pressure', 'units': 'pa',
                                 'band_desc': "Pressure (force per unit area) of the atmosphere on the surface of land, sea and in-land water. It is a measure of the weight of all the air in a column vertically above the area of the Earth's surface represented at a fixed point. Surface pressure is often used in combination with temperature to calculate air density. The strong variation of pressure with altitude makes it difficult to see the low and high pressure systems over mountainous areas, so mean sea level pressure, rather than surface pressure, is normally used for this purpose. The units of this variable are Pascals (Pa). Surface pressure is often measured in hPa and sometimes is presented in the old units of millibars, mb (1 hPa = 1 mb = 100 Pa)."}}
    }


#%%



class ModelObstype(Obstype):
    def __init__(self, obsname, std_unit, description=None, unit_aliases={},
               unit_conversions={},
               model_equivalent_dict={}):

        super().__init__(obsname, std_unit, description, unit_aliases,
                 unit_conversions)

        self.modl_equi_dict = model_equivalent_dict
        self._is_valid()

    def get_info(self):
        databands ={key: item['name'] for key, item in self.modl_equi_dict.items()}
        info_str = f"{self.name} observation with: \n \
    * Known datasetsbands: {databands} \n \
    * standard unit: {self.std_unit} \n \
    * description: {self.description} \n \
    * conversions to known units: {self.conv_table} \n"
        print(info_str)


    def get_mapped_datasets(self):
        return list(self.modl_equi_dict.keys())
    def get_bandname(self, mapname):
        return str(self.modl_equi_dict[mapname]['name'])
    def get_modelunit(self, mapname):
        return str(self.modl_equi_dict[mapname]['units'])

    def has_mapped_band(self, mapname):
        try:
            self.get_bandname(mapname)
            return True
        except KeyError:
            return False

    def add_new_band(self, mapname, bandname, bandunit, band_desc):
        # test if banunit is valid
        if not self.test_if_unit_is_known(bandunit):
            sys.exit(f'{bandunit} is an unknown unit for the {self.name} obstype.')

        if mapname in self.modl_equi_dict.keys():
            # check if band is already knonw
            logger.debug(f'Update {bandname} of (known) map: {mapname}')
        else:
            logger.debug(f'Add new map: {mapname} with band: {bandname}.')
        self.modl_equi_dict[mapname] = {'name': str(bandname),
                                        'units': str(bandunit),
                                         'band_desc': str(band_desc)}


    def _is_valid(self):
        for datasetname in self.modl_equi_dict.keys():
            # Check if unit is available
            if not 'units' in self.modl_equi_dict[datasetname].keys():
                sys.exit(f'No units information is provided for {self.name} for modeldata: {datasetname}')
            # check if the unit is known
            if not self.test_if_unit_is_known(unit_name = self.modl_equi_dict[datasetname]['units']):
                sys.exit(f'Cannot create {self.name} ModelObstype because {self.modl_equi_dict[datasetname]["units"]} is a unknown unit.')




# =============================================================================
# creator
# =============================================================================

def create_model_obstype(obstype, model_equiv_dict):
    """Construct a ModelObstype from a given Obstype.

    A Modelobstype is an inherited child of Obstype, with extra information
    of modeldata bands that represent the obstype.

    Parameters
    ----------
    obstype : metobs_toolkit.Obstype
        The observation type to add modeldata information to.
    model_equiv_dict : dict
        A dictionary with information of how the observation type is found in
        modeldata. A example for pressure is:

            {'ERA5_hourly': {'name': 'surface_pressure', 'units': 'pa',
                                         'band_desc': "Pressure (force per ....

    Returns
    -------
    mod_obstype : metobs_toolkit.ModelObstype
        The Obstype with extra modeldata attributes and methods.

    """
    mod_obstype = ModelObstype(obsname = obstype.name,
                               std_unit = obstype.std_unit,
                               description = None, #description of model is different than observations
                               unit_aliases = obstype.units_aliases,
                               unit_conversions = obstype.conv_table,
                               model_equivalent_dict = model_equiv_dict)
    return mod_obstype


# =============================================================================
# Define obstypes
# =============================================================================


temp_model = create_model_obstype(temperature, model_equiv_dict=tlk_std_modeldata_obstypes['temp'])
pressure_model = create_model_obstype(pressure, model_equiv_dict=tlk_std_modeldata_obstypes['pressure'])

# =============================================================================
# Create obstype dict
# =============================================================================
model_obstypes = {'temp': temp_model,
                   'pressure': pressure_model,
                  }
