#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 13 13:43:15 2023

@author: thoverga
"""
import sys
import math
import numpy as np
import logging
from metobs_toolkit.obstypes import Obstype

from metobs_toolkit.obstypes import (temperature,
                                     pressure,
                                     wind,
                                     direction_aliases)

logger = logging.getLogger(__name__)

# =============================================================================
# Standard modeldata equivalents
# =============================================================================
tlk_std_modeldata_obstypes = {
    'temp': {"ERA5_hourly": {'name': 'temperature_2m', 'units': 'Kelvin',
                             'band_desc': "Temperature of air at 2m above the surface of land, sea or in-land waters. 2m temperature is calculated by interpolating between the lowest model level and the Earth's surface, taking account of the atmospheric conditions."}},

    'pressure': {'ERA5_hourly': {'name': 'surface_pressure', 'units': 'pa',
                                 'band_desc': "Pressure (force per unit area) of the atmosphere on the surface of land, sea and in-land water. It is a measure of the weight of all the air in a column vertically above the area of the Earth's surface represented at a fixed point. Surface pressure is often used in combination with temperature to calculate air density. The strong variation of pressure with altitude makes it difficult to see the low and high pressure systems over mountainous areas, so mean sea level pressure, rather than surface pressure, is normally used for this purpose. The units of this variable are Pascals (Pa). Surface pressure is often measured in hPa and sometimes is presented in the old units of millibars, mb (1 hPa = 1 mb = 100 Pa)."}},
    'u_wind': {'ERA5_hourly': {'name': 'u_component_of_wind_10m', 'units': 'm/s',
                                 'band_desc': "Eastward component of the 10m wind. It is the horizontal speed of air moving towards the east, at a height of ten meters above the surface of the Earth, in meters per second. Care should be taken when comparing this variable with observations, because wind observations vary on small space and time scales and are affected by the local terrain, vegetation and buildings that are represented only on average in the ECMWF Integrated Forecasting System. This variable can be combined with the V component of 10m wind to give the speed and direction of the horizontal 10m wind."}},
    'v_wind': {'ERA5_hourly': {'name': 'v_component_of_wind_10m', 'units': 'm/s',
                                  'band_desc': "Northward component of the 10m wind. It is the horizontal speed of air moving towards the north, at a height of ten meters above the surface of the Earth, in meters per second. Care should be taken when comparing this variable with observations, because wind observations vary on small space and time scales and are affected by the local terrain, vegetation and buildings that are represented only on average in the ECMWF Integrated Forecasting System. This variable can be combined with the U component of 10m wind to give the speed and direction of the horizontal 10m wind."}},

}


class ModelObstype(Obstype):
    """Extension of the Obstype class specific for the obstypes of Modeldata."""

    def __init__(self, obstype,
                 model_equivalent_dict={}):
        """Initiate an Modelobservation type.

        A ModelObstype has the same properties as an Obstype but with some
        extra attributes and methods.

        Parameters
        ----------
        obsname : str
            The name of the new observation type (i.g. 'sensible_heat_flux').
        std_unit : str
            The standard unit for the observation type (i.g. 'J/m²')
        obstype_description : str, ptional
            A more detailed description of the obstype (i.g. '2m SE inside
            canopy'). The default is None.
        unit_aliases : dict, optional
            A dictionary containing unit alias names. Keys represent a unit and
            values are lists with aliases for the units at the keys. The default is {}.
        unit_conversions : dict, optional
            A dictionary containing the conversion information to map to the
            standard units. Here an example of for temperatures (with Celcius
            as standard unit):

                {'Kelvin': ["x - 273.15"], #result is in tlk_std_units
                'Farenheit' : ["x-32.0", "x/1.8"]}, # -->execute from left to write  = (x-32)/1.8

                The default is {}.

        model_equiv_dict : dict
            A dictionary with information of how the observation type is found in
            modeldata. A example for pressure is:

                {'ERA5_hourly': {'name': 'surface_pressure', 'units': 'pa',
                                             'band_desc': "Pressure (force per ....

        Returns
        -------
        None.

        """
        super().__init__(obsname=obstype.name,
                         std_unit=obstype.std_unit,
                         description=obstype.description,
                         unit_aliases=obstype.units_aliases,
                         unit_conversions=obstype.conv_table)

        self.modl_equi_dict = model_equivalent_dict
        self._is_valid()

    def get_info(self):
        """Print out detailed information of the observation type.

        Returns
        -------
        None.

        """
        databands = {key: item['name'] for key, item in self.modl_equi_dict.items()}
        info_str = f"{self.name} observation with: \n \
    * Known datasetsbands: {databands} \n \
    * standard unit: {self.std_unit} \n \
    * description: {self.description} \n \
    * conversions to known units: {self.conv_table} \n"
        print(info_str)

    def get_mapped_datasets(self):
        """Return all gee datasets with a representing band for this obstype."""
        return list(self.modl_equi_dict.keys())

    def get_bandname(self, mapname):
        """Return the representing bandname of the obstype from a given gee dataset."""
        return str(self.modl_equi_dict[mapname]['name'])

    def get_bandname_mapper(self, mapname):
        """Return the representing bandname with tlk standard name as a dict."""
        return {str(self.modl_equi_dict[mapname]['name']): self.name}

    def get_plot_y_label(self, mapname):
        """Return a string to represent the vertical axes of a plot."""
        return f'{self.name} ({self.std_unit}) \ {mapname}: {self.modl_equi_dict[mapname]["name"]}'

    def get_modelunit(self, mapname):
        """Return the units of the representing bandname of the obstype from a given gee dataset."""
        return str(self.modl_equi_dict[mapname]['units'])

    def has_mapped_band(self, mapname):
        """Test is a gee dataset has a representing band."""
        try:
            self.get_bandname(mapname)
            return True
        except KeyError:
            return False

    def add_new_band(self, mapname, bandname, bandunit, band_desc=None):
        """Add a new representing dataset/bandname to the obstype.

        Parameters
        ----------
        mapname : str
            name of the known gee dataset.
        bandname : str
            the name of the representing band.
        bandunit : str
            the unit of the representing band.
        band_desc : str, optional
            A detailed description of the band.

        Returns
        -------
        None.

        """
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
        """Test if all attributes are valid among each other."""
        for datasetname in self.modl_equi_dict.keys():
            # Check if unit is available
            if 'units' not in self.modl_equi_dict[datasetname].keys():
                sys.exit(f'No units information is provided for {self.name} for modeldata: {datasetname}')
            # check if the unit is known
            if not self.test_if_unit_is_known(unit_name=self.modl_equi_dict[datasetname]['units']):
                sys.exit(f'Cannot create {self.name} ModelObstype because {self.modl_equi_dict[datasetname]["units"]} is a unknown unit.')


class ModelObstype_Vectorfield(Obstype):
    def __init__(self, obstype,
                 u_comp_model_equivalent_dict={},
                 v_comp_model_equivalent_dict={}):

        super().__init__(obsname=obstype.name,
                         std_unit=obstype.std_unit,
                         description=obstype.description,
                         unit_aliases=obstype.units_aliases,
                         unit_conversions=obstype.conv_table)


        if set(u_comp_model_equivalent_dict.keys()) != set(v_comp_model_equivalent_dict.keys()):
            sys.exit(f'The mapped gee dataset are not equal for the vector components of {obstype.name}.')

        mod_comp_dict = {}
        for geedataset in u_comp_model_equivalent_dict.keys():
            mod_comp_dict[geedataset] = {'u_comp': u_comp_model_equivalent_dict[geedataset],
                                         'v_comp': v_comp_model_equivalent_dict[geedataset]}


        self.modl_comp_dict = mod_comp_dict
        self._is_valid()

    def get_info(self):
        """Print out detailed information of the observation type.

        Returns
        -------
        None.

        """
        u_databands = {key: item['u_comp']['name'] for key, item in self.modl_comp_dict.items()}
        v_databands = {key: item['v_comp']['name'] for key, item in self.modl_comp_dict.items()}
        info_str = f"{self.name} observation with: \n \
    * Known Vector-East-component datasetsbands: {u_databands} \n \
    * Known Vector-North-component datasetsbands: {v_databands} \n \
    * standard unit: {self.std_unit} \n \
    * description: {self.description} \n \
    * conversions to known units: {self.conv_table} \n"
        print(info_str)

    def get_mapped_datasets(self):
        """Return all gee datasets with a representing band for this obstype."""
        return list(self.modl_comp_dict.keys())

    # def get_bandname(self, mapname):
    #     """Return the representing bandname of the obstype from a given gee dataset."""
    #     return str(self.modl_equi_dict[mapname]['name'])

    def get_bandname_mapper(self, mapname):
        """Return the representing bandname with tlk standard name as a dict."""
        mapper = {str(self.modl_comp_dict[mapname]['u_comp']['name']): f'u_comp_{self.name}',
                  str(self.modl_comp_dict[mapname]['v_comp']['name']): f'v_comp_{self.name}'}

        return mapper

    def get_modelunit(self, mapname):
        """Return the units of the representing bandname of the obstype from a given gee dataset."""
        # u and v comp must have the same units, this is tested in the _is_valid()
        return str(self.modl_comp_dict[mapname]['u_comp']['units'])

    def has_mapped_band(self, mapname):
        """Test is a gee dataset has a representing band."""
        if mapname in self.modl_comp_dict.keys():
            return True
        else:
            return False

    def get_plot_y_label(self, mapname):
        """Return a string to represent the vertical axes of a plot."""
        return f'{self.name} ({self.std_unit}) \ {mapname}: {self.modl_equi_dict[mapname]["u_comp"]["name"]} and {self.modl_equi_dict[mapname]["v_comp"]["name"]}'

    def get_u_column(self):
        return f'u_comp_{self.name}'
    def get_v_column(self):
        return f'v_comp_{self.name}'

    def compute_amplitude(self, df):
        # Compute the data
        data =((df[self.get_u_column()].pow(2)) + (df[self.get_v_column()].pow(2))).pow(1./2)
        # Create a new obstype for the amplitude
        amplitude_obstype = Obstype(obsname = f'{self.name}_amplitude',
                                    std_unit = self.std_unit,
                                    description=f'2D-vector amplitde of {self.name} components.',
                                    unit_aliases = self.units_aliases,
                                    unit_conversions = self.conv_table)
        # convert to model obstype
        mod_equi = {}
        for key, val in self.modl_comp_dict.items():
            mod_equi[key] = val['u_comp']
            mod_equi[key]['name'] = f"{val['u_comp']['name']} and {val['v_comp']['name']}"

        amplitude_obstype = ModelObstype(amplitude_obstype,
                                         model_equivalent_dict=mod_equi
                                         )
        # amplitude_obstype.plotlabel=f"Amplitude of 2D-{self.name} field ({self.std_unit})"
        return data, amplitude_obstype

    def compute_angle(self, df):

        def unit_vector(vector):
            """ Returns the unit vector of the vector.  """
            return vector / np.linalg.norm(vector)

        def angle_between(u_comp, v_comp):
            """ Returns the angle in radians between vectors 'v1' and 'v2'::

                    >>> angle_between((1, 0, 0), (0, 1, 0))
                    1.5707963267948966
                    >>> angle_between((1, 0, 0), (1, 0, 0))
                    0.0
                    >>> angle_between((1, 0, 0), (-1, 0, 0))
                    3.141592653589793
            """

            v2 = (u_comp, v_comp)
            v1_u = unit_vector((0,1)) #North unit arrow
            v2_u = unit_vector(v2)

            angle_rad = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
            angle_degrees = angle_rad * ((180.0/math.pi))
            # return angle_degrees
            # fix the quadrants
            if (v2[0] >= 0) & (v2[1] >= 0):
                # N-E quadrant
                return angle_degrees
            if (v2[0] >= 0) & (v2[1] < 0):
                # S-E quadrant
                return angle_degrees
            if (v2[0] < 0) & (v2[1] < 0):
                # S-W quadrant
                return 180.0 + (180. - angle_degrees)
            if (v2[0] < 0) & (v2[1] >= 0):
                # N-W quadrant
                return 360.0 - angle_degrees

        u_column = self.get_u_column()
        v_column = self.get_v_column()

        data = df.apply(lambda x: angle_between(x[u_column], x[v_column]), axis=1)
        # Create a new obstype for the amplitude
        direction_obstype = Obstype(obsname = f'{self.name}_direction',
                                    std_unit = '° from north (CW)',
                                    description=f'Direction of 2D-vector of {self.name} components.',
                                    unit_aliases = direction_aliases,
                                    unit_conversions = {})
        # convert to model obstype
        mod_equi = {}
        for key, val in self.modl_comp_dict.items():
            mod_equi[key] = val['u_comp']
            mod_equi[key]['name'] = f"{val['u_comp']['name']} and {val['v_comp']['name']}"
            mod_equi[key]['units'] = '° from north (CW)'

        direction_obstype = ModelObstype(direction_obstype,
                                         model_equivalent_dict=mod_equi
                                         )
        # direction_obstype.plotlabel=f"Direction of 2D-{self.name} field ({self.std_unit})"
        return data, direction_obstype






    def add_new_band(self, mapname, bandname_u_comp, bandname_v_comp, bandunit,
                     band_desc_u_comp=None, band_desc_v_comp=None):
        """Add a new representing dataset/bandname to the obstype.

        Parameters
        ----------
        mapname : str
            name of the known gee dataset.
        bandname_u_comp : str
            the name of the representing the Eastwards component band.
        bandname_v_comp : str
            the name of the representing the Northwards component band.
        bandunit : str
            the unit of the representing bands.
        band_desc_u_comp : str, optional
            A detailed description of the Eastwards component of the band.
        band_desc_v_comp : str, optional
            A detailed description of the Northwards component of the band.

        Returns
        -------
        None.

        """
        # test if banunit is valid
        if not self.test_if_unit_is_known(bandunit):
            sys.exit(f'{bandunit} is an unknown unit for the {self.name} obstype.')


        if mapname in self.modl_comp_dict.keys():
            # check if band is already knonw
            logger.debug(f'Update {bandname} of (known) map: {mapname}')
        else:
            logger.debug(f'Add new map: {mapname} with band: {bandname}.')


        self.modl_comp_dict[mapname] = {}
        self.modl_comp_dict[mapname]['u_comp'] = {'name': str(bandname_u_comp),
                                                  'units': str(bandunit),
                                                  'band_desc': str(band_desc_u_comp)}
        self.modl_comp_dict[mapname]['v_comp'] = {'name': str(bandname_v_comp),
                                                  'units': str(bandunit),
                                                  'band_desc': str(band_desc_v_comp)}

    def _is_valid(self):
        """Test if all attributes are valid among each other."""
        for datasetname in self.modl_comp_dict.keys():
            for comp_str, comp in self.modl_comp_dict[datasetname].items():
                # Check if unit is available
                if 'units' not in comp.keys():
                    sys.exit(f'No units information is provided for {self.name} for {comp_str} modeldata_vectorfield: {datasetname}')
                # check if the unit is known
                if not self.test_if_unit_is_known(unit_name=comp['units']):
                    sys.exit(f'Cannot create {self.name} ModelObstype_Vectorfield because {comp["units"]} is a unknown unit in the {comp_str}.')

            # check if the units of the u and v comp are equal
            if len(set([comp['units'] for comp in self.modl_comp_dict[datasetname].values()])) > 1:
                sys.exit(f'The units of the u and v component for {self.name} in the {datasetname} dataset are not equal.')



# =============================================================================
# Define obstypes
# =============================================================================

temp_model = ModelObstype(temperature, model_equivalent_dict=tlk_std_modeldata_obstypes['temp'])
pressure_model = ModelObstype(pressure, model_equivalent_dict=tlk_std_modeldata_obstypes['pressure'])

# Special obstypes
wind_model = ModelObstype_Vectorfield(wind,
                                      u_comp_model_equivalent_dict=tlk_std_modeldata_obstypes['u_wind'],
                                      v_comp_model_equivalent_dict=tlk_std_modeldata_obstypes['v_wind'])


# =============================================================================
# Create obstype dict
# =============================================================================
model_obstypes = {'temp': temp_model,
                  'pressure': pressure_model,
                  'wind': wind_model,
                  }