#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  6 16:08:24 2023

@author: thoverga
"""

import sys
import logging
from collections.abc import Iterable

import numpy as np
logger = logging.getLogger(__name__)



# =============================================================================
# Hardcode obstypes, conversion tables and aliases
# =============================================================================

tlk_std_units = {
    "temp": 'Celsius',
    "radiation_temp": 'Celcius',
    "humidity": '%',
    "precip": 'mm/m² per hour',
    "precip_sum": 'mm/m² from midnight',
    "wind_speed": 'm/s',
    "wind_gust": 'm/s',
    "wind_direction": '° from north (CW)',
    "pressure": 'pa',
    "pressure_at_sea_level": 'pa'}




# conversion between standard-NAMES and aliases
all_units_aliases = {
    "temp": {
        'Celcius': ['celcius', '°C', '°c'],
        'Kelvin': ['K', 'kelvin'],
        'Farenheit': []
        },
    "radiation_temp": {
        'Celcius': ['celcius', '°C', '°c'],
        'Kelvin': ['K', 'kelvin'],
        'Farenheit': []
        },
    "humidity": {
        "%": ['percent', 'percentage']
        },
    'pressure': {
        'pa': ['Pascal', 'pascal', 'Pa'],
        'hpa': ['hecto pascal', 'hPa'],
        'psi': ['Psi'],
        'bar': ['Bar']
        }
    }



all_conversion_table = {
    'temp': {
        'Kelvin': ["x - 273.15"], #result is in tlk_std_units
        'Farenheit' : ["x-32.0", "x/1.8"]}, # -->execute from left to write  = (x-32)/1.8
    'radiation_temp': {
        'Kelvin': ["x - 273.15"], #result is in tlk_std_units
        'Farenheit' : ["x-32.0", "x/1.8"]},
    'humidity': {},
    'pressure': {
        'hpa': ["x * 100"],
        'psi': ['x * 6894.7573'],
        'bar': ['x * 100000.']
    }

}







# =============================================================================
# Observation type class
# =============================================================================

class Obstype:
    """Object with all info and methods for a specific observation type."""

    def __init__(self, obsname, std_unit, description=None):

        self.name = str(obsname)
        self.std_unit = str(std_unit)
        self.description = str(description)

        self.units_aliases = all_units_aliases[obsname]
        self.conv_table = all_conversion_table[obsname]

    def __repr__(self):
        return f"Obstype instance of {self.name}"

    def __str__(self):
        info_str = f"{self.name} observation with: \n \
    * standard unit: {self.std_unit} \n \
    * description: {self.description} \n \
    * conversions to known units: {self.conv_table} \n"
        return info_str

    def get_standard_unit(self):
        print(self.std_unit)
        return self.std_unit

    def add_unit(self, unit_name, conversion=["x"]):
        # check if unit name is already known
        known = self._test_if_unit_is_known(unit_name)
        if known:
            return

        # convert expression to list if it is a string
        if isinstance(conversion, str):
            conversion = [conversion]

        # add converstion to the table
        self.conv_table[str(unit_name)] = conversion

        # add to alias table (without aliasses)
        self.units_aliases[unit_name] = []

        logger.info(f'{unit_name} is added as a {self.name} unit with coversion: {conversion} to {self.std_unit}')


    def convert_to_standard_units(self, input_data, input_unit):
        # check if input unit is known
        known = self._test_if_unit_is_known(input_unit)

        # error when unit is not know
        if not known:
            sys.exit(f'{input_unit} is an unknonw unit for {self.name}. No coversion possible!')

        # Get conversion
        std_unit_name = self._get_std_unit_name(input_unit)
        if std_unit_name == self.std_unit:
            # No conversion needed because already the standard unit
            return input_data

        conv_expr_list = self.conv_table[std_unit_name]

        # covert data
        data = input_data
        for conv in conv_expr_list:
            data = expression_calculator(conv, data)

        return data


    def _get_std_unit_name(self, unit_name):
        for std_unit_name, aliases in self.units_aliases.items():
            if unit_name == std_unit_name:
                return unit_name
            if unit_name in aliases:
                return std_unit_name
        sys.exit(f"No standard unit name is found for {unit_name} for {self.name}")


    def _test_if_unit_is_known(self, unit_name):
        for std_unit_name, aliases in self.units_aliases.items():
            if unit_name == std_unit_name:
                logger.info(f'{unit_name} is a known unit for {self.name}.')
                return True
            if unit_name in aliases:
                logger.info(f'{unit_name} is a known (alias for {std_unit_name}) unit for {self.name}.')
                return True
        return False


def expression_calculator(equation, x):
    """Convert array by equation."""
    if isinstance(x, Iterable):
        x = np.array(x)

    if "+" in equation:
        y = equation.split("+")
        return x + float(y[1])
    elif "-" in equation:
        y = equation.split("-")
        return x - float(y[1])
    elif "/" in equation:
        y = equation.split("/")
        return x / float(y[1])
    elif "*" in equation:
        y = equation.split("*")
        return x * float(y[1])
    else:
        sys.exit(f"expression {equation}, can not be converted to mathematical.")



# =============================================================================
# Create observation types
# =============================================================================


temperature = Obstype(obsname='temp',
                      std_unit='Celcius',
                      description="2m - temperature")

humidity = Obstype(obsname='humidity',
                      std_unit='%',
                      description="2m - relative humidity")

radiation_temp = Obstype(obsname='radiation_temp',
                      std_unit='Celcius',
                      description="2m - Black globe")

pressure = Obstype(obsname='pressure',
                      std_unit='pa',
                      description="atmospheric pressure (at station)")



# =============================================================================
# Debugging
# =============================================================================



# temperature.add_unit('Koki', "x-100")

# testdata = [300, 320, 340]
# # testdata = 92

# testconv = temperature.convert_to_standard_units(input_data=testdata,
#                                                  input_unit = "Farenheit")

# print(testconv)


pres = [1, 1.2, 1.5]


pres_pa = pressure.convert_to_standard_units(pres, 'psi')

print(pres_pa)



