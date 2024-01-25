#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class defenition for regular observation types. The default observationtypes
are define here aswell.
"""

import sys
import logging
from collections.abc import Iterable

import numpy as np

logger = logging.getLogger(__name__)


# =============================================================================
# Standard toolkit units for each observation type
# =============================================================================

tlk_std_units = {
    "temp": "Celsius",
    "radiation_temp": "Celsius",
    "humidity": "%",
    "precip": "mm/m²",
    "precip_sum": "mm/m² from midnight",
    "wind_speed": "m/s",
    "wind_gust": "m/s",
    "wind_direction": "° from north (CW)",
    "pressure": "pa",
    "pressure_at_sea_level": "pa",
}


# =============================================================================
# Aliases for units
# =============================================================================

temp_aliases = {
    "Celsius": [
        "celsius",
        "°C",
        "°c",
        "celcius",
        "Celcius",
    ],  # for the dyselectic developper..
    "Kelvin": ["K", "kelvin"],
    "Farenheit": ["farenheit"],
}
pressure_aliases = {
    "pa": ["Pascal", "pascal", "Pa"],
    "hpa": ["hecto pascal", "hPa"],
    "psi": ["Psi"],
    "bar": ["Bar"],
}

precip_aliases = {"mm/m²": ["mm", "liter", "liters", "l/m²", "milimeter"]}

wind_aliases = {
    "m/s": ["meters/second", "m/sec"],
    "km/h": ["kilometers/hour", "kph"],
    "mph": ["miles/hour"],
}
direction_aliases = {"° from north (CW)": ["°", "degrees"]}


# conversion between standard-NAMES and aliases
all_units_aliases = {
    "temp": temp_aliases,
    "radiation_temp": temp_aliases,
    "humidity": {"%": ["percent", "percentage"]},
    "pressure": pressure_aliases,
    "pressure_at_sea_level": pressure_aliases,
    "precip": precip_aliases,
    "precip_sum": precip_aliases,
    "wind_speed": wind_aliases,
    "wind_gust": wind_aliases,
    "wind_direction": direction_aliases,
}

# =============================================================================
# Unit conversion expressions
# =============================================================================

all_conversion_table = {
    "temp": {
        "Kelvin": ["x - 273.15"],  # result is in tlk_std_units
        "Farenheit": ["x-32.0", "x/1.8"],
    },  # -->execute from left to write  = (x-32)/1.8
    "radiation_temp": {
        "Kelvin": ["x - 273.15"],  # result is in tlk_std_units
        "Farenheit": ["x-32.0", "x/1.8"],
    },
    "humidity": {},
    "pressure": {"hpa": ["x * 100"], "psi": ["x * 6894.7573"], "bar": ["x * 100000."]},
    "pressure_at_sea_level": {
        "hpa": ["x * 100"],
        "psi": ["x * 6894.7573"],
        "bar": ["x * 100000."],
    },
    "precip": {},
    "precip_sum": {},
    "wind_speed": {"km/h": ["x / 3.6"], "mph": ["x * 0.44704"]},
    "wind_gust": {"km/h": ["x / 3.6"], "mph": ["x * 0.44704"]},
    "wind_direction": {},
}

# =============================================================================
# Observation type class
# =============================================================================


class Obstype:
    """Object with all info and methods for a specific observation type."""

    def __init__(
        self, obsname, std_unit, description=None, unit_aliases={}, unit_conversions={}
    ):
        """Initiate an observation type.

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

        Returns
        -------
        None.

        """
        self.name = str(obsname)  # Standard name for the observation type
        self.std_unit = str(std_unit)  # standard unit fot the observation type
        self.description = str(description)

        # Conversion info and mappers
        self.units_aliases = unit_aliases
        self.conv_table = unit_conversions

        # Original column name and units in the data
        self.original_name = None  # Updated on IO
        self.original_unit = None  # updated on IO

        self._check_attributes()

    def __repr__(self):
        """Instance representation."""
        return f"Obstype instance of {self.name}"

    def __str__(self):
        """Text representation."""
        return f"Obstype instance of {self.name}"

    # -----  Setters -------

    def set_description(self, desc):
        """Set the description of the observation type."""
        self.description = str(desc)

    def set_original_name(self, columnname):
        """Set the original name of the observation type."""
        self.original_name = str(columnname)

    def set_original_unit(self, original_unit):
        """Set the original unit of the observation type."""
        self.original_unit = str(original_unit)

    # ------ Getters --------

    def get_info(self):
        """Print out detailed information of the observation type.

        Returns
        -------
        None.

        """
        info_str = f"{self.name} observation with: \n \
    * standard unit: {self.std_unit} \n \
    * data column as {self.original_name} in {self.original_unit} \n \
    * known units and aliases: {self.units_aliases} \n \
    * description: {self.description} \n \
    * conversions to known units: {self.conv_table} \n\n \
    * originates from data column: {self.original_name} with {self.original_unit} as native unit."
        print(info_str)

    def get_orig_name(self):
        """Return the original name of the observation type."""
        return self.original_name

    def get_description(self):
        """Return the descrition of the observation type."""
        if self.description == str(None):
            return "No description available"
        else:
            return str(self.description)

    def get_all_units(self):
        """Return a list with all the known unit (in standard naming)."""
        units = list(self.units_aliases.keys())
        units.append(self.get_standard_unit())
        return list(set(units))

    def get_standard_unit(self):
        """Return the standard unit of the observation type."""
        return self.std_unit

    def get_plot_y_label(self, mapname=None):
        """Return a string to represent the vertical axes of a plot."""
        return f"{self.name} ({self.std_unit})"

    def add_unit(self, unit_name, conversion=["x"]):
        """Add a new unit to an observation type.

        Parameters
        ----------
        unit_name : str
            The name of the new unit.
        conversion : list, optional
            The conversion description to the standard unit. The default is
            ["x"].

        Returns
        -------
        None.

        """
        # check if unit name is already known
        known = self.test_if_unit_is_known(unit_name)
        if known:
            return

        # convert expression to list if it is a string
        if isinstance(conversion, str):
            conversion = [conversion]

        # add converstion to the table
        self.conv_table[str(unit_name)] = conversion

        # add to alias table (without aliasses)
        self.units_aliases[unit_name] = []

        logger.info(
            f"{unit_name} is added as a {self.name} unit with coversion: {conversion} to {self.std_unit}"
        )

    def convert_to_standard_units(self, input_data, input_unit):
        """Convert data from a knonw unit to the standard unit.

        The data can be a collection of numeric values or a single numeric
        value.

        Parameters
        ----------
        input_data : (collection of) numeric
            The data to convert to the standard unit.
        input_unit : str
            The known unit the inputdata is in.

        Returns
        -------
        data  numeric/numpy.array
            The data in standard units.

        """
        # check if input unit is known
        known = self.test_if_unit_is_known(input_unit)

        # error when unit is not know
        if not known:
            sys.exit(
                f"{input_unit} is an unknown unit for {self.name}. No coversion possible!"
            )

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

    # ------------- Helpers ----------------------------------

    def _check_attributes(self):
        """Add units from the conv_table to the aliases if needed."""
        add_to_aliases = {}
        all_std_unit_names = []
        all_aliases = []
        for std_unit, alias_units in self.units_aliases.items():
            all_std_unit_names.append(std_unit)
            all_aliases.extend(alias_units)

        # add empty alias for all obstype present in conv table if no aliases are given
        for unit in self.conv_table.keys():
            if unit not in all_std_unit_names:
                if unit not in all_aliases:
                    add_to_aliases[unit] = []
        # add std unit to aliases if it is not already present
        if self.get_standard_unit() not in all_std_unit_names:
            add_to_aliases[self.get_standard_unit()] = []

        self.units_aliases.update(add_to_aliases)

    def _get_std_unit_name(self, unit_name):
        """Get standard name for a unit name by scanning trough the aliases."""
        for std_unit_name, aliases in self.units_aliases.items():
            if unit_name == std_unit_name:
                return unit_name
            if unit_name in aliases:
                return std_unit_name
        sys.exit(f"No standard unit name is found for {unit_name} for {self.name}")

    def test_if_unit_is_known(self, unit_name):
        """Test is the unit is known.

        Parameters
        ----------
        unit_name : str
            The unit name to test.

        Returns
        -------
        bool
            True if knonw, False else.

        """
        if unit_name == self.std_unit:
            return True
        for std_unit_name, aliases in self.units_aliases.items():
            if unit_name == std_unit_name:
                return True
            if unit_name in aliases:
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

temperature = Obstype(
    obsname="temp",
    std_unit=tlk_std_units["temp"],
    description="2m - temperature",
    unit_aliases=all_units_aliases["temp"],
    unit_conversions=all_conversion_table["temp"],
)

humidity = Obstype(
    obsname="humidity",
    std_unit=tlk_std_units["humidity"],
    description="2m - relative humidity",
    unit_aliases=all_units_aliases["humidity"],
    unit_conversions=all_conversion_table["humidity"],
)

radiation_temp = Obstype(
    obsname="radiation_temp",
    std_unit=tlk_std_units["radiation_temp"],
    description="2m - Black globe",
    unit_aliases=all_units_aliases["radiation_temp"],
    unit_conversions=all_conversion_table["radiation_temp"],
)

pressure = Obstype(
    obsname="pressure",
    std_unit=tlk_std_units["pressure"],
    description="atmospheric pressure (at station)",
    unit_aliases=all_units_aliases["pressure"],
    unit_conversions=all_conversion_table["pressure"],
)

pressure_at_sea_level = Obstype(
    obsname="pressure_at_sea_level",
    std_unit=tlk_std_units["pressure_at_sea_level"],
    description="atmospheric pressure (at sea level)",
    unit_aliases=all_units_aliases["pressure_at_sea_level"],
    unit_conversions=all_conversion_table["pressure_at_sea_level"],
)

precip = Obstype(
    obsname="precip",
    std_unit=tlk_std_units["precip"],
    description="precipitation intensity",
    unit_aliases=all_units_aliases["precip"],
    unit_conversions=all_conversion_table["precip"],
)

precip_sum = Obstype(
    obsname="precip_sum",
    std_unit=tlk_std_units["precip"],
    description="Cummulated precipitation",
    unit_aliases=all_units_aliases["precip_sum"],
    unit_conversions=all_conversion_table["precip_sum"],
)
wind = Obstype(
    obsname="wind_speed",
    std_unit=tlk_std_units["wind_speed"],
    description="wind speed",
    unit_aliases=all_units_aliases["wind_speed"],
    unit_conversions=all_conversion_table["wind_speed"],
)

windgust = Obstype(
    obsname="wind_gust",
    std_unit=tlk_std_units["wind_gust"],
    description="wind gust",
    unit_aliases=all_units_aliases["wind_gust"],
    unit_conversions=all_conversion_table["wind_gust"],
)

wind_direction = Obstype(
    obsname="wind_direction",
    std_unit=tlk_std_units["wind_direction"],
    description="wind direction",
    unit_aliases=all_units_aliases["wind_direction"],
    unit_conversions=all_conversion_table["wind_direction"],
)

# The order of the dictionary is also the order on how columns in dataset are presetnted
tlk_obstypes = {
    "temp": temperature,
    "humidity": humidity,
    "radiation_temp": radiation_temp,
    "pressure": pressure,
    "pressure_at_sea_level": pressure_at_sea_level,
    "precip": precip,
    "precip_sum": precip_sum,
    "wind_speed": wind,
    "wind_gust": windgust,
    "wind_direction": wind_direction,
}
