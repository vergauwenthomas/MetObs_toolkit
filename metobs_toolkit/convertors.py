#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 24 09:34:00 2023

@author: thoverga
"""

import sys
import numpy as np
from collections.abc import Iterable

# =============================================================================
# Unit defenitions and coversions
# =============================================================================
# Keys are the toolkit-units!! (not persee SI)
# values expresions are of the form x $ val, where x is the numeric value,
# $ an operator (+-*/) and val a concersion value
unit_convertors = {
    "Celsius": {"K": "x - 273.15"},
}


# =============================================================================
# Convert functions
# =============================================================================


def expression_calculator(equation, x):
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


def convert_to_toolkit_units(data, data_unit):
    """
    Convert the data to the toolkit perfered unit. Data can be a numeric value or an iterable.
    Data_unit is the unit of the input data.

    The converted data AND the corresponding toolkit unit is returned.

    Parameters
    ----------
    data : numeric, iterable
        numeric data to be converted.
    data_unit : String
        unit name of the data.

    Returns
    -------
    numeric, numpy.array
        The data in toolkit units.
    String
        Corresponding toolkit unit name.

    """
    # check if unit is already a toolkit unit
    if data_unit in unit_convertors.keys():
        return data, data_unit

    # scan the units to find conversion
    expr = {
        toolk_unit: other_unit[data_unit]
        for toolk_unit, other_unit in unit_convertors.items()
        if data_unit in other_unit.keys()
    }

    if len(expr) == 1:  # unique conversion found
        conv_data = expression_calculator(next(iter(expr.values())), data)
        return conv_data, next(iter(expr.keys()))

    elif len(expr) > 1:
        sys.exit(f" Multiple possible conversions found for {data_unit}")
    else:
        sys.exit(f"No conversion found for {data_unit}")
