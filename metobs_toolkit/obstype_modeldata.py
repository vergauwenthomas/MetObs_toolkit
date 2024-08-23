#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class defenition of model observationtypes. These are regular observationtypes
witht extra attributes and methods for interacting with the google earht engine.
"""
import sys
import copy
import math
import numpy as np
import logging
from metobs_toolkit.obstypes import Obstype, expression_calculator

from metobs_toolkit.obstypes import temperature, pressure, wind_speed, direction_aliases

logger = logging.getLogger(__name__)


class ModelObstype(Obstype):
    def __init__(self, obstype, model_unit, model_band):
        super().__init__(
            obsname=obstype.name,
            std_unit=obstype.std_unit,
            description=obstype.description,
            unit_aliases=obstype.units_aliases,
            unit_conversions=obstype.conv_table,
        )

        self.model_unit = str(model_unit)
        self.model_band = str(model_band)

        self._check_validity()

        # to make link with Obstype class
        self.set_original_name(model_band)
        self.set_original_unit(model_unit)

    def _check_validity(self):
        if self.model_unit not in self.get_all_units():
            raise MetobsModelObstypeHandlingError(
                f"{self.model_unit} is not a known unit of {self}."
            )

    # =============================================================================
    # Getters
    # =============================================================================

    def get_modelunit(self):
        return self.model_unit

    def get_modelband(self):
        return self.model_band

    def get_plot_y_label(self):
        """Return a string to represent the vertical axes of a plot."""
        return f"{self.name} ({self.std_unit})\n originates from {self.original_name}"


class ModelObstype_Vectorfield(Obstype):
    def __init__(
        self,
        obstype,
        model_unit,
        model_band_u,
        model_band_v,
        amplitude_obstype_name,
        direction_obstype_name,
    ):

        super().__init__(
            obsname=obstype.name,
            std_unit=obstype.std_unit,
            description=obstype.description,
            unit_aliases=obstype.units_aliases,
            unit_conversions=obstype.conv_table,
        )

        self.model_unit = str(model_unit)
        self.model_band_u = str(model_band_u)
        self.model_band_v = str(model_band_v)

        self._amp_obs_name = str(amplitude_obstype_name)
        self._dir_obs_name = str(direction_obstype_name)

        self._check_validity()

        # to make link with Obstype class
        self.set_original_name(f"{model_band_u} and {model_band_v}")
        self.set_original_unit(model_unit)

    def _check_validity(self):
        if self.model_unit not in self.get_all_units():
            raise MetobsModelObstypeHandlingError(
                f"{self.model_unit} is not a known unit of {self}."
            )

    def get_modelunit(self):
        return self.model_unit

    def get_modelband_u(self):
        return self.model_band_u

    def get_modelband_v(self):
        return self.model_band_v

    def get_plot_y_label(self):
        """Return a string to represent the vertical axes of a plot."""
        return f"{self.name} ({self.std_unit})\n originates from {self.original_name}"

    def compute_angle(self, df):
        """Compute vector direction of 2D vectorfield components.

        The direction column is added to the dataframe and a new ModelObstype,
        representing the angle is returned. The values represents the angles in
        degrees, from north in clock-wise rotation.

        Parameters
        ----------
        modelobs_vectorfield : ModelObstype_Vectorfield
            The vectorfield observation type to compute the vector directions for.
        df : pandas.DataFrame
            The dataframe with the vector components present as columns.

        Returns
        -------
        data : pandas.DataFrame
            The df with an extra column representing the directions.
        amplitude_obstype : ModelObstype
            The (scalar) Modelobstype representation of the angles.

        """

        def unit_vector(vector):
            """Returns the unit vector of the vector."""
            return vector / np.linalg.norm(vector)

        def angle_between(u_comp, v_comp):
            """Returns the angle in ° from North (CW) from 2D Vector components."""

            v2 = (u_comp, v_comp)
            v1_u = unit_vector((0, 1))  # North unit arrow
            v2_u = unit_vector(v2)

            angle_rad = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
            angle_degrees = angle_rad * ((180.0 / math.pi))
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
                return 180.0 + (180.0 - angle_degrees)
            if (v2[0] < 0) & (v2[1] >= 0):
                # N-W quadrant
                return 360.0 - angle_degrees

        u_column = self.get_modelband_u()
        v_column = self.get_modelband_v()

        data = df.apply(lambda x: angle_between(x[u_column], x[v_column]), axis=1)

        # Create a new obstype for the direction
        direction_obstype = Obstype(
            obsname=self._dir_obs_name,
            std_unit="° from north (CW)",
            description=f"Direction of 2D-vector of {self.name} components.",
            unit_aliases=direction_aliases,
            unit_conversions={},
        )
        # convert to model obstype
        direction_modelobstype = ModelObstype(
            obstype=direction_obstype,
            model_unit="° from north (CW)",  # indep of units
            model_band=None,
        )

        return data, direction_modelobstype

    def compute_amplitude(self, df):
        """Compute amplitude of 2D vectorfield components.

        The amplitude column is added to the dataframe and a new ModelObstype,
        representing the amplitude is returned. All attributes wrt the units are
        inherited from the ModelObstype_vectorfield.

        Parameters
        ----------
        modelobs_vectorfield : ModelObstype_Vectorfield
            The vectorfield observation type to compute the vector amplitudes for.
        df : pandas.DataFrame
            The dataframe with the vector components present as columns.

        Returns
        -------
        data : pandas.DataFrame
            The df with an extra column representing the amplitudes.
        amplitude_obstype : ModelObstype
            The (scalar) Modelobstype representation of the amplitudes.

        """
        # Compute the data
        data = (
            (df[self.get_modelband_u()].pow(2)) + (df[self.get_modelband_v()].pow(2))
        ).pow(1.0 / 2)

        # Create a new Obstype for the amplitude
        amplitude_obstype = Obstype(
            obsname=self._amp_obs_name,
            std_unit=self.std_unit,
            description=f"2D-vector amplitde of {self.name} components.",
            unit_aliases=self.units_aliases,
            unit_conversions=self.conv_table,
        )

        # convert to model obstype
        amplitude_modelobstype = ModelObstype(
            obstype=amplitude_obstype, model_unit=self.get_modelunit(), model_band=None
        )

        return data, amplitude_modelobstype


# =============================================================================
# Define era5 default obstypes
# =============================================================================

temp_era5 = ModelObstype(
    obstype=temperature, model_unit="Kelvin", model_band="temperature_2m"
)
temp_era5.set_description(
    "Temperature of air at 2m above the surface of land, sea or in-land waters. 2m temperature is calculated by interpolating between the lowest model level and the Earth's surface, taking account of the atmospheric conditions."
)


pressure_era5 = ModelObstype(
    obstype=pressure, model_unit="pa", model_band="surface_pressure"
)
pressure_era5.set_description(
    "Pressure (force per unit area) of the atmosphere on the surface of land, sea and in-land water. It is a measure of the weight of all the air in a column vertically above the area of the Earth's surface represented at a fixed point. Surface pressure is often used in combination with temperature to calculate air density. The strong variation of pressure with altitude makes it difficult to see the low and high pressure systems over mountainous areas, so mean sea level pressure, rather than surface pressure, is normally used for this purpose. The units of this variable are Pascals (Pa). Surface pressure is often measured in hPa and sometimes is presented in the old units of millibars, mb (1 hPa = 1 mb = 100 Pa).",
)


# create a new obstype that represent the vectorfield of wind
wind_components = Obstype(
    "wind",
    std_unit=wind_speed.get_standard_unit(),
    description=None,
    unit_aliases=wind_speed.units_aliases,
    unit_conversions=wind_speed.conv_table,
)

wind_era5 = ModelObstype_Vectorfield(
    obstype=wind_components,
    model_unit="m/s",
    model_band_u="u_component_of_wind_10m",
    model_band_v="v_component_of_wind_10m",
    amplitude_obstype_name="wind_speed",
    direction_obstype_name="wind_direction",
)
wind_era5.set_description(
    "2D-vector combined 10m windspeed. Care should be taken when comparing this variable with observations, because wind observations vary on small space and time scales and are affected by the local terrain, vegetation and buildings that are represented only on average in the ECMWF Integrated Forecasting System."
)


default_era5_obstypes = {
    temp_era5.name: temp_era5,
    pressure_era5.name: pressure_era5,
    wind_era5.name: wind_era5,
}


class MetobsModelObstypeHandlingError(Exception):
    """Exception raised for errors in the ModelObstype"""

    pass
