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
        """Initiate an observation type, to link with a GEE dataset band.

        A ModelObstype is specific to a GEE Dataset, and is therefore added
        to a `GeeDynamicModelData` (that facilitates the link with a GEE dataset).

        All methods and attributes are inherited from the `Obstype` class.

        Parameters
        ----------
        obstype : metobs_toolkit.Obstype
            The Obstype that represents the band in the GEE dataset.
        model_unit : str
            The units of the GEE band. This can be found in the details of
            the corresponding GEE dataset. This unit must be known by the obstype (
            add it if it is not known).
        model_band : str
            The name of the band that represents the obstype. This can be found in
            the details of the GEE dataset.

        Returns
        -------
        None

        See also
        ----------
        Obstype: A regular observation type
        ModelObstype_Vectorfield: A vector representation of a ModelObstype.

        Examples
        ---------
        For example, we create a `ModelObstype` for downward solar radiation at the
        surface, to be used with the ERA5-land GEE dataset.

        >>> import metobs_toolkit
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> dataset.gee_datasets['ERA5-land'].modelobstypes
        {'temp': ModelObstype instance of temp (linked to band: temperature_2m), 'pressure': ModelObstype instance of pressure (linked to band: surface_pressure), 'wind': ModelObstype_Vectorfield instance of wind (linked to bands: u_component_of_wind_10m and v_component_of_wind_10m)}

        There is no default solar radiation modeldata `ModelObstype` present for the
        ERA5 `GeeDynamicModelData`. Thus we must create one.

        >>> dataset.obstypes
        {'temp': Obstype instance of temp, 'humidity': Obstype instance of humidity, 'radiation_temp': Obstype instance of radiation_temp, 'pressure': Obstype instance of pressure, 'pressure_at_sea_level': Obstype instance of pressure_at_sea_level, 'precip': Obstype instance of precip, 'precip_sum': Obstype instance of precip_sum, 'wind_speed': Obstype instance of wind_speed, 'wind_gust': Obstype instance of wind_gust, 'wind_direction': Obstype instance of wind_direction}

        We see that there is no default `Obstype` that we can use for
        solar radiation. Therefore we must create a new `Obstype`.

        >>> sol_rad_down_surf = metobs_toolkit.Obstype(
        ...                         obsname="solar_rad_down_at_surface",
        ...                         std_unit='J/m²')

        By looking up the details of ERA5-land (https://developers.google.com/earth-engine/datasets/catalog/ECMWF_ERA5_LAND_HOURLY),
        we find that the corresponding band is *surface_solar_radiation_downwards* and
        the units are in "J/m²" (no coincidence).

        >>> sol_rad_down_for_ERA5 = metobs_toolkit.ModelObstype(
        ...                            obstype=sol_rad_down_surf,
        ...                            model_unit="J/m²",
        ...                            model_band="surface_solar_radiation_downwards")
        >>> sol_rad_down_for_ERA5
        ModelObstype instance of solar_rad_down_at_surface (linked to band: surface_solar_radiation_downwards)

        In practice we add it to the known ModelObstypes of the ERA5 `GeeDynamicModelData`.

        >>> era5_mod = dataset.gee_datasets['ERA5-land']
        >>> era5_mod.add_modelobstype(sol_rad_down_for_ERA5)
        >>> era5_mod.modelobstypes
        {'temp': ModelObstype instance of temp (linked to band: temperature_2m), 'pressure': ModelObstype instance of pressure (linked to band: surface_pressure), 'wind': ModelObstype_Vectorfield instance of wind (linked to bands: u_component_of_wind_10m and v_component_of_wind_10m), 'solar_rad_down_at_surface': ModelObstype instance of solar_rad_down_at_surface (linked to band: surface_solar_radiation_downwards)}

        """
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

        # link with vectorfield
        # The amplitude and direction fields, created by a vectorfield, are
        # Modelobstpes. But the model_band does not realy exist, so special
        # care is needed for them
        self._originates_from_vectorfield = False  #

    def __eq__(self, other):
        is_eq = (
            (self.name == other.name)
            & (self.std_unit == other.std_unit)
            & (self.description == other.description)
            & (self.units_aliases == other.units_aliases)
            & (self.conv_table == other.conv_table)
            & (self.model_unit == other.model_unit)
            & (self.model_band == other.model_band)
            & (self._originates_from_vectorfield == other._originates_from_vectorfield)
        )
        return is_eq

    def _check_validity(self):
        if not self.test_if_unit_is_known(unit_name=self.model_unit):
            raise MetobsModelObstypeHandlingError(
                f"{self.model_unit} is not a known unit of {self}."
            )
        # convert to standard -naming
        self.model_unit = self._get_std_unit_name(self.model_unit)

    def __str__(self):
        return f"{type(self).__name__} instance of {self.name} (linked to band: {self.model_band})"

    def __repr__(self):
        return str(self)

    # =============================================================================
    # Getters
    # =============================================================================

    def get_modelunit(self):
        """Get the (original) unit of the values on the GEE-side."""
        return self.model_unit

    def get_modelband(self):
        """Get the name of the corresponding band."""
        return self.model_band

    def _get_plot_y_label(self):
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
        """Initiate an observation type, to link with GEE dataset bands as vector components.

        A ModelObstype_Vectorfield represents an observationtype, which is not
        directly present in a GEE dataset, but its components are. Wind speed, is
        a common example that is measured by a station, but only the U and V
        components are present in the GEE model.

        When extracting the modeldata, there are two ModelObstypes that are formed:
            * amplitude_field: This is the computed amplitude of the vectorfield
            * direction_field: The orientation of the vectorfield, in degrees
              with the North (if u-v components) as zero and counter-clock-wise
              rotation.

        All methods and attributes are inherited from the `Obstype` class.

        Parameters
        ----------
        obstype : metobs_toolkit.Obstype
            The Obstype that represents the band in the GEE dataset.
        model_unit : str
            The units of the GEE bands. This can be found in the details of
            the corresponding GEE dataset. This unit must be known by the obstype (
            add it if it is not known).
        model_band_u : str
            The name of the band that represents the (U) component. This can be found in
            the details of the GEE dataset.
        model_band_v : str
            The name of the band that represents the (V) component. This can be found in
            the details of the GEE dataset.
        amplitude_obstype_name : str
            The name for the ModelObstype that represents the amplitude of
            the vectorfield.
        dirction_obstype_name : str
            The name for the ModelObstype that represents the direction of the
            vectorfield vectors.

        Returns
        -------
        None

        See also
        ----------
        Obstype: A regular observation type
        ModelObstype_Vectorfield: A vector representation of a ModelObstype.

        Note
        --------
        The U and V component can be defined as any combination of orthogonal
        vectors. There is however no effect on the amplitude of the computed
        vectors, but the zero-position of the vector direction can vary!

        Note
        ----------
        The units of both components of the field must be the same.

        Examples
        ---------

        As example, we create a `ModelObstype_Vectorfield` for 10m wind.
        The windspeed (and direction) is often measured directly, but is stored
        as wind components in a model or GEE dataset.

        >>> import metobs_toolkit
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> dataset.gee_datasets['ERA5-land'].modelobstypes
        {'temp': ModelObstype instance of temp (linked to band: temperature_2m), 'pressure': ModelObstype instance of pressure (linked to band: surface_pressure), 'wind': ModelObstype_Vectorfield instance of wind (linked to bands: u_component_of_wind_10m and v_component_of_wind_10m)}


        Note that there is by default already a `ModelObstype_Vectorfield` for
        the wind. This example is thus purely illustrative.


        We start by creating a (regular) `Obstype` that represents the components
        of the components.

        >>> wind_10m_component = metobs_toolkit.Obstype(
        ...        "wind_10m",
        ...        std_unit='m/s',
        ...        description=' .. component of the 10m wind field ..',
        ...        unit_aliases= {
        ...                    "m/s": ["meters/second", "m/sec"],
        ...                    "km/h": ["kilometers/hour", "kph"],
        ...                    "mph": ["miles/hour"]},
        ...        unit_conversions={
        ...                    "km/h": ["x / 3.6"],
        ...                    "mph": ["x * 0.44704"]})
        >>> wind_10m_component
        Obstype instance of wind_10m

        Now we can create a `ModeObstype_Vectorfield` based on the `wind_10m_component`
        definition.


        >>> wind_10m_era5 = metobs_toolkit.ModelObstype_Vectorfield(
        ...                        obstype=wind_10m_component,
        ...                        model_unit="m/s", # see GEE dataset details
        ...                        model_band_u="u_component_of_wind_10m",  # see GEE dataset details
        ...                        model_band_v="v_component_of_wind_10m",  # see GEE dataset details
        ...                        amplitude_obstype_name="wind_speed",
        ...                        direction_obstype_name="wind_direction")

        >>> wind_10m_era5
        ModelObstype_Vectorfield instance of wind_10m (linked to bands: u_component_of_wind_10m and v_component_of_wind_10m)

        If you want to use it for extracting data, add it to your dataset's Modeldata first.

        >>> dataset.gee_datasets
        {'lcz': GeeStaticModelData instance of lcz  (no metadata has been set) , 'altitude': GeeStaticModelData instance of altitude  (no metadata has been set) , 'worldcover': GeeStaticModelData instance of worldcover  (no metadata has been set) , 'ERA5-land': Empty GeeDynamicModelData instance of ERA5-land }

        >>> dataset.gee_datasets['ERA5-land'].add_modelobstype(wind_10m_era5)
        >>> dataset.gee_datasets['ERA5-land'].get_info()
        Empty GeeDynamicModelData instance of ERA5-land
        ------ Details ---------
        <BLANKLINE>
         * name: ERA5-land
         * location: ECMWF/ERA5_LAND/HOURLY
         * value_type: numeric
         * scale: 2500
         * is_static: False
         * is_image: False
         * is_mosaic: False
         * credentials:
         * time res: 1h
        <BLANKLINE>
         -- Known Modelobstypes --
        <BLANKLINE>
         * temp : ModelObstype instance of temp (linked to band: temperature_2m)
            (conversion: Kelvin --> Celsius)
         * pressure : ModelObstype instance of pressure (linked to band: surface_pressure)
            (conversion: pa --> pa)
         * wind : ModelObstype_Vectorfield instance of wind (linked to bands: u_component_of_wind_10m and v_component_of_wind_10m)
            vectorfield that will be converted to:
              * wind_speed
              * wind_direction
            (conversion: m/s --> m/s)
         * wind_10m : ModelObstype_Vectorfield instance of wind_10m (linked to bands: u_component_of_wind_10m and v_component_of_wind_10m)
            vectorfield that will be converted to:
              * wind_speed
              * wind_direction
            (conversion: m/s --> m/s)
        <BLANKLINE>
         -- Metadata --
        <BLANKLINE>
        No metadf is set.
        <BLANKLINE>
         -- Modeldata --
        <BLANKLINE>
        No model data is set.


        """

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

    def __str__(self):
        return f"{type(self).__name__} instance of {self.name} (linked to bands: {self.get_modelband_u()} and {self.get_modelband_v()})"

    def __repr__(self):
        return str(self)

    def get_modelunit(self):
        return self.model_unit

    def get_modelband_u(self):
        return self.model_band_u

    def get_modelband_v(self):
        return self.model_band_v

    def _get_plot_y_label(self):
        """Return a string to represent the vertical axes of a plot."""
        return f"{self.name} ({self.std_unit})\n originates from {self.original_name}"

    def _compute_angle(self, df):
        """Compute vector direction of 2D vectorfield components.

        The direction column is added to the dataframe and a new ModelObstype,
        representing the angle is returned. The values represent the angles in
        degrees, from north in clockwise rotation.

        Parameters
        ----------
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
            model_band=self._dir_obs_name,  # NOTE: this band does not exist, but column is created with this name by the toolkit
        )
        direction_modelobstype._originates_from_vectorfield = True

        return data, direction_modelobstype

    def compute_amplitude(self, df):
        """Compute amplitude of 2D vectorfield components.

        The amplitude column is added to the dataframe and a new ModelObstype,
        representing the amplitude is returned. All attributes wrt the units are
        inherited from the ModelObstype_vectorfield.

        Parameters
        ------------
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
            obstype=amplitude_obstype,
            model_unit=self.get_modelunit(),
            model_band=self._amp_obs_name,
        )  # NOTE: this band does not exist, but column is created with this name by the toolkit
        amplitude_modelobstype._originates_from_vectorfield = True

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


default_era5_obstypes = [temp_era5, pressure_era5, wind_era5]


class MetobsModelObstypeHandlingError(Exception):
    """Exception raised for errors in the ModelObstype"""

    pass


# =============================================================================
# Docstring test
# =============================================================================
if __name__ == "__main__":
    from metobs_toolkit.doctest_fmt import setup_and_run_doctest

    setup_and_run_doctest()
