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
from metobs_toolkit.obstypes import Obstype
from metobs_toolkit.obstypes import expression_calculator

from metobs_toolkit.obstypes import temperature, pressure, wind, direction_aliases

logger = logging.getLogger(__name__)


class ModelObstype(Obstype):
    """Extension of the Obstype class specific for the obstypes of Modeldata."""

    def __init__(self, obstype, band_name, band_unit, band_description):
        """TODO fix update docstring
        Initiate an Modelobservation type.

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
        super().__init__(
            obsname=obstype.name,
            std_unit=obstype.std_unit,
            description=obstype.description,
            unit_aliases=obstype.units_aliases,
            unit_conversions=obstype.conv_table,
        )

        # self.modl_equi_dict = model_equivalent_dict
        self.band_name = str(band_name)
        self.band_unit = str(band_unit)
        self.band_desc = str(band_description)
        self._is_valid()

    def __repr__(self):
        """Instance representation."""
        return f"ModelObstype instance of {self.name}"

    def __str__(self):
        """Text representation."""
        return f"ModelObstype instance of {self.name}"

    def get_info(self):
        """Print out detailed information of the observation type.

        Returns
        -------
        None.

        """
        databands = {key: item["name"] for key, item in self.modl_equi_dict.items()}
        info_str = f"{self.name} observation with: \n \
    * Known datasetsbands: {databands} \n \
    * standard unit: {self.std_unit} \n \
    * description: {self.description} \n \
    * conversions to known units: {self.conv_table} \n"
        print(info_str)

    # def get_mapped_datasets(self):
    #     """Return all gee datasets with a representing band for this obstype."""
    #     return list(self.modl_equi_dict.keys())
    def get_band_description(self):
        """Return the band description."""
        return str(self.band_desc)

    def get_bandname(self):
        """Return the representing bandname of the obstype"""
        return str(self.band_name)

    def get_bandname_mapper(self):
        """Return the representing bandname with tlk standard name as a dict."""
        return {self.get_bandname(): self.name}

    def get_plot_y_label(self):
        """Return a string to represent the vertical axes of a plot."""
        return f"{self.name} ({self.std_unit}) \n {self.get_bandname()}"

    def get_bandunit(self):
        """Return the units of the representing bandname of the obstype."""
        return str(self.band_unit)

    def _is_valid(self):
        """Test if all attributes are valid among each other."""
        # check if the band unit is a knonw unit
        assert (
            self.band_unit in self.get_all_units()
        ), f"{self.band_unit} is not a known unit of {self}"


class ModelObstype_Vectorfield(Obstype):
    def __init__(
        self,
        obstype,
        band_name_u,
        band_unit_u,
        band_description_u,
        band_name_v,
        band_unit_v,
        band_description_v,
    ):

        super().__init__(
            obsname=obstype.name,
            std_unit=obstype.std_unit,
            description=obstype.description,
            unit_aliases=obstype.units_aliases,
            unit_conversions=obstype.conv_table,
        )

        # mod_comp_dict = {}
        # for geedataset in u_comp_model_equivalent_dict.keys():
        #     mod_comp_dict[geedataset] = {
        #         "u_comp": u_comp_model_equivalent_dict[geedataset],
        #         "v_comp": v_comp_model_equivalent_dict[geedataset],
        #     }

        self.band_unit_u = str(band_unit_u)
        self.band_unit_v = str(band_unit_v)

        self.band_name_u = str(band_name_u)
        self.band_name_v = str(band_name_v)

        self.band_desc_u = str(band_description_u)
        self.band_desc_v = str(band_description_v)

        # self.modl_comp_dict = mod_comp_dict
        self._is_valid()

    def __repr__(self):
        """Instance representation."""
        return f"ModelObstype_Vectorfield instance of {self.name}"

    def __str__(self):
        """Text representation."""
        return f"ModelObstype_Vectorfield instance of {self.name}"

    def get_info(self):
        """Print out detailed information of the observation type.

        Returns
        -------
        None.

        """
        info_str = f"{self.name} observation with: \n \
    * Known Vector-East-component databand: {self.band_name_u} \n \
    * Known Vector-North-component datasetsbands: {self.band_name_v} \n \
    * standard unit: {self.std_unit} \n \
    * description: {self.description} \n \
    * conversions to known units: {self.conv_table} \n"
        print(info_str)

    def get_bandname_mapper(self):
        """Return the representing bandname with tlk standard name as a dict."""
        mapper = {
            self.band_name_u: f"u_comp_{self.name}",
            self.band_name_v: f"v_comp_{self.name}",
        }

        return mapper

    def get_bandunit(self):
        """Return the units of the representing bandname of the obstype from a given gee dataset."""
        # u and v comp must have the same units, this is tested in the _is_valid()
        return self.band_unit_u

    def get_plot_y_label(self):
        """Return a string to represent the vertical axes of a plot."""
        return f"{self.name} ({self.std_unit}) \n {self.band_name_u} and {self.band_name_v}"

    def get_u_column(self):
        return f"u_comp_{self.name}"

    def get_v_column(self):
        return f"v_comp_{self.name}"

    def create_the_scalar_modelobstypes(self):
        """Create regular ModelObstypes for the u and v components"""

        u_obstype = ModelObstype(
            obstype=Obstype(
                obsname=self.get_bandname_mapper()[self.band_name_u],
                std_unit=self.std_unit,
                description=self.description,
                unit_aliases=self.units_aliases,
                unit_conversions=self.conv_table,
            ),
            band_name=self.band_name_u,
            band_unit=self.band_unit_u,
            band_description=self.band_desc_u,
        )

        v_obstype = ModelObstype(
            obstype=Obstype(
                obsname=self.get_bandname_mapper()[self.band_name_v],
                std_unit=self.std_unit,
                description=self.description,
                unit_aliases=self.units_aliases,
                unit_conversions=self.conv_table,
            ),
            band_name=self.band_name_v,
            band_unit=self.band_unit_v,
            band_description=self.band_desc_v,
        )

        return u_obstype, v_obstype

    # def add_new_band(
    #     self,
    #     mapname,
    #     bandname_u_comp,
    #     bandname_v_comp,
    #     bandunit,
    #     band_desc_u_comp=None,
    #     band_desc_v_comp=None,
    # ):
    #     """Add a new representing dataset/bandname to the obstype.

    #     Parameters
    #     ----------
    #     mapname : str
    #         name of the known gee dataset.
    #     bandname_u_comp : str
    #         the name of the representing the Eastwards component band.
    #     bandname_v_comp : str
    #         the name of the representing the Northwards component band.
    #     bandunit : str
    #         the unit of the representing bands.
    #     band_desc_u_comp : str, optional
    #         A detailed description of the Eastwards component of the band.
    #     band_desc_v_comp : str, optional
    #         A detailed description of the Northwards component of the band.

    #     Returns
    #     -------
    #     None.

    #     """
    #     # test if banunit is valid
    #     if not self.test_if_unit_is_known(bandunit):
    #         sys.exit(f"{bandunit} is an unknown unit for the {self.name} obstype.")

    #     if mapname in self.modl_comp_dict.keys():
    #         # check if band is already knonw
    #         logger.debug(f"Update {bandname} of (known) map: {mapname}")
    #     else:
    #         logger.debug(f"Add new map: {mapname} with band: {bandname}.")

    #     self.modl_comp_dict[mapname] = {}
    #     self.modl_comp_dict[mapname]["u_comp"] = {
    #         "name": str(bandname_u_comp),
    #         "units": str(bandunit),
    #         "band_desc": str(band_desc_u_comp),
    #     }
    #     self.modl_comp_dict[mapname]["v_comp"] = {
    #         "name": str(bandname_v_comp),
    #         "units": str(bandunit),
    #         "band_desc": str(band_desc_v_comp),
    #     }

    def _is_valid(self):
        """Test if all attributes are valid among each other."""
        # check if the band unit is a knonw unit
        assert (
            self.band_unit_u in self.get_all_units()
        ), f"{self.band_unit_u} (U- band unit) is not a known unit of {self}"
        assert (
            self.band_unit_v in self.get_all_units()
        ), f"{self.band_unit_v} (V- band unit) is not a known unit of {self}"

        # Check if band units are equal
        assert (
            self.band_unit_u == self.band_unit_v
        ), f"The band units of the U and V components are not equal: {self.band_unit_u} != {self.band_unit_v}"

    def convert_to_standard_units(self, input_df, input_unit):
        """Convert data from a known unit to the standard unit.

        The data must be a pandas dataframe with both the u and v component
        prensent as columns.

        Parameters
        ----------
        input_data : pandas.DataFrame
            The dataframe containig the u and v data as columns with the name
            as definde by the get_u_column and get_v_column methods.
        input_unit : str
            The known unit the inputdata is in.

        Returns
        -------
        data_u_component :  pandas.Series
            The u component of the data in standard units.
        data_v_component : pandas.Series
            The v component of the data in standard units.

        """
        # check if input unit is known
        assert self.test_if_unit_is_known(
            input_unit
        ), f"{input_unit} is an unknown unit for {self.name}. No coversion possible!"

        # Get conversion
        std_unit_name = self._get_std_unit_name(input_unit)
        if std_unit_name == self.std_unit:
            # No conversion needed because already the standard unit
            return input_df[self.get_u_column()], input_df[self.get_v_column()]

        conv_expr_list = self.conv_table[std_unit_name]

        # covert data u component
        data_u = input_df[self.get_u_column()]
        data_v = input_df[self.get_v_column()]
        for conv in conv_expr_list:
            data_u = expression_calculator(conv, data_u)
            data_v = expression_calculator(conv, data_v)

        return data_u, data_v


#%% New obs creator functions
def compute_amplitude(modelobs_vectorfield, df):
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
        (df[modelobs_vectorfield.get_u_column()].pow(2))
        + (df[modelobs_vectorfield.get_v_column()].pow(2))
    ).pow(1.0 / 2)

    # Create a new obstype for the amplitude
    amplitude_obstype = Obstype(
        obsname=f"{modelobs_vectorfield.name}_amplitude",
        std_unit=modelobs_vectorfield.std_unit,
        description=f"2D-vector amplitde of {modelobs_vectorfield.name} components.",
        unit_aliases=modelobs_vectorfield.units_aliases,
        unit_conversions=modelobs_vectorfield.conv_table,
    )

    # Create a Modelobstype
    amp_obstype = ModelObstype(
        obstype=amplitude_obstype,
        band_name=f"dummy for {modelobs_vectorfield.band_name_u} and {modelobs_vectorfield.band_name_v}",
        band_unit=modelobs_vectorfield.get_bandunit(),
        band_description=f"Pytagorean amplitude of {modelobs_vectorfield.band_name_u} and {modelobs_vectorfield.band_name_v}",
    )

    return data, amp_obstype


def compute_angle(modelobs_vectorfield, df):
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

    u_column = modelobs_vectorfield.get_u_column()
    v_column = modelobs_vectorfield.get_v_column()

    data = df.apply(lambda x: angle_between(x[u_column], x[v_column]), axis=1)

    # Create a new obstype for the amplitude
    direction_obstype = Obstype(
        obsname=f"{modelobs_vectorfield.name}_direction",
        std_unit="° from north (CW)",
        description=f"Direction of 2D-vector of {modelobs_vectorfield.name} components.",
        unit_aliases=direction_aliases,
        unit_conversions={},
    )

    # Create a Modelobstype
    dir_obstype = ModelObstype(
        obstype=direction_obstype,
        band_name=f"dummy for {modelobs_vectorfield.band_name_u} and {modelobs_vectorfield.band_name_v}",
        band_unit="° from north (CW)",
        band_description=f"direction in ° from north counter-clockwise of {modelobs_vectorfield.band_name_u} and {modelobs_vectorfield.band_name_v}",
    )

    return data, dir_obstype


# =============================================================================
# Define default ERA5 obstypes
# =============================================================================

temp_era5 = ModelObstype(
    obstype=temperature,
    band_name="temperature_2m",
    band_unit="Kelvin",
    band_description="Temperature of air at 2m above the surface of land, sea or in-land waters. 2m temperature is calculated by interpolating between the lowest model level and the Earth's surface, taking account of the atmospheric conditions.",
)


pressure_era5 = ModelObstype(
    obstype=pressure,
    band_name="surface_pressure",
    band_unit="pa",
    band_description="Pressure (force per unit area) of the atmosphere on the surface of land, sea and in-land water. It is a measure of the weight of all the air in a column vertically above the area of the Earth's surface represented at a fixed point. Surface pressure is often used in combination with temperature to calculate air density. The strong variation of pressure with altitude makes it difficult to see the low and high pressure systems over mountainous areas, so mean sea level pressure, rather than surface pressure, is normally used for this purpose. The units of this variable are Pascals (Pa). Surface pressure is often measured in hPa and sometimes is presented in the old units of millibars, mb (1 hPa = 1 mb = 100 Pa).",
)


# Special obstypes
wind.name = "wind"  # otherwise it is windspeed, which is confusing for vectorfield
wind_era5 = ModelObstype_Vectorfield(
    obstype=wind,
    band_name_u="u_component_of_wind_10m",
    band_unit_u="m/s",
    band_description_u="Eastward component of the 10m wind. It is the horizontal speed of air moving towards the east, at a height of ten meters above the surface of the Earth, in meters per second. Care should be taken when comparing this variable with observations, because wind observations vary on small space and time scales and are affected by the local terrain, vegetation and buildings that are represented only on average in the ECMWF Integrated Forecasting System. This variable can be combined with the V component of 10m wind to give the speed and direction of the horizontal 10m wind.",
    band_name_v="v_component_of_wind_10m",
    band_unit_v="m/s",
    band_description_v="Northward component of the 10m wind. It is the horizontal speed of air moving towards the north, at a height of ten meters above the surface of the Earth, in meters per second. Care should be taken when comparing this variable with observations, because wind observations vary on small space and time scales and are affected by the local terrain, vegetation and buildings that are represented only on average in the ECMWF Integrated Forecasting System. This variable can be combined with the U component of 10m wind to give the speed and direction of the horizontal 10m wind.",
)


# =============================================================================
# Create obstype dict
# =============================================================================
era5_default_model_obstypes = {
    "temp": temp_era5,
    "pressure": pressure_era5,
    "wind": wind_era5,
}
