import math
from typing import Union
import logging
import numpy as np
import pandas as pd
import pint

# Use logger with name "<metobs_toolkit>"
logger = logging.getLogger("<metobs_toolkit>")

# ------------------------------------------
#    Setup unit registry
# ------------------------------------------
ureg = pint.UnitRegistry(system="SI")
ureg.formatter.default_format = ".3f"
pint.set_application_registry(ureg)


def fmt_unit_to_str(unit) -> str:
    """
    Format a pint unit or quantity to string.

    Parameters
    ----------
    unit : pint.Unit or pint.Quantity
        The unit or quantity to format.

    Returns
    -------
    str
        The formatted unit as string.
    """
    logger.debug(f"{fmt_unit_to_str.__name__} called")
    if isinstance(unit, pint.Unit):
        return str(unit)
    if isinstance(unit, pint.Quantity):
        if unit.magnitude == 1:
            return str(unit.u)
        else:
            # a non-trivial quantity
            return str(unit)
    return str(unit)


class Obstype:
    """
    Class representing an observation type with standard unit and description.

    Parameters
    ----------
    obsname : str
        Name of the observation type.
    std_unit : str or pint.Unit
        Standard unit for the observation type.
    description : str
        Description of the observation type.
    """

    def __init__(self, obsname: str, std_unit: Union[str, pint.Unit], description: str):
        # set name
        self._name = str(obsname)
        # set standard unit
        self._std_unit = _fmtunit(std_unit)
        # set description
        self._description = str(description)
        # open slots
        self._original_name = None
        self._original_unit = None

    def __eq__(self, other) -> bool:
        """Check equality with another Obstype object."""
        if not isinstance(other, Obstype):
            return False
        return (
            self._name == other._name
            and self._std_unit == other._std_unit
            and self._description == other._description
        )

    @property
    def name(self) -> str:
        """Return the name of the observation type."""
        return str(self._name)

    @property
    def std_unit(self) -> str:
        """Return the standard unit as string."""
        return fmt_unit_to_str(self._std_unit)

    @property
    def description(self) -> str:
        """Return the description of the observation type."""
        return str(self._description)

    @description.setter
    def description(self, value: str):
        """Set the description of the observation type."""
        self._description = str(value)

    @property
    def orginal_name(self): 
        """Return the original name of the observation type."""
        return str(self._original_name)

    @orginal_name.setter
    def original_name(self, value): 
        """Set the original name of the observation type."""
        self._original_name = str(value)

    def get_compatible_units(self) -> list:
        """
        Get all units compatible with the standard unit.

        Returns
        -------
        list
            List of compatible units as strings.
        """
        logger.debug(f"{self.__class__.__name__}.get_compatible_units called for {self}")
        # Get all units related to the same dimension
        compunits = list(ureg.get_compatible_units(self._std_unit.dimensionality))
        # return the units as strings
        return [fmt_unit_to_str(uni) for uni in compunits]

    def is_compatible_with(self, other: "Obstype") -> bool:
        """
        Test if the other Obstype is compatible with this one.

        Parameters
        ----------
        other : Obstype
            The other Obstype to test compatibility with.

        Returns
        -------
        bool
            True if compatible, False otherwise.
        """
        logger.debug(f"{self.__class__.__name__}.is_compatible_with called for {self}")
        return self._std_unit.is_compatible_with(other._std_unit)

    @property
    def original_unit(self) -> str:
        """Return the original unit as string."""
        return fmt_unit_to_str(self._original_unit)

    @original_unit.setter
    def original_unit(self, value):
        """Set the original unit and check compatibility with standard unit."""
        self._original_unit = _fmtunit(value)
        # test if it is a compatible unit wrt the standard unit
        if not self._original_unit.is_compatible_with(self._std_unit):
            raise MetObsUnitsIncompatible(
                f"{self._original_unit} is not compatible with the standard unit ({self.std_unit} of {self}) "
            )

    def __repr__(self):
        """Instance representation."""
        return f"{type(self).__name__} instance of {self.name}"

    def __str__(self):
        """Text representation."""
        return f"{type(self).__name__} instance of {self.name}"

    def get_info(self, printout: bool = True) -> Union[None, str]:
        """
        Print or return detailed information of the observation type.

        Parameters
        ----------
        printout : bool, optional
            If True, print the information. If False, return as string. Default is True.

        Returns
        -------
        None or str
            None if printout is True, otherwise the information string.
        """
        logger.debug(f"{self.__class__.__name__}.get_info called for {self}")
        info_str = f"{self.name} observation with: \n \
    * standard unit: {self.std_unit} \n \
    * data column as {self.original_name} in {self.original_unit} \n \
    * description: {self.description}"
        if printout:
            print(info_str)
        else:
            return info_str

    def _get_plot_y_label(self) -> str:
        """Return a string to represent the vertical axes of a plot."""
        return f"{self.name} ({self.std_unit})"

    def convert_to_standard_units(self, input_data, input_unit):
        """
        Convert input data to the standard unit.

        Parameters
        ----------
        input_data : array-like, pd.Series, int, or float
            The data to convert.
        input_unit : str or pint.Unit
            The unit of the input data.

        Returns
        -------
        array-like, pd.Series, int, or float
            The converted data in standard units.
        """
        logger.debug(f"{self.__class__.__name__}.convert_to_standard_units called for {self}")
        # format input unit
        input_unit = _fmtunit(input_unit)
        # Test if inputunit is compatible with std unit
        if not input_unit.is_compatible_with(self._std_unit):
            raise MetObsUnitUnknown(
                f"{input_unit} is not compatible with the standard unit ({self.std_unit} of {self}) "
            )
        # convert data
        return convert_units(
            records=input_data, cur_unit=input_unit, trg_unit=self._std_unit
        )


class ModelObstype(Obstype):
    """
    Class representing a model observation type, extending Obstype.

    Parameters
    ----------
    obstype : Obstype
        The base Obstype.
    model_unit : str or pint.Unit
        The model's unit.
    model_band : str
        The model's band name.
    """

    def __init__(self, obstype: Obstype, model_unit: Union[str, pint.Unit], model_band: str):
        
        # set regular obstype
        super().__init__(
            obsname=obstype.name,
            std_unit=obstype.std_unit,
            description=obstype.description,
        )
        # Set modelunit
        self._model_unit = _fmtunit(model_unit)
        # test if it is a compatible unit wrt the standard unit
        if not self._model_unit.is_compatible_with(self._std_unit):
            raise MetObsUnitUnknown(
                f"{self._model_unit} is not compatible with the standard unit ({self.std_unit} of {self}) "
            )
        # Set modelband
        self._model_band = str(model_band)

    @property
    def model_unit(self) -> str:
        """Return the model unit as string."""
        return str(self._model_unit)

    @property
    def model_band(self) -> str:
        """Return the model band as string."""
        return str(self._model_band)

    def get_info(self, printout: bool = True) -> Union[None, str]:
        """
        Print or return detailed information of the model observation type.

        Parameters
        ----------
        printout : bool, optional
            If True, print the information. If False, return as string. Default is True.

        Returns
        -------
        None or str
            None if printout is True, otherwise the information string.
        """
        logger.debug(f"{self.__class__.__name__}.get_info called for {self}")
        standardinfo = super().get_info(False)
        standardinfo += (
            f"\n     * corresponding bandname {self.model_band} in {self.model_unit}"
        )
        if printout:
            print(standardinfo)
        else:
            return standardinfo


class ModelObstype_Vectorfield(Obstype):
    """
    Class representing a model observation type for vector fields.

    Parameters
    ----------
    obstype : Obstype
        The base Obstype.
    model_unit : str or pint.Unit
        The model's unit.
    model_band_u : str
        The model's U-component band name.
    model_band_v : str
        The model's V-component band name.
    amplitude_obstype_name : str
        Name for the amplitude obstype.
    direction_obstype_name : str
        Name for the direction obstype.
    """

    def __init__(
        self,
        obstype: Obstype,
        model_unit: Union[str, pint.Unit],
        model_band_u: str,
        model_band_v: str,
        amplitude_obstype_name: str,
        direction_obstype_name: str,
    ):
        # set regular obstype
        super().__init__(
            obsname=obstype.name,
            std_unit=obstype.std_unit,
            description=obstype.description,
        )
        # Set modelunit
        self._model_unit = _fmtunit(model_unit)
        # test if it is a compatible unit wrt the standard unit
        if not self._model_unit.is_compatible_with(self._std_unit):
            raise MetObsUnitUnknown(
                f"{self._model_unit} is not compatible with the standard unit ({self.std_unit} of {self}) "
            )
        # Set bandnames
        self._model_band_u = str(model_band_u)
        self._model_band_v = str(model_band_v)
        self._amp_obs_name = str(amplitude_obstype_name)
        self._dir_obs_name = str(direction_obstype_name)

    @property
    def model_unit(self) -> str:
        """Return the model unit as string."""
        return str(self._model_unit)

    @property
    def model_band_u(self) -> str:
        """Return the model U-component band as string."""
        return str(self._model_band_u)

    @property
    def model_band_v(self) -> str:
        """Return the model V-component band as string."""
        return str(self._model_band_v)

    @property
    def amplitude_obstype_name(self) -> str:
        """Return the amplitude obstype name as string."""
        return str(self._amp_obs_name)

    @property
    def direction_obstype_name(self) -> str:
        """Return the direction obstype name as string."""
        return str(self._dir_obs_name)

    def get_info(self, printout: bool = True) -> Union[None, str]: 
        """
        Print or return detailed information of the vectorfield model observation type.

        Parameters
        ----------
        printout : bool, optional
            If True, print the information. If False, return as string. Default is True.

        Returns
        -------
        None or str
            None if printout is True, otherwise the information string.
        """
        logger.debug(f"{self.__class__.__name__}.get_info called for {self}")
        standardinfo = super().get_info(False)
        standardinfo += f"\n     * U-component bandname {self.model_band_u} in {self.model_unit} \n \
    * V-component bandname {self.model_band_v} in {self.model_unit}"
        if printout:
            print(standardinfo)
        else:
            return standardinfo

    def _get_plot_y_label(self) -> str:
        """Return a string to represent the vertical axes of a plot."""
        return f"{self.name} ({self.std_unit})\n originates from {self.original_name}"

    def _compute_angle(self, df: pd.DataFrame) -> Union[pd.DataFrame, ModelObstype]:
        """
        Compute vector direction of 2D vectorfield components.

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
            The (scalar) ModelObstype representation of the angles.
        """
        logger.debug(f"{self.__class__.__name__}._compute_angle called for {self}")

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

        u_column = self.model_band_u
        v_column = self.model_band_v

        data = df.apply(lambda x: angle_between(x[u_column], x[v_column]), axis=1)

        # Create a new obstype for the direction
        direction_obstype = Obstype(
            obsname=self.direction_obstype_name,
            std_unit=ureg.degree,
            description=f"Direction of 2D-vector of {self.name} components.",
        )
        # convert to model obstype
        direction_modelobstype = ModelObstype(
            obstype=direction_obstype,
            model_unit=ureg.degree,  # indep of units
            model_band=self.direction_obstype_name,  # NOTE: this band does not exist, but column is created with this name by the toolkit
        )
        direction_modelobstype._originates_from_vectorfield = True

        return data, direction_modelobstype

    def compute_amplitude(self, df: pd.DataFrame) -> Union[pd.DataFrame, ModelObstype]:
        """
        Compute amplitude of 2D vectorfield components.

        The amplitude column is added to the dataframe and a new ModelObstype,
        representing the amplitude is returned. All attributes with respect to the units are
        inherited from the ModelObstype_Vectorfield.

        Parameters
        ----------
        df : pandas.DataFrame
            The dataframe with the vector components present as columns.

        Returns
        -------
        data : pandas.DataFrame
            The df with an extra column representing the amplitudes.
        amplitude_obstype : ModelObstype
            The (scalar) ModelObstype representation of the amplitudes.
        """
        logger.debug(f"{self.__class__.__name__}.compute_amplitude called for {self}")
        # Compute the data
        data = ((df[self.model_band_u].pow(2)) + (df[self.model_band_v].pow(2))).pow(
            1.0 / 2
        )

        # Create a new Obstype for the amplitude
        amplitude_obstype = Obstype(
            obsname=self.amplitude_obstype_name,
            std_unit=self.std_unit,
            description=f"2D-vector amplitude of {self.name} components.",
        )

        # convert to model obstype
        amplitude_modelobstype = ModelObstype(
            obstype=amplitude_obstype,
            model_unit=self.model_unit,
            model_band=self.amplitude_obstype_name,
        )  # NOTE: this band does not exist, but column is created with this name by the toolkit
        amplitude_modelobstype._originates_from_vectorfield = True

        return data, amplitude_modelobstype


# ------------------------------------------
#    Helpers
# ------------------------------------------


def is_known_unit(unit: str) -> bool:
    """
    Check if a unit string is known to the unit registry.

    Parameters
    ----------
    unit : str
        The unit string to check.

    Returns
    -------
    bool
        True if known, False otherwise.
    """
    logger.debug(f"is_known_unit called")
    try:
        ureg.parse_expression(unit)
        return True
    except pint.errors.UndefinedUnitError:
        return False


def _fmtunit(value) -> pint.Unit:
    """
    Convert unit input to pint.Unit.

    Parameters
    ----------
    value : str, pint.Unit, or pint.Quantity
        The value to convert.

    Returns
    -------
    pint.Unit
        The formatted unit.

    Raises
    ------
    MetObsUnitUnknown
        If the value cannot be converted to a known unit.
    """
    logger.debug(f"_fmtunit called")
    if isinstance(value, pint.Unit):
        return value
    elif isinstance(value, pint.Quantity):
        if value.magnitude == 1:
            return value.u
        else:
            raise MetObsUnitUnknown(
                f"{value} is a pint.Quantity with non-1 magnitude and cannot be converted to a unit."
            )
    elif isinstance(value, str):
        if not is_known_unit(value):
            raise MetObsUnitUnknown(
                f"{value} is not a known unit. (See https://github.com/hgrecco/pint/blob/master/pint/default_en.txt for all known units.)"
            )
        else:
            return ureg.parse_expression(value)
    else:
        raise MetObsUnitUnknown(f"{value} is not a string or pint.Unit.")


def convert_units(records, cur_unit, trg_unit):
    """
    Convert records from one unit to another.

    Parameters
    ----------
    records : pd.Series, np.ndarray, int, or float
        The data to convert.
    cur_unit : str or pint.Unit
        The current unit of the data.
    trg_unit : str or pint.Unit
        The target unit to convert to.

    Returns
    -------
    pd.Series, np.ndarray, int, or float
        The converted data.
    """
    logger.debug(f"Converting data from {cur_unit} --> {trg_unit} ")
    if isinstance(records, pd.Series):
        trgvalues = ureg.Quantity(records.to_numpy(), cur_unit).to(trg_unit)
        return pd.Series(
            index=records.index, data=trgvalues.magnitude, name=records.name
        )
    elif isinstance(records, np.ndarray):
        return ureg.Quantity(records, cur_unit).to(trg_unit).magnitude
    elif isinstance(records, int):
        return ureg.Quantity(records, cur_unit).to(trg_unit).magnitude
    elif isinstance(records, float):
        return ureg.Quantity(records, cur_unit).to(trg_unit).magnitude
    else:
        raise NotImplementedError(f"{records} is not a supported input type.")


# ------------------------------------------
#    Errors
# ------------------------------------------
class MetObsUnitsIncompatible(Exception):
    """Raised when an incompatible unit is set."""
    pass

class MetObsUnitUnknown(Exception):
    """Raised when an invalid unit is set."""
    pass

# ------------------------------------------
#    Default obstypes
# ------------------------------------------

temperature = Obstype(
    obsname="temp", std_unit=ureg.degC, description="2m - temperature"
)

humidity = Obstype(
    obsname="humidity",
    std_unit=ureg.percent,
    description="2m - relative humidity",
)

radiation_temp = Obstype(
    obsname="radiation_temp",
    std_unit=ureg.degC,
    description="2m - Black globe",
)

pressure = Obstype(
    obsname="pressure",
    std_unit=ureg.hectopascal,
    description="atmospheric pressure (at station)",
)

pressure_at_sea_level = Obstype(
    obsname="pressure_at_sea_level",
    std_unit=ureg.hectopascal,
    description="atmospheric pressure (at sea level)",
)

precip = Obstype(
    obsname="precip",
    std_unit=ureg.mm / (ureg.meter * ureg.meter),
    description="precipitation intensity",
)

precip_sum = Obstype(
    obsname="precip_sum",
    std_unit=ureg.mm / (ureg.meter * ureg.meter),
    description="Cummulated precipitation",
)
wind_speed = Obstype(
    obsname="wind_speed",
    std_unit=ureg.meter / ureg.second,
    description="wind speed",
)

windgust = Obstype(
    obsname="wind_gust",
    std_unit=ureg.meter / ureg.second,
    description="wind gust",
)

wind_direction = Obstype(
    obsname="wind_direction",
    std_unit=ureg.degree,
    description="wind direction",
)

tlk_obstypes = {
    "temp": temperature,
    "humidity": humidity,
    "radiation_temp": radiation_temp,
    "pressure": pressure,
    "pressure_at_sea_level": pressure_at_sea_level,
    "precip": precip,
    "precip_sum": precip_sum,
    "wind_speed": wind_speed,
    "wind_gust": windgust,
    "wind_direction": wind_direction,
}

temp_era5 = ModelObstype(
    obstype=temperature, model_unit=ureg.degK, model_band="temperature_2m"
)
temp_era5.description = "Temperature of air at 2m above the surface of land, sea or in-land waters. 2m temperature is calculated by interpolating between the lowest model level and the Earth's surface, taking account of the atmospheric conditions."

pressure_era5 = ModelObstype(
    obstype=pressure, model_unit="pascal", model_band="surface_pressure"
)
pressure_era5.description = (
    "Pressure (force per unit area) of the atmosphere on the surface of land, sea and in-land water. It is a measure of the weight of all the air in a column vertically above the area of the Earth's surface represented at a fixed point. Surface pressure is often used in combination with temperature to calculate air density. The strong variation of pressure with altitude makes it difficult to see the low and high pressure systems over mountainous areas, so mean sea level pressure, rather than surface pressure, is normally used for this purpose. The units of this variable are Pascals (Pa). Surface pressure is often measured in hPa and sometimes is presented in the old units of millibars, mb (1 hPa = 1 mb = 100 Pa).",
)

# create a new obstype that represent the vectorfield of wind
wind_components = Obstype(
    "wind",
    std_unit=wind_speed.std_unit,
    description=None,
)

wind_era5 = ModelObstype_Vectorfield(
    obstype=wind_components,
    model_unit=ureg.meter / ureg.second,
    model_band_u="u_component_of_wind_10m",
    model_band_v="v_component_of_wind_10m",
    amplitude_obstype_name="wind_speed",
    direction_obstype_name="wind_direction",
)
wind_era5.description = "2D-vector combined 10m windspeed. Care should be taken when comparing this variable with observations, because wind observations vary on small space and time scales and are affected by the local terrain, vegetation and buildings that are represented only on average in the ECMWF Integrated Forecasting System."

default_era5_obstypes = [temp_era5, pressure_era5, wind_era5]
