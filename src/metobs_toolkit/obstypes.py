import math
from typing import Union
import logging
import numpy as np
import pandas as pd
import pint


from metobs_toolkit.backend_collection.errorclasses import (
    MetObsUnitsIncompatible,
    MetObsUnitUnknown,
    MetObsAdditionError,
)
import metobs_toolkit.backend_collection.printing_collection as printing

# Use logger with name "<metobs_toolkit>"
from metobs_toolkit.backend_collection.loggingmodule import log_entry

logger = logging.getLogger("<metobs_toolkit>")

# ------------------------------------------
#    Setup unit registry
# ------------------------------------------
ureg = pint.UnitRegistry(system="SI")
ureg.formatter.default_format = ".3f"
pint.set_application_registry(ureg)


@log_entry
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


def _join_str_attr(a: str, b: str) -> str:
    if a != b:
        return f"{a} +++ {b}"
    else:
        return a


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

    def _id(self) -> str:
        """A physical unique id.

        In the __add__ methods, if the id of two instances differs, adding is
        a regular concatenation.
        """
        return f"{self.name}_{self.std_unit}"

    def __add__(self, other: "Obstype") -> "Obstype":
        # This function is called when other instances,
        # that hold Obstype 's are joined.

        if self._id() == other._id():
            # combine is valid
            trg_description = _join_str_attr(self.description, other.description)
            trg_orig_name = _join_str_attr(self.original_name, other.original_name)

            # Since trg_orig_unit must be convertible, to pint unit a string join
            # is not applicable. We use the orig_unit from self (it is not used
            # after raw data import, so when __add__ is called, we are passed
            # that stage. )
            if self.original_unit != other.original_unit:
                # If better solution can be made, go ahead.
                trg_orig_unit = self.original_unit
            else:
                trg_orig_unit = self.original_unit

            trg_obs = Obstype(
                obsname=self.name, std_unit=self.std_unit, description=trg_description
            )

            trg_obs.original_name = trg_orig_name
            trg_obs.original_unit = trg_orig_unit

        else:
            raise MetObsAdditionError(
                f"{self} + {other} could not be executes since they do not have the same id ({self._id()} != {other._id()})"
            )

        return trg_obs

    def __eq__(self, other) -> bool:
        """Check equality with another Obstype object."""
        if not isinstance(other, Obstype):
            return False
        return (
            self._name == other._name
            and self.std_unit == other.std_unit
            and self.description == other.description
        )

    def __repr__(self):
        """Instance representation."""
        return f"{type(self).__name__}(id={self._id()})"

    def __str__(self):
        """Text representation."""
        return f"{type(self).__name__} instance of {self.name}"

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

    @property
    def original_name(self):
        """Return the original name of the observation type."""
        return str(self._original_name)

    @property
    def original_unit(self) -> str:
        """Return the original unit as string."""
        return fmt_unit_to_str(self._original_unit)

    @name.setter
    @log_entry
    def name(self, value: str):
        """Set the name of the observation type."""
        self._name = str(value)

    @description.setter
    @log_entry
    def description(self, value: str):
        """Set the description of the observation type."""
        self._description = str(value)

    @original_name.setter
    @log_entry
    def original_name(self, value):
        """Set the original name of the observation type."""
        self._original_name = str(value)

    @original_unit.setter
    @log_entry
    def original_unit(self, value):
        """Set the original unit and check compatibility with standard unit."""

        # If a tempalate is used, with more mapped columns then present in the
        # imported raw data, then these obstypes (without data), have original_unit
        # set to None. This is a valid value
        if (value is None) | (str(value) == "None"):
            self._original_unit = None
            return

        self._original_unit = _fmtunit(value)
        # test if it is a compatible unit wrt the standard unit
        if not self._original_unit.is_compatible_with(self._std_unit):
            raise MetObsUnitsIncompatible(
                f"{self._original_unit} is not compatible with the standard unit ({self.std_unit} of {self}) "
            )

    @log_entry
    def get_compatible_units(self) -> list:
        """
        Get all units compatible with the standard unit.

        Returns
        -------
        list
            List of compatible units as strings.
        """
        logger.debug(
            f"{self.__class__.__name__}.get_compatible_units called for {self}"
        )
        # Get all units related to the same dimension
        compunits = list(ureg.get_compatible_units(self._std_unit.dimensionality))
        # return the units as strings
        return [fmt_unit_to_str(uni) for uni in compunits]

    @log_entry
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

    def _get_plot_y_label(self) -> str:
        """Return a string to represent the vertical axes of a plot."""
        return f"{self.name} ({self.std_unit})"

    def _get_info_core(self) -> str:
        infostr = ""
        infostr += printing.print_fmt_line(f"{self.name} observation with:", 0)
        infostr += printing.print_fmt_line(f"standard unit: {self.std_unit}")
        infostr += printing.print_fmt_line(f"description: {self.description}")
        return infostr

    @log_entry
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
        infostr = ""
        infostr += printing.print_fmt_title("General info of Obstype")
        infostr += self._get_info_core()
        if printout:
            print(infostr)
        else:
            return infostr

    @log_entry
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
        logger.debug(
            f"{self.__class__.__name__}.convert_to_standard_units called for {self}"
        )
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

    def __init__(
        self, obstype: Obstype, model_unit: Union[str, pint.Unit], model_band: str
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
        # Set modelband
        self._model_band = str(model_band)

    def _id(self) -> str:
        """A physical unique id.

        In the __add__ methods, if the id of two instances differs, adding is
        a regular concatenation.
        """
        return f"{super()._id()}_{self.model_band}"

    def __add__(self, other: "ModelObstype") -> "ModelObstype":
        if self._id() != other._id():
            raise MetObsAdditionError(
                f"{self} + {other} could not be executes since they do not have the same id ({self._id()} != {other._id()})"
            )

        # Use the super().__add__ to combine the base Obstype attributes
        combined_obstype = super().__add__(other)

        # Since model_unit must be convertible, to pint unit a string join
        # is not applicable. We use the model_unit from self (it is not used
        # after import, so when __add__ is called, we are passed
        # that stage. )
        if self.model_unit != other.model_unit:
            # Note: this statement is highly unlikely, the same band of a model
            # will not change units!
            # If better solution can be made, go ahead.
            trg_model_unit = self.model_unit
        else:
            trg_model_unit = self.model_unit

        # Create new ModelObstype with combined attributes
        combined = ModelObstype(
            obstype=combined_obstype,
            model_unit=trg_model_unit,
            model_band=self.model_band,  # this is part of the id, so self and other are the same
        )
        return combined

    @property
    def model_unit(self) -> str:
        """Return the model unit as string."""
        return fmt_unit_to_str(self._model_unit)

    @property
    def model_band(self) -> str:
        """Return the model band as string."""
        return str(self._model_band)

    @model_unit.setter
    @log_entry
    def model_unit(self, value: Union[str, pint.Unit]):
        self._model_unit = _fmtunit(value)
        # test if it is a compatible unit wrt the standard unit
        if not self._model_unit.is_compatible_with(self._std_unit):
            raise MetObsUnitsIncompatible(
                f"{self._model_unit} is not compatible with the standard unit ({self.std_unit} of {self}) "
            )

    @model_band.setter
    @log_entry
    def model_band(self, value: str):
        self._model_band = str(value)

    @log_entry
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
        infostr = ""
        infostr += printing.print_fmt_title("General info of ModelObstype")
        infostr += printing.print_fmt_section("Obstype info")
        infostr += super()._get_info_core()
        infostr += printing.print_fmt_section("Model related info")
        infostr += printing.print_fmt_line(f"corresponding bandname: {self.model_band}")
        infostr += printing.print_fmt_line(
            f"original modeldata unit: {self.model_unit}"
        )
        if printout:
            print(infostr)
        else:
            return infostr


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

    def _id(self) -> str:
        """A physical unique id.

        In the __add__ methods, if the id of two instances differs, adding is
        a regular concatenation.
        """
        return f"{super()._id()}_{self.model_band_u}_{self.model_band_v}"

    def __add__(self, other: "ModelObstype_Vectorfield") -> "ModelObstype_Vectorfield":
        if self._id() != other._id():
            raise MetObsAdditionError(
                f"{self} + {other} could not be executes since they do not have the same id ({self._id()} != {other._id()})"
            )

        # Use the super().__add__ to combine the base Obstype attributes
        combined_obstype = super().__add__(other)

        # Since model_unit must be convertible, to pint unit a string join
        # is not applicable. We use the model_unit from self (it is not used
        # after import, so when __add__ is called, we are passed
        # that stage. )
        if self.model_unit != other.model_unit:
            # Note: this statement is highly unlikely, the same band of a model
            # will not change units!
            # If better solution can be made, go ahead.
            trg_model_unit = self.model_unit
        else:
            trg_model_unit = self.model_unit

        # Create new ModelObstype with combined attributes
        combined = ModelObstype_Vectorfield(
            obstype=combined_obstype,
            model_unit=trg_model_unit,
            model_band_u=self.model_band_u,  # in ID
            model_band_v=self._model_band_v,  # in ID
            amplitude_obstype_name=self.amplitude_obstype_name,  # This seems oke
            direction_obstype_name=self.direction_obstype_name,  # This seems oke
        )
        return combined

    @property
    def model_unit(self) -> str:
        """Return the model unit as string."""
        return fmt_unit_to_str(self._model_unit)

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

    @log_entry
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
        infostr = ""
        infostr += printing.print_fmt_title("General info of ModelObstype_Vectorfield")
        infostr += printing.print_fmt_section("Obstype info")
        infostr += super()._get_info_core()
        infostr += printing.print_fmt_section("Model related info")
        infostr += printing.print_fmt_line(f"U-component bandname: {self.model_band_u}")
        infostr += printing.print_fmt_line(f"in {self.model_unit}", 2)
        infostr += printing.print_fmt_line(f"V-component bandname: {self.model_band_v}")
        infostr += printing.print_fmt_line(f"in {self.model_unit}", 2)

        if printout:
            print(infostr)
        else:
            return infostr

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

        @log_entry
        def unit_vector(vector):
            """Returns the unit vector of the vector."""
            return vector / np.linalg.norm(vector)

        @log_entry
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

    @log_entry
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


@log_entry
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
    logger.debug("is_known_unit called")
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
    logger.debug("_fmtunit called")
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


@log_entry
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
