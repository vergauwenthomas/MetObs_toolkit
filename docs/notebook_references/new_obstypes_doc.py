#!/usr/bin/env python
# coding: utf-8

# # Creating a new observation type

# ## Observation types for Datasets
# The toolkit comes with a set of predefined observation types. Each observation type has a standard-toolkit-unit,
# this is the unit the toolkit will store and display the values.
#
# An overview can be found on [Mapping to the toolkit](./template_doc.html#Mapping to the toolkit) page.
#
# Each observation type is represented by an instance of the `metobs_toolkit.Obstype` class.
#
# As an example, here is the definition of the temperature observation type:

# In[5]:


import metobs_toolkit

temperature = metobs_toolkit.Obstype(
    obsname="temp",  # The name of the observation type
    std_unit="Celsius",  # The standard unit
    description="2m - temperature",  # A more detailed description (optional)
    unit_aliases={
        # Common units and a list of aliases for them.
        "Celsius": ["celsius", "°C", "°c", "celcius", "Celcius"],
        "Kelvin": ["K", "kelvin"],
        "Farenheit": ["farenheit"],
    },
    # Conversion schemes for common units to the standard unit.
    unit_conversions={
        "Kelvin": ["x - 273.15"],  # result is in tlk_std_units (aka Celsius)
        "Farenheit": ["x-32.0", "x/1.8"],
    },  # -->execute from left to write  = (x-32)/1.8},
)

temperature


# You can use `Obstype.get_info()` to print out an overview of the observation.

# In[6]:


temperature.get_info()


# In the same manner, we can create a new observationtype by using the `Dataset.add_new_observationtype()` method.

# In[7]:


import metobs_toolkit

# create an new observation type
wind_component_east = metobs_toolkit.Obstype(
    obsname="wind_u_comp",  # The name of the observation type
    std_unit="m/s",  # The standard unit
    description="2m - u component of the wind (5min averages)",  # A more detailed description (optional)
    unit_aliases={
        # Common units and a list of aliases for them.
        "m/s": ["meter/s"]
    },
    # Conversion schemes for common units to the standard unit.
    unit_conversions={"km/s": ["x / 3.6"]},  # result is in tlk_std_units (aka m/s)
)

wind_component_east.get_info()


# add your observation type to a dataset
your_dataset = metobs_toolkit.Dataset()
your_dataset.add_new_observationtype(Obstype=wind_component_east)


# If you want to add a new unit to an existing observation type you can do so by using the `Dataset.add_new_unit()` method.
