
Creating a new observation type
==================================

Observation types for Datasets
--------------------------------

The toolkit comes with a set of predefined observation types. Each observation type has a standard-toolkit-unit,
this is the unit the toolkit will store and display the values.

An overview can be found on `this <./template_mapping.html#toolkit-standards>`_ page.

Each observation type is represented by an instance of the :py:meth:`Obstype<metobs_toolkit.obstypes.Obstype>` class.

As an example, here is the definition of the temperature observation type:

.. code-block:: python

   temperature = Obstype(obsname='temp', #The name of the observation type
                         std_unit= 'Celsius', #The standard unit
                         description="2m - temperature", #A more detailed description (optional)
                         unit_aliases={
                            # Common units and a list of aliases for them.
                             'Celsius': ['celsius', '°C', '°c', 'celcius', 'Celcius'],
                             'Kelvin': ['K', 'kelvin'],
                             'Farenheit': ['farenheit']},
                            # Conversion schemes for common units to the standard unit.
                         unit_conversions={
                             'Kelvin': ["x - 273.15"], #result is in tlk_std_units (aka Celsius)
                             'Farenheit' : ["x-32.0", "x/1.8"]}, # -->execute from left to write  = (x-32)/1.8},
                         )

Similar as this example a user can create a new observation type and add it to a :py:meth:`Dataset<metobs_toolkit.dataset.Dataset>`,
using the :py:meth:`add_new_observationtype()<metobs_toolkit.dataset.Dataset.add_new_observationtype>` method.

.. code-block:: python

   import metobs_toolkit

   #create an new observationtype
   wind_component_east = metobs_toolkit.Obstype(
                         obsname='wind_u_comp', #The name of the observation type
                         std_unit= 'm/s', #The standard unit
                         description="2m - u component of the wind (5min averages)", #A more detailed description (optional)
                         unit_aliases={
                            # Common units and a list of aliases for them.
                             'm/s': ['meter/s'],
                            # Conversion schemes for common units to the standard unit.
                         unit_conversions={'km/s': ["x / 3.6"]} #result is in tlk_std_units (aka m/s)
                         )

   #add your observation type to a dataset
   your_dataset = metobs_toolkit.Dataset()
   your_dataset.add_new_observationtype(Obstype=wind_component_east)

   # Now you can import a datafile with wind_u_comp data!


If you want to add a new unit to an existing observation type you can do so by
using the :py:meth:`add_new_unit()<metobs_toolkit.dataset.Dataset.add_new_unit>` method.


Observation types for (ERA5) Modeldata
----------------------------------------
Modeldata objects also holds a similar set of observation types. But in addition
to the observation types stored in the Dataset, extra information is stored
on where which (ERA5) band and unit the observation type represents. Here is an
example on how to create a new observation type for a :py:meth:`Modeldata<metobs_toolkit.modeldata.Modeldata>` instance.

.. code-block:: python

   import metobs_toolkit

   #create an new observationtype
   wind_component_east = metobs_toolkit.Obstype(
                         obsname='wind_u_comp', #The name of the observation type
                         std_unit= 'm/s', #The standard unit
                         description="10m - east component of the wind ", #A more detailed description (optional)
                         unit_aliases={
                            # Common units and a list of aliases for them.
                             'm/s': ['meter/s'],
                            # Conversion schemes for common units to the standard unit.
                         unit_conversions={'km/s': ["x / 3.6"]} #result is in tlk_std_units (aka m/s)
                         )
   # create a modeldata instance
   model_data = metobs_toolkit.Modeldata("ERA5_hourly")

   # add new obstype to model_data
   model_data.add_obstype(Obstype=wind_component_east,
                          bandname='u_component_of_wind_10m', #See: https://developers.google.com/earth-engine/datasets/catalog/ECMWF_ERA5_LAND_HOURLY#bands
                          band_units='m/s',
                          )

   # Collect the U-wind component for your stations:
   model_data = your_dataset.get_modeldata(modeldata=model_data,
                                           obstype = 'wind_u_comp')
