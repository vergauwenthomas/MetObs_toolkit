***********************
Mapping to the toolkit
***********************

The MetObs-toolkit uses standard names and formats for your data. To use the toolkit,
your observational data must be converted to the toolkit standards this is referred to as **mapping**.

To specify how the mapping must be done a **template** is used. This template contains
all the information on how to convert your tabular data to the toolkit standards.
Since the structure of data files differs for different networks, this template is
unique for each data file. A template is saved as a tabular .csv file to reuse and share them.

On this page, you can find information on how to construct a template.


..  _link-target:

Toolkit Standards
====================

The toolkit has standard names for observation types and metadata. Here these standards are presented and described.


.. list-table:: Standard observation types
   :widths: 25 25 15
   :header-rows: 1

   * - Standard name
     - Toolkit description
     - Type
   * - temp
     - temperature
     - numeric
   * - humidity
     - Relative humidity
     - numeric
   * - precip
     - precipitation intensity
     - numeric
   * - precip_sum
     - accumulated precipitation
     - numeric
   * - pressure
     - air pressure (measured)
     - numeric
   * - pressure_at_sea_level
     - air pressure (corrected to sea level)
     - numeric
   * - wind_speed
     - wind speed
     - numeric
   * - wind_gust
     - wind gust
     - numeric
   * - wind_direction
     - wind direction as ° from the north, clock-wise
     - numeric
   * - radiation_temp
     - radiation temperature (black globe observations)
     - numeric


.. list-table:: Standard Metadata
   :widths: 20 25 15
   :header-rows: 1

   * - Standard name
     - Toolkit description
     - Type
   * - name
     - the name of the stations (must be unique for each station)
     - string
   * - lat
     - the latitude of the station
     - numeric
   * - lon
     - the longitude of the station
     - numeric
   * - location
     - location (the city/region of the stations) (OPTIONAL)
     - string
   * - call_name
     - call_name (an informal name of the stations) (OPTIONAL)
     - string
   * - network
     - network (the name of the network the stations belong to) (OPTIONAL)
     - string


In the template, you map your observations and metadata to one of these standards. What is not mapped, will not be used in the toolkit.


Data structures
=======================

To make a template you must be aware of which format your data is in. The toolkit can handle the following data structures:

**long-format**
   Observations are stacked in rows per station. One column represents the station names.

   .. list-table:: long-format example
      :widths: 15 15 15 15
      :header-rows: 1

      * - timestamp
        - 2mT-passive
        - 2m-rel-hum
        - ID
      * - 2022-06-07 13:20:00
        - 16.4
        - 77.3
        - station_A
      * - 2022-06-07 13:30:00
        - 16.7
        - 75.6
        - station_A
      * - 2022-06-07 13:20:00
        - 18.3
        - 68.9
        - station_B
      * - 2022-06-07 13:30:00
        - 18.6
        - 71.9
        - station_B

**Wide-format**
   Columns represent different stations. The data represents one observation type.

   .. list-table:: Wide-format example (temperature)
      :widths: 15 15 15
      :header-rows: 1

      * - timestamp
        - station_A
        - station_B
      * - 2022-06-07 13:20:00
        - 16.4
        - 18.3
      * - 2022-06-07 13:30:00
        - 16.7
        - 18.6

**Single-station-format**
   The same as a long format but without a column indicating the station names.
   Be aware that the toolkit interprets it as **observations coming from one station**.

   .. list-table:: Single-station-format example
      :widths: 15 15 15
      :header-rows: 1

      * - timestamp
        - 2mT-passive
        - 2m-rel-hum
      * - 2022-06-07 13:20:00
        - 16.4
        - 77.3
      * - 2022-06-07 13:30:00
        - 16.7
        - 75.6
      * - 2022-06-07 13:40:00
        - 17.2
        - 77.0
      * - 2022-06-07 13:50:00
        - 17.2
        - 76.9

Metadata structures
=======================
The metadata **must be in a Wide-format**. Here an example

.. list-table:: Metadata example
   :widths: 15 15 15 15
   :header-rows: 1

   * - ID
     - Northening
     - Eastening
     - Networkname
   * - station_A
     - 51.3664
     - 4.67785
     - demo-network
   * - station_B
     - 51.6752
     - 5.1332
     - demo-network


Template creation
=======================

Once you have converted your tabular data files to either long-, wide-, or single-station-format, and saved them as a .csv file, a template can be made.

.. Note::
   If you want to use a metadata file, make sure it is converted to a wide-format and saved as a .csv file.

The fastest and simplest way to make a template is by using the *metobs_toolkit.build_template_prompt()* function.

.. code-block:: python

   import metobs_toolkit

   #create a template
   metobs_toolkit.build_template_prompt()


This function will prompt questions and build a template that matches your data file (and metadata) file.
The *template.csv* file will be stored at a location of your choice.

To use this template, feed the path to the *template.csv* file to the data_template_file (and metadata_template_file)
arguments of the :py:meth:`update_settings()<metobs_toolkit.dataset_settings_updater.Dataset.update_settings>` method.


.. note::
   When the prompt ask's if you need further help, and you type yes, some more questions are prompted.
   Once all information is given to the prompt, it will print out a piece of code that you have to run to load your data into the toolkit.
