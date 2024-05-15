***********************
Mapping to the toolkit
***********************

The MetObs-toolkit uses standard names and formats for your data. To use the toolkit,
your observational data must be converted to the toolkit standards this is referred to as **mapping**.

To specify how the mapping must be done a **template** is used. This template contains
all the information on how to convert your tabular data to the toolkit standards.
Since the structure of data files differs for different networks, this template is
unique for each data file.

A template is saved as a json file (>v0.2.1) to reuse and share. Previous versions of the toolkit used a .csv file type for the template.


On this page, you can find information on data formats and how to construct a template.


Template creation
=======================

Once you have converted your tabular data files to either long-, wide-, or single-station-format, and saved them as a .csv file, a template can be made. See below for more details on these data structures.

.. Note::
   If you want to use a metadata file, make sure it is converted to a long-format (=each station takes up a row) and saved as a .csv file.

The fastest and simplest way to make a template is by using the :py:meth:`metobs_toolkit.build_template_prompt()<metobs_toolkit.data_templates.template_build_prompt.build_template_prompt>` function.

.. code-block:: python

   import metobs_toolkit

   #create a template
   metobs_toolkit.build_template_prompt()


This function will prompt questions and build a template that matches your data file (and metadata) file.
The *template.json* file will be stored at a location of your choice. If you plan to use multiple datafiles, make sure to rename your templates accordingly.

To use this template, feed the path to the *template.csv* file to the data_template_file (and metadata_template_file)
arguments of the :py:meth:`update_settings()<metobs_toolkit.dataset_settings_updater.Dataset.update_settings>` method.

.. code-block:: python

   import metobs_toolkit

   #1. Define the paths to your files:
   data_file = r"... <path_to_your_datafile> ..."
   meta_data_file = r"... <path_to_your_metadatafile> ..."
   template = r"... <path_to_your_templatefile> ..."

   #2. initiate a dataset:
   your_dataset = metobs_toolkit.Dataset()

   #3. Update the paths to your files:
   your_dataset.update_settings(
       input_data_file = data_file,
       input_metadata_file=meta_data_file,
       template_file = template,
       )

   #4. Import your data :
   your_dataset.import_data_from_file()



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
The metadata **must be in a long-format**. Here an example

.. list-table:: Metadata example
   :widths: 15 15 15 15
   :header-rows: 1

   * - ID
     - Northing
     - Easting
     - Networkname
   * - station_A
     - 51.3664
     - 4.67785
     - demo-network
   * - station_B
     - 51.6752
     - 5.1332
     - demo-network

.. note::
   All CSV data files must be in UTF-8 encoding. For most CSV files, this condition is already met. To make sure, in Microsoft Excel (or similar), you can specify to export as **`CSV UTF-8`**.
   If you encounter an error, mentioning a `"/ueff..."` tag in a CSV file, it is solved by converting the CSV to UTF-8.
