#!/usr/bin/env python
# coding: utf-8

# # Mapping to the toolkit
#
# The MetObs-toolkit uses standard names and formats for your data. To use the toolkit,
# your observational data must be converted to the toolkit standards this is referred to as **mapping**.
#
# To specify how the mapping must be done a **template** is used. This template contains
# all the information on how to convert your tabular data to the toolkit standards.
# Since the structure of data files differs for different networks, this template is
# unique for each data file. A template is saved as a tabular .json file to reuse and share them.
#
# On this page, you can find information on how to construct a template.

# # Toolkit Standards
#
# The toolkit has standard names for observation types and metadata. Here these standards are presented and described.

# In[9]:


# This codeblock is for illustration, it has no practical use.
from metobs_toolkit.miscellaneous import _tlk_print_standard_obstypes

_tlk_print_standard_obstypes()


# ## Data Structures
#
# To make a template you must be aware of which format your data is in. The toolkit can handle the following data structures:

# ### Long-format
# Observations are stacked in rows per station. One column represents the station names.

# | Timestamp   | 2m Temperature | 2m Humidity | ID |
# | -------- | ------- |  ------- | ------- |
# | 2022-06-07 13:20:00  | 16.4 |  77.3 | Station_A |
# | 2022-06-07 13:30:00  | 16.7 |  75.6 | Station_A |
# | 2022-06-07 13:20:00  | 18.3 |  68.9 | Station_B |
# | 2022-06-07 13:30:00  | 18.6 |  71.9 | Station_B |
#

# ### Single-statio-format
# The same as a long format but without a column indicating the station names. Be aware that the toolkit interprets it as observations coming from one station.

# | Timestamp   | 2m Temperature | 2m Humidity |
# | -------- | ------- |  ------- |
# | 2022-06-07 13:20:00  | 16.4 |  77.3 |
# | 2022-06-07 13:30:00  | 16.7 |  75.6 |

# ### Wide-format
# Columns represent different stations. The data represents one observation type.

# | Timestamp   | Station_A | Station_B |
# | -------- | ------- |  ------- |
# | 2022-06-07 13:20:00  | 16.4 |  18.3 |
# | 2022-06-07 13:30:00  | 16.7 |  18.6 |

# ## Template creation
#
# Once you have converted your tabular data files to either long-, wide-, or single-station-format, and saved them as a .csv file, a template can be made.

#
# <div class="alert alert-block alert-info">
# <b>Note:</b> If you want to use a metadata file, make sure it is converted to a Wide-format and saved as a .csv file.
#
# </div>

# The fastest and simplest way to make a template is by using the `metobs_toolkit.build_template_prompt()` function.

# ```python
# import metobs_toolkit
#
# #create a template
# metobs_toolkit.build_template_prompt()
# ```

# <div class="alert alert-block alert-info">
# <b>Note:</b> When the prompt asks if you need further help, and you type yes, some more questions are prompted. Once all information is given to the prompt, it will print out a piece of code that you have to run to load your data into the toolkit.
# </div>

# <div class="alert alert-block alert-warning">
# <b>Warning:</b> All CSV data files must be in UTF-8 encoding. For most CSV files, this condition is already met. To make sure, in Microsoft Excel (or similar), you can specify to export as `CSV UTF-8`.
#    If you encounter an error, mentioning a `"/ueff..."` tag in a CSV file, it is solved by converting the CSV to UTF-8.
# </div>

# This function will prompt questions and build a template that matches your data file (and metadata) file. The *template.json* file will be stored at a location of your choice.
#
# To use this template, add its file path to the arguments of the `update_settings()` method.

# In[10]:


import metobs_toolkit

your_dataset = metobs_toolkit.Dataset()  # initiate an empty dataset
your_dataset.update_settings(
    input_data_file=metobs_toolkit.demo_datafile,  # Path to your data (csv) file
    input_metadata_file=metobs_toolkit.demo_metadatafile,  # Path to your metadata (csv) file
    template_file=metobs_toolkit.demo_template,
)  # Path to your template (json) file.


# The template (file) is read when calling the `Dataset.import_data_from_file()` method, and converted to a `metobs_toolkit.Template` which is accesible for each dataset.

# In[11]:


your_dataset.import_data_from_file()  # will read the data, metadata and template.

your_dataset.template


# An overview of the template can be printed using the `show()` on the `Template` instance:

# In[12]:


your_dataset.template.show()
