#!/usr/bin/env python
# coding: utf-8

# # Logging
#
# The toolkit generates logging messages, that can warn and help the users. By default logs are not written or stored. If you want to read the logs, you can do this by
# * forwarding the logs to a file --> this file will be updated with log messages (live)
# * forwarding the logs to stout (i.g. your coding envrionment) --> this will print out the logs in your coding environment (live)
#
# There are convenient functions for both cases. Specify the logginglevel to specify how detailed logs you want.

# ## Stream Logs
#
# To stream logs, call the `add_StreamHandler()` function.

# In[1]:


import metobs_toolkit

# Forward the logs to the coding environment
metobs_toolkit.add_StreamHandler(setlvl="INFO")

# Import a dataset to generate logs
dataset = metobs_toolkit.Dataset()
data_file = metobs_toolkit.demo_datafile
metadata_file = metobs_toolkit.demo_metadatafile
template_file = metobs_toolkit.demo_template
dataset.update_settings(
    input_data_file=data_file,
    input_metadata_file=metadata_file,
    template_file=template_file,
)
dataset.import_data_from_file()


# ## File logs
#
# To write logs to a file, call the `add_FileHandler()` function, and sepcify the path.

# In[4]:


metobs_toolkit.add_FileHandler(
    trglogfile="...",  # path to the logfile (can be nonexisting)
    setlvl="DEBUG",
    clearlog=True,
)


# In[ ]:
