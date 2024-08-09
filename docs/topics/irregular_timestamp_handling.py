#!/usr/bin/env python
# coding: utf-8

# # Working with irregular timestamps
#
#
# Some datasets have irregular time frequencies of the observations. These datasets
# come with some extra challenges. Here is some information on how to deal with them.
#
# A common problem that can arise is that most observations are **not present** and
# that **a lot of missing observations** (and gaps) are introduced. This is because
# the toolkit assumes that each station has observations at a constant frequency. So the toolkit expects
# perfectly regular timestamp series. The toolkit will hence ignore observations
# that are not on the frequency, so observations get lost. Also, it looks for observations
# on perfectly regular time intervals, so when a timestamp is not present, it is assumed to be missing.
#
#
# To avoid these problems you can **synchronize** your observations. Synchronizing will
# convert your irregular dataset **to a regular dataset** and an **easy origin** is chosen if possible.
# (The origin is the first timestamp of your dataset.) Converting your dataset to a regular dataset is performed
# by shifting the timestamp of an observation. For example, if a frequency of 5 minutes is assumed and the observation
# has a timestamp at 54 minutes and 47 seconds, the timestamp is shifted to 55 minutes. A certain
# maximal threshold needs to be set to avoid observations being shifted too much. This threshold is
# called the tolerance and it indicates what the **maximal time-translation** error can be for one
# observation timestamp.
#
#
# Synchronizing your observations can be performed with he :py:meth:`sync_observations()<metobs_toolkit.dataset.Dataset.sync_observations>`
# method. As an argument of this function you must provide a tolerance.
#

# ## Example
#
# In this example, we use a small dataset with three stations that have irregular and unsynchronized timestamps. The dataset can be found in [here](https://github.com/vergauwenthomas/MetObs_toolkit/blob/master/tests/test_data/wide_test_data.csv) and the template file can be found [here](https://github.com/vergauwenthomas/MetObs_toolkit/blob/master/tests/test_data/wide_test_template.json).
#
#
#

# In[1]:


import pandas as pd

datafile = "https://raw.githubusercontent.com/vergauwenthomas/MetObs_toolkit/master/tests/test_data/wide_test_data.csv"
templatefile = "https://raw.githubusercontent.com/vergauwenthomas/MetObs_toolkit/master/tests/test_data/wide_test_template.json"


# As an example, here is how the data looks like
df = pd.read_csv(datafile)
df


# It can clearly be seen that timestamps are not synchronized over the stations and that the time resolution is not perfect for the stations. We can fix these timestamps in the toolkit when we import the datafile.
#
# When a datafile is imported, the toolkit will convert the time series to "perfect" timestamps (= equally spaced timestamps). To do this, the toolkit must make an estimate of the **frequency** (=time resolution) for each station. In addition, the toolkit will estimate an **origin** and a **last timestamp** for each station. Once these three parameters (freq, origin, and last timestamp) are computed, a perfect time series is created. At last, the toolkit will map the records in the datafile to these perfect timestamps.
#
#
#
# When the timestamps in the datafile are irregular, we can fix them by specifying tolerances and simplifications for mapping to perfect timestamps when importing the data from file: ``Dataset.import_data_from_file()``

# In[2]:


import metobs_toolkit

dataset = metobs_toolkit.Dataset()
dataset.update_settings(input_data_file=datafile, template_file=templatefile)


dataset.import_data_from_file(
    freq_estimation_method="highest",  # highest or median
    freq_estimation_simplify_tolerance="2min",  # Try to simplify the frequency, the maximum simplification tolerance is 2 minutes
    origin_simplify_tolerance="5min",  # try to simplify the origin, the maximum simplification tolerance is 5 minutes
    timestamp_tolerance="4min",  # The maximum tolerance for mapping records to perfect timestamps is 5 minutes
    templatefile_is_url=True,
)


# The freqency, origin and latest timestamp are stored per station in the `Dataset.metadf` in the `dataset_resolution`, `dt_start` and `dt_end` columns:

# In[3]:


dataset.metadf


# Because of the wide-structured data, the toolkit assumes a 1-minute-frequency. What we can do is to coarsen the time resolution to hourly.
#
# **Note**: to avoid propagation of errors and tolerances, it is best to shift timestamps only once. For this example, this can be done by importing the data without simplifications or tolerances and coarsening it with tolerances.
#
# **Note**: In this wide datafile example, the toolkit interprets the timestamps without a value as gaps, and is the assumption of 1-minute-frequency the correct one. This issue is typically related to wide-structured-datasets.

# In[4]:


import datetime

# import without simlifications and tolerances (to avoid error propagation/cumulation)
dataset.import_data_from_file(
    freq_estimation_method="highest",  # highest or median
    freq_estimation_simplify_tolerance="0min",  # Try to simplify the frequency, the maximum simplification tolerance is 2 minutes
    origin_simplify_tolerance="0min",  # try to simplify the origin, the maximum simplification tolerance is 5 minutes
    timestamp_tolerance="0min",  # The maximum tolerance for mapping records to perfect timestamps is 5 minutes
    templatefile_is_url=True,
)

dataset.sync_records(
    timestamp_shift_tolerance="6min",
    freq_shift_tolerance="0min",
    fixed_origin=None,
    fixed_enddt=None,
    fixed_freq="1h",
    direction="nearest",
)

dataset.make_plot()


# We can see that the timestamps in the Dataset are "perfect", and that the stations are syncronized.

# In[5]:


dataset.get_full_status_df()["temp"]["value"].unstack().transpose()
