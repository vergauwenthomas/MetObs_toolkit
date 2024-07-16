=====================================
Using irregular timestamp datasets
=====================================

Some datasets have irregular time frequencies of the observations. These datasets
come with some extra challenges. Here is some information on how to deal with them.

A common problem that can arise is that most observations are **not present** and
that **a lot of missing observations** (and gaps) are introduced. This is because
the toolkit assumes that each station has observations at a constant frequency. So the toolkit expects
perfectly regular timestamp series. The toolkit will hence ignore observations
that are not on the frequency, so observations get lost. Also, it looks for observations
on perfectly regular time intervals, so when a timestamp is not present, it is assumed to be missing.


To avoid these problems you can **synchronize** your observations. Synchronizing will
convert your irregular dataset **to a regular dataset** and an **easy origin** is chosen if possible.
(The origin is the first timestamp of your dataset.) Converting your dataset to a regular dataset is performed
by shifting the timestamp of an observation. For example, if a frequency of 5 minutes is assumed and the observation
has a timestamp at 54 minutes and 47 seconds, the timestamp is shifted to 55 minutes. A certain
maximal threshold needs to be set to avoid observations being shifted too much. This threshold is
called the tolerance and it indicates what the **maximal time-translation** error can be for one
observation timestamp.


Synchronizing your observations can be performed with he :py:meth:`sync_observations()<metobs_toolkit.dataset.Dataset.sync_observations>`
method. As an argument of this function you must provide a tolerance.

Example
---------
Let's take an example dataset with Netatmo(*) data. These data are known for having irregular
timestamps. On average the time resolution is 5 minutes. In the data file,
we can see that there are 4320 observational records. However, when we import it
into the toolkit, only 87 observational records remain:

(*) `Netatmo <https://www.netatmo.com/nl-be/smart-weather-station>`_ is a commercial company that sells automatic weather stations
for personal use.


.. code-block:: python

   #code illustration

   #initialize dataset
   your_dataset = metobs_toolkit.Dataset()

   #specify paths
   dataset.update_settings(
                           input_data_file=' .. path to netatmo data ..',
                           data_template_file=' .. template file .. ',
                           )
   #import the data
   dataset.import_data_from_file()

   print(dataset)

   Dataset instance containing:
        *1 stations
        *['temp', 'humidity'] observation types
        *87 observation records
        *0 records labeled as outliers
        *85 gaps
        *0 missing observations
        *records range: 2021-02-27 08:56:22+00:00 --> 2021-03-13 18:45:56+00:00 (total duration:  14 days 09:49:34)

The toolkit assumes a certain value for the frequency for each station. We can find this in the .metadf attribute:

.. code-block:: python

   print(dataset.metadf['dataset_resolution'])

   name
   netatmo_station   0 days 00:05:00
   Name: dataset_resolution, dtype: timedelta64[ns]



We can synchronize the dataset using this code example:

.. code-block:: python

   # Code illustration

   # Initialize dataset
   your_dataset = metobs_toolkit.Dataset()

   # Specify paths
   dataset.update_settings(
                           input_data_file=' .. path to netatmo data ..',
                           data_template_file=' .. template file .. ',
                           )
   # Import the data
   dataset.import_data_from_file(**testdata[use_dataset]['kwargs'])

   # Syncronize the data with a tolerance of 3 minutes
   dataset.sync_observations(tolerance='3min')

   print(dataset)

   Dataset instance containing:
        *1 stations
        *['temp', 'humidity'] observation types
        *4059 observation records
        *938 records labeled as outliers
        *0 gaps
        *92 missing observations
        *records range: 2021-02-27 08:55:00+00:00 --> 2021-03-13 18:45:00+00:00 (total duration:  14 days 09:50:00)


   # Note: the frequency is not changed
   print(dataset.metadf['dataset_resolution'])

   name
   netatmo_station   0 days 00:05:00
   Name: dataset_resolution, dtype: timedelta64[ns]


The :py:meth:`sync_observations()<metobs_toolkit.dataset.Dataset.sync_observations>` method can also
be used to synchronize the time series of multiple stations. In that case, the method works by trying to find stations with similar
resolutions, finding an origin that works for all stations in this group, and creating a regular time series.
