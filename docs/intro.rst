
*******************
Introduction
*******************
This package is designed for handling meteorological observations for urban or non traditional networks. It includes tools to clean up and analyse your data.



How to install
=======================

In order to use the package python 3.7 or higher is required.
To install the package one can use pip:

.. code-block:: console

   pip3 install metobs-toolkit

.. note::

   For now this is a development version, so to install you need to specify the latest version explicitly: pip3 install metobs-toolkit==0.0.2ax (where x is the latest version).

.. note::
   To install the package in a notebook, one has to add ! in front of the pip install command.

and import it in python

.. code-block:: python

   import metobs_toolkit

   #Check your version
   metobs_toolkit.__version__


How to use this toolkit
=========================

This toolkit is a python package based on object-oriented programming (OOP). Here you can find a short discription of the classes that are directyly used by the users:


Dataset()
-----------
This class is at the heart of the toolkit and it holds all observation and metadata. The dataset class has three 'containers' to hold observations:

Dataset.df :
    All(*) records will start in the *df-container*. This container contains the observations that we assume to be correct.



    (*): One exception is the observatios with a duplicated timestamp, these will be passed to the outliersdf-container directly.

Dataset.outliersdf :
    When applying quality control, some observations may be labeled as outliers. When an observation is labeled as outlier, it is added to the *outliersdf-container*.
    The record labeled as outlier are still kept inside the df-container but the observation value is removed (set to Nan).



Dataset.missing_obs :
    When importing a datafile, an observation frequency is estimated for eacht station. A missing observation is a record that is not in the observations, but is assumed by the station frequency.
    A missing observation is thus a record, without an observation value. These records are stored in the *missing_obs-container*.

Dataset.gaps:
    When a sequence of (repeating) missing observations is found, a test is performed to check if the length(*) of the series is larger than a threshold (i.e. the gap defenition).
    If the series is larger than the threshold, we interpret it as a *gap* and it is removed from the missing_obs-container.

    (*): Note that the defenition of a gap is based on a number of consecutive repeating missing records! The minimal gapsize is therefore depending on the observational frequency of each station.






Note :
-------
A **record** refers to a unique combination of timestamp, corresponding station and observation type.