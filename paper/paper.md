---
title: 'MetObs-Toolkit - a Python toolkit for using non-traditional meteorological observations'
tags:
  - Python
  - Meteorology
  - Urban climate
  - Observations
authors:
  - name: Adrian M. Price-Whelan
    orcid: 0000-0000-0000-0000
    equal-contrib: true
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Author Without ORCID
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: 2
  - name: Author with no affiliation
    corresponding: true # (This is how to denote the corresponding author)
    affiliation: 3
affiliations:
 - name: Lyman Spitzer, Jr. Fellow, Princeton University, USA
   index: 1
 - name: Institution Name, Country
   index: 2
 - name: Independent Researcher, Country
   index: 3
date: 13 August 2017
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Astrophysical Journal <- The name of the AAS journal.
---


# Summary
Meteorological networks consist of automated weather stations (AWS) that monitor atmospheric variables in a specific region. These networks are increasingly established by research institutes and universities as sensor technology becomes more affordable. The main objectives of these networks are to examine the systematic influences of land use and land cover on atmospheric processes. Urban meteorological networks, for instance, investigate the urbanization effects on the local climate.

However, these networks often face challenges due to the lack of standardization and quality control. Some of these challenges are:

* AWS are usually low-cost and may have measurement errors or biases.
* Siting errors and biases
* Data storage formats are not consistent or compatible across different networks.
* Metadata, such as station location, elevation, and instrument specifications, are not easily accessible or documented.

The MetObs-toolkit is a Python package that aims to address these issues and facilitate the use of non-traditional observations for further research. The package provides automated quality control (AQC) techniques to identify and flag erroneous observations and to fill in missing observations or gaps. The package also offers a set of tools for analyzing the climatological impact of land use and land cover based on the cleaned observations.


# Statement of need
Non-traditional observations are known to be affected by errors or missing values. Some researchers may choose to discard such observations from their datasets to ensure data quality. However, this may not be feasible or desirable when the data span long periods or when the network size is small. Moreover, it may not be obvious whether an observation is valid or erroneous. Therefore, quality control methods are essential for non-traditional observations. 

There are sophisticated software tools that perform various quality checks based on statistical outlier detection (TITAN paper, CrowdQC+ paper). However, these tools require a specific format and structure of the data, which may entail some data preprocessing by the user. When the time series of different stations are not synchronized, this becomes quite challenging and thus automated quality control is often neglected.

Datasets from non-traditional sources often suffer from missing observations, which are mainly caused by technical failures of weather stations or their communication systems. These missing observations create gaps in the data that reduce the availability and quality of the data for analysis. Since these networks usually have a limited number of stations, and there is little overlap in the landcover representation among different stations, filling these gaps is essential for obtaining reliable results. Therefore, a method for gap-filling is needed.

Furthermore, processing geospatial datasets can be challenging, especially when automating GIS processes. Scientists who lack GIS expertise often download entire datasets, use specific GIS software to visualize the data, and manually annotate landcover information for each station. This approach is time-consuming, error-prone, and limited to point extraction of categorical data. Hence, there is a demand for some GIS tools, within the same Python framework, that automates GIS extractions. By integrating the Google Earth Engine API into this toolkit, we can achieve this and also avoid downloading large geospatial datasets by using Google Cloud.

A common challenge for researchers is to create scripts that can generate analyses that are applicable to different networks. This limits the possibilities of code-sharing or comparison studies among networks. The metobs-toolkit dataset structure provides a standard format that enables the reuse of analysis scripts across networks.


# Technical implementation

The MetObs-Toolkit provides a comprehensive framework for processing raw datasets and preparing them for analysis. This involves applying AQC and filling in missing observations. The user has the option to fill in observations marked as outliers as well. 



isert figure


The MetObs-Toolkit consists of three main classes: Dataset, Analysis, and Modeldata. The Dataset class contains the data, which are categorized into 'good observations', outliers, missing observations and gaps, as well as metadata. The Modeldata class stores time series from models at the station locations, while the Analysis class stores observations for analysis and methods to perform various analyses. For improving user experience an instance of the Settings class, with default values, is passed and stored on the creation of these classes.

# Features

At the time of writing, the most prominent features of the MetObs-Toolkit are:

* All objects have user-friendly methods for plotting and obtaining summary information with get_info().
* The time frequency for each station is estimated using a specific method.. The user can then resample the time resolutions, synchronize stations or convert frequencies to a more convenient frequency within a given tolerance.
* A list of quality control checks is implemented including more advanced spatial QC methods by making use of the Titanlib package.
* Techniques for filling in missing observations and gaps. The user can decide whether to treat observations labeled as outliers as missing observations.
* The user can access the Google Earth Engine, which hosts a large collection of geospatial datasets. These interactions are typically used for
   * Extracting time series from model data (ERA5 recommendations), which can be used for (debiased) gap filling.
   * Extracting geospatial metadata at the station locations. Commonly used categories are altitudes, Local climate zones (LCZ), and landcover fractions in circular buffers.
* The Analysis provides the following functionalities:
   * data filtering based on query or time-derivatives criteria
   * data aggregation into time-derivatives groups for cycle analysis (e.g., diurnal, annual, seasonal, etc.) and visualization
   * landcover effect investigation by using correlation and correlation variation across time-derivatives cycles.


# Futur lookout

The MetObs-Toolkit is an evolving project that responds to the community's needs and the developers' enhancements. The standardized format of data in the toolkit facilitates the application of these enhancements to various datasets. We present some of the features and requests that have been proposed by the community and developers to improve the toolkit:
* Developing a graphical user interface (GUI) for the toolkit. A GUI would enable users to create templates, adjust QC settings and plot data interactively. The MetObs-GUI is currently under construction.
* Using satellite products (on GEE) to filter data by overcast conditions.
* Implementation of spatiotemporal methods for Gap filling like (link Italiaanse paper)
* Interaction with the COST FAIRNESS knowledge sharing platform (KSP). This platform aims to standardize data and metadata for meteorological networks and to make them FAIR. Users can then directly download datasets from the KSP through an API without creating a template for each network.
* Incorporating open data storage platforms (WOW, Netatmo, etc.) to allow users to include crowdsourced observations in their analysis.
* Extending the analysis with spatial methods (i.g. Kringing interpolation) 

# Acknowledgments

The authors gratefully acknowledge the contributions of the Atmospheric Physics group of Ghent University for their valuable feedback and support in developing this package. Their efforts in creating well-documented exercises and demos using the MetObs-Toolkit are greatly appreciated.

The authors also express their gratitude to Michiel Vieijra and Andrei Covaci for their code contributions and insights. Moreover, the authors thank Amber Jacobs for her work on gap-filling techniques as part of her Master's thesis.

Furthermore, the authors would like to thank all participants of the COST FAIRNESS summer school 2023 in Ghent for their role as Beta testers. This group of scientists, from various European countries, working on urban climate represented the target audience of this package. Their input, ideas, and recommendations were instrumental in improving the MetObs-Toolkit.

Finally, the authors extend their special thanks to Prof. Steven Caluwaerts for his guidance and trust in developing this toolkit.


