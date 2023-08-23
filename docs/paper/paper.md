---
title: 'MetObs-Toolkit - a Python toolkit for using non-traditional meteorological observations'
tags:
  - Python
  - Meteorology
  - Urban climate
  - Observations
authors:
  - name: Thomas Vergauwen
    orcid: 0000-0003-2899-9218
    equal-contrib: false
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Michiel Vieijra
    orcid: 0000-0003-0817-2846
    equal-contrib: False # (This is how you can denote equal contributions between multiple authors)
    affiliation: 1
  - name: Andrei Covaci
    orcid: 0000-0001-5147-2460
    equal-contrib: False # (This is how you can denote equal contributions between multiple authors)
    affiliation: 3
  - name: Amber Jacobs
    equal-contrib: False # (This is how you can denote equal contributions between multiple authors)
    affiliation: 1
  - name: Sara Top
    orcid: 0000-0003-1281-790X
    equal-contrib: False # (This is how you can denote equal contributions between multiple authors)
    affiliation: 1
  - name: Wout Dewettinck
    orcid: 0000-0002-0728-5331
    equal-contrib: False # (This is how you can denote equal contributions between multiple authors)
    affiliation: 1
  - name: Kobe Vandelanotte
    orcid: 0009-0001-1252-7315
    equal-contrib: False # (This is how you can denote equal contributions between multiple authors)
    affiliation: "1, 2"
  - name: Steven Caluwaerts
    orcid: 0000-0001-7456-3891
    equal-contrib: False # (This is how you can denote equal contributions between multiple authors)
    affiliation: "1, 2"
affiliations:
 - name: Royal Meteorological Institute of Belgium, Brussels, Belgium
   index: 1
 - name: University of Ghent, Ghent, Belgium
   index: 2
 - name: University of Brussels (VUB), Brussels, Belgium
   index: 3
 - name: Flemisch Institute for Technological Reseach (VITO), Mol, Belgium
   index: 4
date: 13 August 2023
bibliography: paper.bib
nocite: |
  @*

---


# Summary
In-situ weather station observations are highly important for meteorological and climatological research. The evolution towards more affordable sensor technology and data communication has resulted in the emergence of novel meteorological networks alongside the traditional high-quality measurement networks of meteorological institutions. Examples include urban measurement networks intended to study the impact of cities [@mocca] and networks consisting of devices of weather enthusiasts e.g. [@crowdsourcing_status]. However, exploiting the data of such non-traditional networks comes with significant challenges [@crowdsourcing]. Firstly, sensors and data communication protocols are usually low-cost and this in general results in an increase of measurement errors, biases and data gaps. Secondly, data storage formats and temporal measurement frequencies are not consistent or compatible across different networks. Finally, metadata, such as station location, elevation, and instrument specifications, are not easily accessible or documented.
The MetObs-toolkit is a Python package developed to address these issues and facilitate the use of non-traditional observations. The package provides automated quality control (QC) techniques to identify and flag erroneous observations and includes methods to fill data gaps. Additionally, the package also offers a set of tools for analyzing the data e.g. linkage with popular land use datasets [@worldcover; @lcz_map] is included such that microclimate effects can be investigated by the MetObs toolkit.


# Statement of need
The primary objective of the MetObs toolkit is  to enable scientists to process meteorological observations into datasets ready for analysis. The data cleaning process involves three steps: 

1.  resampling the time resolution if necessary,
2.  identifying erroneous or missing records, and
3.  filling the missing records.
   
Sophisticated software such as TITAN/CrowdQC+  [@titan2020; @CrowdQC] exists for identifying erroneous observations (QC), which is one aspect of cleaning a dataset. These packages are designed specifically for this purpose. Moreover, researchers often face the challenge of coding scripts that can generate analyses, particularly when using geographical datasets such as landcover datasets. Traditionally, this requires the installation of numerous packages, storage of geographical datasets, and GIS manipulations (often manually done with specific GIS software). The toolkit implements one user-friendly framework for creating various plots, generating analysis statistics, and incorporating GIS data through the use of the Google Earth engine.
By using the toolkit, scientists can set up a pipeline to process raw data into analysis in an easy-to-use (and install) manner. Additionally, the developed pipeline can be directly applied to other datasets without any formatting issues.

![A schematic overview of the main MetObs Toolkit functionalities.\label{fig:overview_fig}](overview_fig.png)

# Technical implementation

The MetObs-Toolkit provides a comprehensive framework for scientists to process raw meteorological data for analysis. The process consists of the following steps, visualized in the figure below:

Mapping the raw data to the toolkit standards by creating a template. Once the raw data is imported into the Toolkit Dataset, missing observations are identified and methods to resample and synchronize observations can be used.
Quality control is performed in the form of a series of checks. Advanced quality control methods are available through the implementation of TITAN into the toolkit. The user can choose to keep the outliers or convert them to missing records (which can be filled).
Gap filling is applied by using interpolation methods and/or importing ERA5 reanalysis time series to fill the gaps. The latter is stored as a Toolkit Modeldata, which has a set of methods to directly import the required time series through the use of the Google Earth engine API.
The user obtains a cleaned-up dataset ready for analysis. A set of typical analysis techniques such as filters, aggregation schemes, and landcover correlation estimates are implemented in the Toolkit-Analysis class.

\autoref{fig:overview_fig} gives an overview of the main framework of the MetObs-toolkit, but it is an evolving project that responds to the community's needs and input. As an example, the development of a graphical user interface (GUI) for the toolkit is planned. A GUI would increase the ease of use by enabling to create templates, adjust QC settings and plot data interactively. 

It is to be expected that some novel functionalities will be added in the coming years and the documentation will be updated accordingly. 




# Acknowledgments

The authors would like to thank all participants of the [COST FAIRNESS](https://www.fairness-ca20108.eu/) (CA20108) summer school 2023 in Ghent for their role as Beta testers. This group of scientists, from various European countries, working on urban climate represented the target audience of this package. Their input, ideas, and recommendations were instrumental in improving the MetObs-Toolkit.

No specific funding has been obtained to build the MetObs-Toolkit, but the authors have been supported by different Belgian and Flemish scientific research grants.

FWO: Sara (fellowship 1270723N) and Wout (fellowship 1157523N)
BELSPO: Kobe (B2/223/P1/CORDEX.be II), Thomas (B2/202/P1/CS-MASK), Michiel (B2/212/P2/CLIMPACTH) and Steven (FED-tWIN Prf-2020-018_AURA)
VITO: Ian (VITO PhD fellowship, I will verify if he has a code)

# References
