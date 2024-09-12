#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module contains the Dataset class and all its methods.

A Dataset holds all observations and is at the center of the
MetObs-toolkit.
"""

import os
import sys
import copy
from datetime import timedelta
import pytz
import logging
import pandas as pd
import numpy as np
import pickle


from metobs_toolkit.data_import import import_data_from_csv, import_metadata_from_csv
from metobs_toolkit.printing import print_dataset_info


from metobs_toolkit.settings_files.default_formats_settings import (
    label_def,
)

from metobs_toolkit.writing_files import (
    write_df_to_csv,
    _does_trg_file_exist,
    _remove_file,
    MetobsWritingError,
)

from metobs_toolkit.gap import (
    find_gaps,
    get_station_gaps,
)


from metobs_toolkit.df_helpers import (
    empty_outliers_df,
    metadf_to_gdf,
    # xs_save,
    # concat_save,
    _simplify_time,
)

from metobs_toolkit.obstypes import Obstype as Obstype_class


# dataset extensions
from metobs_toolkit.datasetbase import DatasetBase
from metobs_toolkit.dataset_settings_updater import DatasetSettingsCore
from metobs_toolkit.dataset_visuals import DatasetVisuals
from metobs_toolkit.dataset_qc_handling import DatasetQCCore
from metobs_toolkit.dataset_gap_handling import DatasetGapCore
from metobs_toolkit.dataset_modeldata import DatasetModelData

logger = logging.getLogger(__name__)


# =============================================================================
# Dataset class
# =============================================================================


class Dataset(
    DatasetBase,
    DatasetSettingsCore,
    DatasetVisuals,
    DatasetQCCore,
    DatasetGapCore,
    DatasetModelData,
):
    """Objects holding observations and methods on observations."""

    def __init__(self):
        """Construct all the necessary attributes for Dataset object."""
        logger.info("Initialise dataset")

        DatasetBase.__init__(self)  # holds df, metadf, obstypes and settings
        self._istype = "Dataset"

    def __add__(
        self, other, timestamp_tolerance="0min", freq_simplify_tolerance="0min"
    ):
        """
        Add another Dataset to self.

        Parameters
        ----------
        other : metobs_toolkit.Dataset()
            The dataset to add.
        timestamp_tolerance : Timedelta or str, optional
            The tolerance string or object represents the maximum translation
            in time to apply on a timestamp for conversion to an ideal set of timestamps.
            Ex: '5min' is 5 minutes, '1H', is one hour. A zero-tolerance (thus no
            simplification) can be set by '0min'. The default is '0min'.
        freq_simplify_tolerance : Timedelta or str, optional
            The tolerance string or object represents the maximum translation
            in time to form a simplified frequency estimation.
            Ex: '5min' is 5 minutes, '1H', is one hour. A zero-tolerance (thus no
            simplification) can be set by '0min'. The default is '0min'.

        Returns
        -------
        new : metobs_toolkit.Dataset()
            The combined dataset.

        Warning
        ---------
        Gaps have to be recomputed, so all gaps (and filling information) will
        be lost.

        Warning
        ---------
        In case of duplicates (in records, outliers, metadate) the reference of
        self will be used.

        Warning
        ----------
        The order of QC applied checks, is dubious when different QC methods or
        orders are used between self and other. The impact of this is limited
        to a possible false representation of the individual QC check effectiveness
        when `Dataset.get_qc_stats()` is called.

        Examples
        --------
        Combining two ``Datasets`` together:

        >>> Combined = dataset_A + dataset_B # doctest: +SKIP

        """
        # check if there is data
        self._data_is_required_check()
        other._data_is_required_check()

        logger.info(f"adding {str(other)} to {str(self)}")

        timestamp_tolerance = self._timedelta_arg_check(timestamp_tolerance)
        freq_simplify_tolerance = self._timedelta_arg_check(freq_simplify_tolerance)

        # the toolkit makes a new dataframe, and assumes the df from self and other
        # to be the input data.
        # This means that any progress in gaps is removed an gaps are located
        # again.

        new = Dataset()

        #  ---- df ----
        # combine both long dataframes (without gaps !)
        newdf = pd.concat(
            [
                self.df.drop(self._get_gaps_df_for_stacking().index),
                other.df.drop(other._get_gaps_df_for_stacking().index),
            ]
        )

        if newdf.index.duplicated().any():
            logger.warning(
                "Duplicate records are found between self and other. The records of self are used!"
            )
            # drop duplicated indeces, keep the records from self
            newdf = newdf[~newdf.index.duplicated(keep="first")]

        new._set_df(newdf)

        #  ----- outliers df ---------

        newoutliersdf = pd.concat([self.outliersdf, other.outliersdf])
        if newoutliersdf.index.duplicated().any():
            logger.warning(
                "Duplicate outliers are found between self and other. For duplicates, the outliers of self are used!"
            )
            # drop duplicated indeces, keep the records from self
            newoutliersdf = newoutliersdf[~newoutliersdf.index.duplicated(keep="first")]

        new._set_outliersdf(newoutliersdf)

        # ------- meta df ----------
        newmetadf = pd.concat([self.metadf, other.metadf])
        if newmetadf.index.duplicated().any():
            logger.warning(
                "Duplicate metadata is found between self and other. For duplicates, the outliers of self are used!"
            )
            # drop duplicated indeces, keep the records from self
            newmetadf = newmetadf[~newmetadf.index.duplicated(keep="first")]
        # newmetadf = metadf_to_gdf(newmetadf)
        new._set_metadf(newmetadf)

        # ---- template -----
        # once the mapping is done the template is of no use anymore (obstypes,
        # are updated with descrptions etc). We use the template of self, for
        # sake of having a complete Dataset
        new.template = self.template

        # ----- obstypes ------

        # check if obstypes have the same defenition
        for obstypename in self._get_present_obstypes():
            if obstypename in other._get_present_obstypes():
                if not self.obstypes[obstypename] == other.obstypes[obstypename]:
                    logger.warning(
                        f"The obstype {obstypename} is different between self and other, the defenition of self is used"
                    )

        new_obstypes = other.obstypes.copy()
        new_obstypes.update(self.obstypes)
        new._set_obstypes(new_obstypes)

        # ----- Settings -------
        # overload from self
        new._set_settings(self.settings)

        # ------ QC order --------

        # The QC order of applied qc is (for now) only used to make an
        # estimate on how effective a specific check is --> get_qc_stats.

        if not self._applied_qc.reset_index(drop=True).equals(
            other._applied_qc.reset_index(drop=True)
        ):
            logger.warning(
                "The order and applied QC differs between self and other. QC check effectiveness can not be estimated correctly when calling get_qc_stats()."
            )
        new_applied_qc = pd.concat([self._applied_qc, other._applied_qc]).copy()
        new_applied_qc = new_applied_qc.drop_duplicates(subset=["obstype", "checkname"])
        new._applied_qc = new_applied_qc

        #  ------- Gaps -------------
        # Gaps have to be recaluculated using a frequency assumtion from the
        # combination of self.df and other.df, thus NOT the native frequency if
        # their is a coarsening allied on either of them.
        # find the start, end timestamps and frequency for each station + write it to the metadf
        new._get_timestamps_info(
            freq_estimation_method="highest",  # does not matter on perfect timeseries
            freq_simplify_tolerance=freq_simplify_tolerance,  # Do no chain error oropagation by default
            origin_simplify_tolerance="0min",
        )

        # Convert the records to clean equidistanced records for both the df and outliersdf
        new._construct_equi_spaced_records(
            timestamp_mapping_tolerance=timestamp_tolerance
        )

        # # Find gaps on Import resolution
        gaps = find_gaps(
            df=new.df,
            metadf=new.metadf,
            outliersdf=new.outliersdf,
            obstypes=new.obstypes,
        )
        new._set_gaps(gaps)

        return new

    def show(self, show_all_settings=False, max_disp_n_gaps=5):
        """Show detailed information of the Dataset.

        A function to print out a detailed overview of information about the
        Dataset.

        Parameters
        ----------
        show_all_settings : bool, optional
            If True all the settings are printed out. The default is False.
        max_disp_n_gaps: int, optional
            The maximum number of gaps to display detailed information.

        Returns
        -------
        None.

        See Also
        --------
        get_info: Alias of show()

        Examples
        --------
        Create a ``Dataset`` and fill it with data (and metadata).

        >>> import metobs_toolkit
        >>>
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> dataset.import_data_from_file(
        ...                         input_data_file=metobs_toolkit.demo_datafile,
        ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                         template_file=metobs_toolkit.demo_template,
        ...                         )

        Apply `show` on your Dataset

        >>> # Print out details
        >>> dataset.show()
        --------  General ---------
        ...

        """
        logger.info("Show basic info of dataset.")

        print_dataset_info(self, show_all_settings)

    def get_info(self, show_all_settings=False, max_disp_n_gaps=5):
        """Alias of show().

        A function to print out a detailed overview information about the Dataset.

        Parameters
        ----------
        show_all_settings : bool, optional
            If True all the settings are printed out. The default is False.
        max_disp_n_gaps: int, optional
            The maximum number of gaps to display detailed information of.

        Returns
        -------
        None.

        See Also
        --------
        show:  The same method.

        Examples
        --------
        Create a ``Dataset`` and fill it with data (and metadata).

        >>> import metobs_toolkit
        >>>
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> dataset.import_data_from_file(
        ...                         input_data_file=metobs_toolkit.demo_datafile,
        ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                         template_file=metobs_toolkit.demo_template,
        ...                         )

        Apply `show` on your Dataset

        >>> # Print out details
        >>> dataset.get_info()
        --------  General ---------
        ...

        """
        self.show(show_all_settings, max_disp_n_gaps)

    def save_dataset(
        self, outputfolder=None, filename="saved_dataset.pkl", overwrite=False
    ):
        """Save a Dataset instance to a (pickle) file.

        Parameters
        ----------
        outputfolder : str or None, optional
            The path to the folder to save the file. If None, the outputfolder
            from the Settings is used. The default is None.
        filename : str, optional
            The name of the output file. The default is 'saved_dataset.pkl'.
        overwrite : bool, optional
            If True, the target file will be overwritten if it exist. The
            default is False.

        Returns
        -------
        None.

        See Also
        --------
        import_dataset: Import a dataset from a pickle file.
        write_to_csv: Write data to a csv file.

        Examples
        --------
        Create a ``Dataset`` and fill it with data (and metadata).

        >>> import metobs_toolkit
        >>>
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> dataset.import_data_from_file(
        ...                         input_data_file=metobs_toolkit.demo_datafile,
        ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                         template_file=metobs_toolkit.demo_template,
        ...                         )

        Save it as a pickle file. As a demo, the pickle will be saved in the
        current working directory (`os.getcwd()`)

        >>> # Save dataset to a .pkl file
        >>> dataset.save_dataset(outputfolder=os.getcwd(),
        ...                     filename='your_saved_dataset.pkl',
        ...                     overwrite=True)
        Dataset saved in ...

        """
        # check if outputfolder is known and exists
        if outputfolder is None:
            outputfolder = self.settings.IO["output_folder"]
            assert (
                outputfolder is not None
            ), "No outputfolder is given, and no outputfolder is found in the settings."

        assert os.path.isdir(outputfolder), f"{outputfolder} is not a directory!"

        # check file extension in the filename:
        if filename[-4:] != ".pkl":
            filename += ".pkl"

        full_path = os.path.join(outputfolder, filename)
        if (os.path.isfile(full_path)) & overwrite:
            logger.info(f"The file {full_path} will be overwritten!")
            os.remove(full_path)

        # check if file exists
        assert not os.path.isfile(full_path), f"{full_path} is already a file!"

        with open(full_path, "wb") as outp:
            pickle.dump(self, outp, pickle.HIGHEST_PROTOCOL)

        print(f"Dataset saved in {full_path}")
        logger.info(f"Dataset saved in {full_path}")

    def get_full_status_df(self, return_as_wide=True):
        """Combine all records, outliers, and gaps in one Dataframe

        Records, outliers, and gaps are separately stored in each Dataset. This
        method will combine them into one pandas.DataFrame.

        The full dataframe displays the observation values, a label, and how
        the records is stored in the Dataset (as a good observation, an outlier or gap).

        Parameters
        ----------

        return_as_wide : bool, optional
            If True, the dataset is wide-structured (observationtypes are spread
            over different columns). If False, all records are stacked in
            a long-format. The default is True.

        Returns
        ---------
        combdf : pandas.DataFrame()
            A dataframe containing a continuous time resolution of records, where each
            record is labeled.

        Examples
        --------
        Create a ``Dataset`` and fill it with data (and metadata).

        >>> import metobs_toolkit
        >>>
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> dataset.import_data_from_file(
        ...                         input_data_file=metobs_toolkit.demo_datafile,
        ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                         template_file=metobs_toolkit.demo_template,
        ...                         )

        To get a dataframe with all records, outliers, and gaps, use the
        `get_full_status_df()` method.

        >>> # Combine all into one Dataframe
        >>> combined_df = dataset.get_full_status_df()
        >>> combined_df[['temp', 'humidity']].head() #select only temperature and humidity
        obstype                              temp                              humidity
                                            value label toolkit_representation    value label toolkit_representation
        name      datetime
        vlinder01 2022-09-01 00:00:00+00:00  18.8    ok            observation     65.0    ok            observation
                  2022-09-01 00:05:00+00:00  18.8    ok            observation     65.0    ok            observation
                  2022-09-01 00:10:00+00:00  18.8    ok            observation     65.0    ok            observation
                  2022-09-01 00:15:00+00:00  18.7    ok            observation     65.0    ok            observation
                  2022-09-01 00:20:00+00:00  18.7    ok            observation     65.0    ok            observation

        If you want it in a long structure:

        >>> combined_df = dataset.get_full_status_df(return_as_wide=False)
        >>> combined_df.head()
                                              value label toolkit_representation
        name      obstype  datetime
        vlinder01 humidity 2022-09-01 00:00:00+00:00   65.0    ok            observation
                           2022-09-01 00:05:00+00:00   65.0    ok            observation
                           2022-09-01 00:10:00+00:00   65.0    ok            observation
                           2022-09-01 00:15:00+00:00   65.0    ok            observation
                           2022-09-01 00:20:00+00:00   65.0    ok            observation
        """
        # check if there is data
        self._data_is_required_check()

        # =============================================================================
        # Stack observations and outliers
        # =============================================================================
        df = self.df
        # note: df is a pointer, and adding these colmns will add them
        # also in the self.df
        df["label"] = label_def["goodrecord"]["label"]
        df["toolkit_representation"] = "observation"

        # =============================================================================
        # Stack outliers
        # =============================================================================

        outliersdf = self.outliersdf
        outliersdf["toolkit_representation"] = "outlier"

        combdf = pd.concat([df, outliersdf])  # combine the two

        # Since outliers are present records in the df (as NaN's) we introduce
        # duplicats in the index of combdf. We drop the duplicates and keep,
        # the records comming from outliersdf (=last)

        combdf = combdf[~combdf.index.duplicated(keep="last")]

        # =============================================================================
        # Stack gaps
        # =============================================================================

        gapsdf = (
            self._get_gaps_df_for_stacking()
        )  # get a gapdf in the long (similar as outliersdf) structure
        # map labels to known labels (must have a color def in the settings)
        if not gapsdf.empty:
            gapsdf["label"] = gapsdf["fill_method"].replace(
                {"not filled": label_def["regular_gap"]["label"]}
            )

        gapsdf = gapsdf[["value", "label"]]

        gapsdf["toolkit_representation"] = "gap"

        combdf = pd.concat([combdf, gapsdf])  # combine

        # Since gaps are present records in the df (as NaN's, because of the
        # ideal freq structure in the df) we introduce
        # duplicats in the index of combdf. We drop the duplicates and keep,
        # the records comming from outliersdf (=last)

        combdf = combdf[~combdf.index.duplicated(keep="last")]
        # =============================================================================
        # Formatting the combineddf
        # =============================================================================

        assert (
            not combdf.index.duplicated().any()
        ), "Duplicates found in the combdf --> report bug."

        # for some reason the dtype of the datetime index-level is 'obstype' and
        # thus not a datetimeindex. This must be fixed
        combdf = combdf.reset_index()
        combdf["datetime"] = pd.to_datetime(combdf["datetime"])
        combdf = combdf.set_index(["name", "obstype", "datetime"]).sort_index()

        if return_as_wide:
            combdf = combdf.unstack(level="obstype").reorder_levels(
                order=[1, 0], axis=1
            )

        # pointer issue
        self.df = self.df[["value"]]
        self.outliersdf = self.outliersdf[["value", "label"]]
        return combdf

    def add_new_observationtype(self, Obstype):
        """Add a new observation type to the known observation types.

        The observation can only be added if it is not already present in the
        known observation types. If that is the case that you probably need to
        use use the Dataset.add_new_unit() method.

        Parameters
        ----------
        Obstype: metobs_toolkit.Obstype
            The new Obstype to add.

        Returns
        -------
        None.

        See Also
        --------
        Obstype: The Obstype class.
        add_new_unit : Add a new unit to a known obstype.

        Examples
        --------
        Start by creating a new ``Obstype``

        >>> import metobs_toolkit
        >>> co2_concentration = metobs_toolkit.Obstype(obsname='co2',
        ...                                            std_unit='ppm')
        >>> #add other units to it (if needed)
        >>> co2_concentration.add_unit(unit_name='ppb',
        ...                            conversion=['x / 1000'], #1 ppb = 0.001 ppm
        ...                           )
        >>> #Set a description
        >>> co2_concentration.set_description(desc='The CO2 concentration measured at 2m above surface')

        Create a ``Dataset`` and fill it with data (and metadata).

        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> dataset.import_data_from_file(
        ...                         input_data_file=metobs_toolkit.demo_datafile,
        ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                         template_file=metobs_toolkit.demo_template,
        ...                         )

        Now add the newly created Obtype to the Dataset.

        >>> #Add it to a Dataset
        >>> dataset = metobs_toolkit.Dataset()
        >>> dataset.add_new_observationtype(co2_concentration)

        By inspecting the `Dataset.obstypes` we can see that the CO2 concentration
        is added.

        >>> dataset.obstypes
        {'temp': Obstype instance of temp, 'humidity': Obstype instance of humidity, 'radiation_temp': Obstype instance of radiation_temp, 'pressure': Obstype instance of pressure, 'pressure_at_sea_level': Obstype instance of pressure_at_sea_level, 'precip': Obstype instance of precip, 'precip_sum': Obstype instance of precip_sum, 'wind_speed': Obstype instance of wind_speed, 'wind_gust': Obstype instance of wind_gust, 'wind_direction': Obstype instance of wind_direction, 'co2': Obstype instance of co2}
        """
        # Test if the obstype is of the correct class.
        if not isinstance(Obstype, Obstype_class):
            sys.exit(
                f"{Obstype} is not an instance of metobs_toolkit.obstypes.Obstype."
            )

        # Test if the obsname is already in use
        if Obstype.name in self.obstypes.keys():
            logger.warning(
                f"{Obstype.name} is already a known observation type: {self.obstypes[Obstype.name]}"
            )
            return

        # Update the known obstypes
        logger.info(f"Adding {Obstype} to the list of known observation types.")
        self.obstypes[Obstype.name] = Obstype

    def add_new_unit(self, obstype, new_unit, conversion_expression=[]):
        """Add a new unit to a known observation type.

        Parameters
        ----------
        obstype : str
            The observation type to add the new unit to.
        new_unit : str
            The new unit name.
        conversion_expression : list or str, optional
            The conversion expression to the standard unit of the observation
            type. The expression is a (list of) strings with simple algebraic
            operations, where x represents the value in the new unit, and the
            result is the value in the standard unit. Two examples for
            temperature (with a standard unit in Celsius):

                ["x - 273.15"] #if the new_unit is Kelvin
                ["x-32.0", "x/1.8"] #if the new unit is Fahrenheit

            The default is [].

        Returns
        -------
        None.

        See Also
        --------
        Obstype: The Obstype class.
        add_new_observationtype : Add a new observationtype to the Dataset.

        Examples
        --------

        Create an empty ``Dataset``.

        >>> import metobs_toolkit
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset

        There are default Observationtypes present in all Datasets.

        >>> dataset.obstypes
        {'temp': Obstype instance of temp,
         'humidity': Obstype instance of humidity,
         'radiation_temp': Obstype instance of radiation_temp,
         'pressure': Obstype instance of pressure,
         'pressure_at_sea_level': Obstype instance of pressure_at_sea_level,
         'precip': Obstype instance of precip,
         'precip_sum': Obstype instance of precip_sum,
         'wind_speed': Obstype instance of wind_speed,
         'wind_gust': Obstype instance of wind_gust,
         'wind_direction': Obstype instance of wind_direction}

        We can add a new unit to an existing obstype by using ``add_new_unit()``
        and specifying how it relates to the standard unit.

        You can get the standard unit of an Obsytype by calling the ``get_standard_unit()`` on it.

        >>> dataset.obstypes['temp'].get_standard_unit()
        'Celsius'

        To see all known units for an Obstype, use the `get_all_units()` on it.

        >>> dataset.obstypes['temp'].get_all_units()
        ['Celsius', 'Fahrenheit', 'Kelvin']

        Now we add a new unit

        >>> dataset.add_new_unit(obstype = 'temp',
        ...                      new_unit= 'your_new_unit',
        ...                      conversion_expression = ['x+3', 'x * 2'])
        >>> # The conversion means: 1 [your_new_unit] = (1 + 3) * 2 [°C]
        >>> dataset.obstypes['temp'].get_info()
        temp observation with:
         * standard unit: Celsius
         * data column as Temperatuur in Celsius
         * known units and aliases: {'Celsius': ['celsius', '°C', '°c', 'celcius', 'Celcius'], 'Kelvin': ['K', 'kelvin'], 'Fahrenheit': ['fahrenheit'], 'your_new_unit': []}
         * description: 2mT passive
         * conversions to known units: {'Kelvin': ['x - 273.15'], 'Fahrenheit': ['x-32.0', 'x/1.8'], 'your_new_unit': ['x+3', 'x * 2']}
         * originates from data column: Temperatuur with Celsius as native unit.
        """
        # test if observation is present
        if not obstype in self.obstypes.keys():
            logger.warning(f"{obstype} is not a known obstype! No unit can be added.")
            return

        # check if the unit is already present
        is_present = self.obstypes[obstype].test_if_unit_is_known(new_unit)
        if is_present:
            logger.info(
                f"{new_unit} is already a known unit of {self.obstypes[obstype]}"
            )
            return

        self.obstypes[obstype].add_unit(
            unit_name=new_unit, conversion=conversion_expression
        )

    def show_settings(self):
        """Show detailed information of the stored Settings.

        A function that prints out all the settings, structured per thematic.

        Returns
        -------
        None.

        See Also
        --------
        get_info: Print out basic info of the Dataset.
        Settings: Settings class
        Dataset.settings: The attribute holding the `Settings`

        Examples
        --------

        Create an empty ``Dataset``.

        >>> import metobs_toolkit
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset

        Print out all the current (in this case the default) settings.

        >>> dataset.show_settings()
        All settings:...

        """
        self.settings.show()

    def get_station(self, stationname):
        """Filter out one station of the Dataset.

        Extract a metobs_toolkit.Station object from the dataset by name.

        Parameters
        ----------
        stationname : str
            The name of the station.

        Returns
        -------
        station: metobs_toolkit.Station
            The station object.

        See Also
        --------
        Station: The Station class

        Examples
        --------
        Create a ``Dataset`` and fill it with data (and metadata).

        >>> import metobs_toolkit
        >>>
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> dataset.import_data_from_file(
        ...                         input_data_file=metobs_toolkit.demo_datafile,
        ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                         template_file=metobs_toolkit.demo_template,
        ...                         )

        To see which stations are present in your `Dataset`, you can use the
        metadata or by using the 'name' index of the records:

        >>> #Using the metadata (index):
        >>> dataset.metadf.index
        Index(['vlinder01', 'vlinder02', 'vlinder03', 'vlinder04', 'vlinder05',
               'vlinder06', 'vlinder07', 'vlinder08', 'vlinder09', 'vlinder10',
               'vlinder11', 'vlinder12', 'vlinder13', 'vlinder14', 'vlinder15',
               'vlinder16', 'vlinder17', 'vlinder18', 'vlinder19', 'vlinder20',
               'vlinder21', 'vlinder22', 'vlinder23', 'vlinder24', 'vlinder25',
               'vlinder26', 'vlinder27', 'vlinder28'],
              dtype='object', name='name')

        >>> #Using the records
        >>> dataset.df.index.get_level_values('name').unique()
        Index(['vlinder01', 'vlinder02', 'vlinder03', 'vlinder04', 'vlinder05',
               'vlinder06', 'vlinder07', 'vlinder08', 'vlinder09', 'vlinder10',
               'vlinder11', 'vlinder12', 'vlinder13', 'vlinder14', 'vlinder15',
               'vlinder16', 'vlinder17', 'vlinder18', 'vlinder19', 'vlinder20',
               'vlinder21', 'vlinder22', 'vlinder23', 'vlinder24', 'vlinder25',
               'vlinder26', 'vlinder27', 'vlinder28'],
              dtype='object', name='name')

        Now we can select a station by name

        >>> favorite_station = dataset.get_station('vlinder05')
        >>> print(favorite_station)
        Station instance containing:
         *1 stations
         *['humidity', 'temp', 'wind_direction', 'wind_speed'] observation types present
         *17280 observation records (not Nan's)
         *0 records labeled as outliers
         *0 gaps
         *records range: 2022-09-01 00:00:00+00:00 --> 2022-09-15 23:55:00+00:00 (total duration:  14 days 23:55:00)
         *time zone of the records: UTC
         *Coordinates are available for all stations.
         *Known GEE datasets for: ['lcz', 'altitude', 'worldcover', 'ERA5-land']

        """
        from metobs_toolkit.station import Station  # isn't this tricky !!??

        # check if there is data
        self._data_is_required_check()

        logger.info(f"Extract {stationname} from dataset.")

        if self.df.empty:
            raise MetobsDatasetError(
                f"Cannot get station {stationname} from an empty Dataset."
            )

        # important: make sure all station attributes are of the same time as dataset.
        # so that all methods can be inherited.

        try:
            sta_df = self.df.xs(stationname, level="name", drop_level=False)
            sta_metadf = self.metadf.loc[stationname].to_frame().transpose()
            sta_metadf.index.name = "name"
        except KeyError:
            logger.warning(f"{stationname} not found in the dataset.")
            return None

        try:
            sta_outliers = self.outliersdf.xs(
                stationname, level="name", drop_level=False
            )
        except KeyError:
            sta_outliers = empty_outliers_df()

        sta_gaps = get_station_gaps(self.gaps, stationname)

        return Station(
            name=stationname,
            df=sta_df,
            outliersdf=sta_outliers,
            gaps=sta_gaps,
            gee_datasets=list(self.gee_datasets.values()),
            metadf=sta_metadf,
            obstypes=self.obstypes,
            template=self.template,
            settings=self.settings,
            _applied_qc=self._applied_qc.copy(),
        )

    def write_to_csv(
        self,
        data_file,
        metadata_file=None,
        overwrite=False,
        all_outliers_as_nan=False,
        kwargs_to_csv={"index": False},
    ):
        """Write Dataset to a csv file.

        Write all present records to a file. This is done by combining the
        good records, outliers and gaps into one dataframe. An extra
        'label-column', for each observationtype is added.

        The metadata can be written to a separate file.


        Parameters
        ----------
        data_file : str or None
            The path to a csv location to write the data to. If the path does
            not have an ".csv" extension, it will be added. If None, no data is
            written.
        metadata_file : str or None, optional
            The path to a csv location to write the metadata to. If the path does
            not have an ".csv" extension, it will be added. If None, no metadata is
            written.
        overwrite : bool, optional
            If True, the target files will be written even if they already exist.
            The default is False.
        all_outliers_as_nan : bool, optional
            If True, all records flagged as outliers are represented by Nan.
            The default is False.
        kwargs_to_csv : dict, optional
            Kwargs that are passed to the pandas.DataFrame.to_csv() method.
            The default is {'index': False}.

        Returns
        -------
        None.

        See Also
        --------
        save_dataset: Save a dataset as a pickle file.


        Examples
        --------
        Create a ``Dataset`` and fill it with data (and metadata).

        >>> import metobs_toolkit
        >>>
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> dataset.import_data_from_file(
        ...                         input_data_file=metobs_toolkit.demo_datafile,
        ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                         template_file=metobs_toolkit.demo_template,
        ...                         )

        Now we will create a file (.csv) with all records and a file with the metadata.

        >>> # create paths to target files
        >>> import os
        >>> target_data_file = os.path.join(os.getcwd(), 'your_data.csv')
        >>> target_metadatafile = os.path.join(os.getcwd(), 'your_metadata.csv')
        >>>
        >>> dataset.write_to_csv(data_file=target_data_file,
        ...                      metadata_file=target_metadatafile,
        ...                      overwrite=True)
        write to file: /home/.../your_data.csv
        write to file: /home/.../your_metadata.csv
        """

        logger.info("Writing the dataset to a csv file")

        # check if there is data
        self._data_is_required_check()

        # check path to target files
        if data_file is None:
            # do not write data
            write_data = False
        else:
            write_data = True
            # make sure it has a .csv filetype extension
            if not str(data_file).endswith(".csv"):
                # add .csv
                logger.warning(
                    f'{data_file} has no ".csv" extension, this will be added to the path.'
                )
                data_file = f"{data_file}.csv"
            # check if path exist
            if _does_trg_file_exist(str(data_file)):
                if not overwrite:
                    raise MetobsWritingError(
                        f"{data_file} already exists. Use overwrite=True, or change the path."
                    )
                else:
                    _remove_file(str(data_file))

        # check path to target files
        if metadata_file is None:
            # do not write data
            write_metadata = False
        else:
            write_metadata = True
            if not str(metadata_file).endswith(".csv"):
                # add .csv
                logger.warning(
                    f'{metadata_file} has no ".csv" extension, this will be added to the path.'
                )
                metadata_file = f"{metadata_file}.csv"
            # check if path exist
            if _does_trg_file_exist(str(metadata_file)):
                if not overwrite:
                    raise MetobsWritingError(
                        f"{metadata_file} already exists. Use overwrite=True, or change the path."
                    )
                else:
                    _remove_file(str(metadata_file))

        if (not write_data) & (not write_metadata):
            raise MetobsDatasetError(
                f"Cannot write data to a file since no target datafile or target metadatafile is given."
            )

        # =============================================================================
        # Write data
        # =============================================================================

        if write_data:
            combdf = self.get_full_status_df()

            present_obs = combdf.columns.get_level_values(0).unique()
            # Too a single columnindex
            combdf.columns = combdf.columns.map("_".join)

            # drop the toolkit representation columns
            combdf = combdf.drop(
                columns=[
                    col
                    for col in combdf.columns
                    if col.endswith("_toolkit_representation")
                ]
            )

            for ob in present_obs:
                if all_outliers_as_nan:
                    combdf.loc[combdf[f"{ob}_label"] != "ok", f"{ob}_value"] = np.nan
            # write to file
            write_df_to_csv(
                df=combdf,
                trgfile=str(data_file),
                to_csv_kwargs=kwargs_to_csv,
            )

        if write_metadata:
            metadf = self.metadf

            # add present obstype info to the metadf
            for obs in self.df.index.get_level_values("obstype").unique():
                obstype = self.obstypes[obs]
                metadf[f"{obs}_unit"] = obstype.get_standard_unit()
                metadf[f"{obs}_description"] = obstype.get_description()

            write_df_to_csv(
                df=metadf,
                trgfile=str(metadata_file),
                to_csv_kwargs=kwargs_to_csv,
            )
        return

    def coarsen_time_resolution(
        self,
        origin=None,
        freq="60min",
        direction="nearest",
        timestamp_shift_tolerance="4min",
        # limit=1,
    ):
        """Resample the observations to coarser timeresolution.

        This method is used to convert the time resolution of the Dataset. This
        will affect the records (.df), the outliers (.outliersdf) and gaps (.gaps).



        The assumed dataset resolution (stored in the metadf attribute) will be
        updated.

        Parameters
        ----------
        origin : datetime.datetime, optional
            Define the origin (first timestamp) for the observations. The origin
            is timezone naive and is assumed to have the same timezone as the
            obervations. If None, the earliest occurring timestamp is used as
            origin. The default is None.
        freq : DateOffset, Timedelta or str, optional
            The target frequency to coarsen all records to.
            Ex: '15min' is 15 minutes, '1h', is one hour. The default is '60min'.
        direction : 'backward', 'forward', or 'nearest'
            Whether to search for prior, subsequent, or closest matches for
            mapping to ideal timestamps. The default is 'nearest'.
        timestamp_shift_tolerance : Timedelta or str
            The tolerance string or object representing the maximum translation
            (in time) to map a timestamp to a target timestamp.
            Ex: '5min' is 5 minutes. The default is '4min'.


        Returns
        -------
        None.

        Note
        ------
        This method is actually a call to `Dataset.sync_records()` with a fixed
        frequency and a zero-tolerance on frequency shift.

        Warning
        ---------
        Since the gaps depend on the record's frequency and origin, all gaps are
        removed and re-located. All progress in gap(filling) will be lost.

        Note
        -------
        It is technically possible to increase the time resolution. This will
        however not result in an information increase, more gaps are
        created instead.

        See Also
        --------
        sync_records: Syncronize records in time.

        Examples
        --------
        Create a ``Dataset`` and fill it with data (and metadata).

        >>> import metobs_toolkit
        >>>
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> dataset.import_data_from_file(
        ...                         input_data_file=metobs_toolkit.demo_datafile,
        ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                         template_file=metobs_toolkit.demo_template,
        ...                         )

        When importing a data from file, a frequency assumption is made for
        each station (separatly) in the Dataset. This is stored in the `metadf`

        >>> dataset.metadf['dataset_resolution']
        name
        vlinder01   0 days 00:05:00
        vlinder02   0 days 00:05:00
        vlinder03   0 days 00:05:00
        vlinder04   0 days 00:05:00
        vlinder05   0 days 00:05:00
                          ...
        vlinder24   0 days 00:05:00
        vlinder25   0 days 00:05:00
        vlinder26   0 days 00:05:00
        vlinder27   0 days 00:05:00
        vlinder28   0 days 00:05:00
        Name: dataset_resolution, Length: 28, dtype: timedelta64[ns]

        We can coarsen the time resolution to 15 minutes with a maximum
        allowd timestamp shift of 4 minutes.

        >>> dataset.coarsen_time_resolution(freq='15min',
        ...                                 timestamp_shift_tolerance='4min')
        >>> dataset.metadf['dataset_resolution']
        name
        vlinder01   0 days 00:15:00
        vlinder02   0 days 00:15:00
        vlinder03   0 days 00:15:00
        vlinder04   0 days 00:15:00
        vlinder05   0 days 00:15:00
                          ...
        vlinder24   0 days 00:15:00
        vlinder25   0 days 00:15:00
        vlinder26   0 days 00:15:00
        vlinder27   0 days 00:15:00
        vlinder28   0 days 00:15:00
        Name: dataset_resolution, Length: 28, dtype: timedelta64[ns]

        """
        # check if there is data
        self._data_is_required_check()
        # coarsening/resamplin is the same as sync with a fixed frequency.
        self.sync_records(
            timestamp_shift_tolerance=timestamp_shift_tolerance,
            freq_shift_tolerance="0min",
            fixed_origin=origin,
            fixed_enddt=None,
            fixed_freq=freq,  # fixed freq
            direction=direction,
        )

    def sync_records(
        self,
        timestamp_shift_tolerance="2min",
        freq_shift_tolerance="1min",
        fixed_origin=None,
        fixed_enddt=None,
        fixed_freq=None,
        direction="nearest",
    ):
        """Synchronize observation records in time.

        Synchronization is shifting timestamps by litle so that the same
        timestamps are found across multiple stations. The maximum allowed shift
        of a timestamp is set by timestamp_shift_tolerance.

        Synchronized observations must have the same frequency. This method will
        therefore try to alter the frequency a bit, so that it is in sync with other
        observations (from other stations).

        The timestamp-related columns in the Dataset.metadf are updated, and
        gaps are reconstructed! All current info stored in the gaps will be lost.

        In practice, the synchronisation is applied at the start of your workflow
        on datasets with irregular or unsynchronized timestamps.

        Parameters
        ----------
        timestamp_shift_tolerance : str or Timedelta, optional
            The tolerance string or object represents the maximum translation
            in time for a timestamp. The default is "2min".
        freq_shift_tolerance : str or Timedelta, optional,
            The tolerance string or object representing the maximum translation
            in time for a frequency. The default is "1min".
        fixed_origin : datetime.datetime, optional
            Define the origin (first timestamp) for the observations. This is
            applied on all the stations. If the origin is timezone naive, it is
            assumed to have the same timezone as the observations. If None, the
            (shifted) earliest occurring timestamp is used as the origin. The
            default is None.
        fixed_enddt : datetime.datetime, optional
            Define the latest timestamp for the obervations. This is
            applied on all the stations. If the fixed_enddt is timezone naive,
            it is assumed to have the same timezone as the obervations. If
            None, the (shifted) latest occuring timestamp is used as enddt
            The default is None.
        fixed_freq : str or Timedelta, optional,
            The target frequency applied on all records. If None,
            (a simplified version of the) present frequency is used. The default is None.
        direction : 'backward', 'forward', or 'nearest'
            Whether to search for prior, subsequent, or closest matches for
            mapping to ideal timestamps. The default is 'nearest'.


        Returns
        -------
        None.

        Note
        --------
        Keep in mind that this method will overwrite the df, outliersdf and gaps.

        See Also
        --------
        coarsen_time_resolution: To coarsen the time resolution to a fixed freq.

        Examples
        --------
        Create a ``Dataset`` and fill it with data (and metadata).

        >>> import metobs_toolkit
        >>>
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> dataset.import_data_from_file(
        ...                         input_data_file=metobs_toolkit.demo_datafile,
        ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                         template_file=metobs_toolkit.demo_template,
        ...                         )

        When importing a data from file, a frequency, origin and enddt assumption
        is made for each station (separatly) in the Dataset. This is stored in the `metadf`.

        >>> dataset.metadf[['dataset_resolution', 'dt_start', 'dt_end']].head()
                  dataset_resolution                  dt_start                    dt_end
        name
        vlinder01    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder02    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder03    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder04    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder05    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00

        When you have irregular/unsynchronized timestamp, you can sync them

        >>> dataset.sync_records(
        ...     timestamp_shift_tolerance="7min",
        ...     freq_shift_tolerance="4min")

        """
        # check if there is data
        self._data_is_required_check()

        # checking arguments
        timestamp_shift_tolerance = self._timedelta_arg_check(timestamp_shift_tolerance)
        freq_shift_tolerance = self._timedelta_arg_check(freq_shift_tolerance)
        fixed_origin = self._datetime_arg_check(fixed_origin)
        fixed_enddt = self._datetime_arg_check(fixed_enddt)
        fixed_freq = self._timedelta_arg_check(fixed_freq)

        logger.info(f"Syncronizing records")
        # 1. find common freqencies

        # TODO: This method is technically not the most suitable method, but it is
        # straight forward and will work on many cases.

        # The ideal method is to create groups base on the dataset resolution (unsimplified),
        # grouping by target simplified frequencies within the tolerance. It may sound
        # easy, but it is not straight forward to implement IMO
        if fixed_freq is not None:
            simple_freqs = pd.Series(
                data=pd.Timedelta(fixed_freq), index=self.metadf.index
            )
        else:
            simple_freqs = self.metadf["dataset_resolution"].apply(
                lambda x: _simplify_time(
                    time=x,
                    max_simplyfi_error=freq_shift_tolerance,
                    zero_protection=True,
                )
            )

        logger.debug(f"Simplified target resolutions: {simple_freqs}")
        # 2. find common origins
        if fixed_origin is not None:
            simple_origins = pd.Series(data=fixed_origin, index=self.metadf.index)
        else:
            simple_origins = self.metadf["dt_start"].apply(
                lambda x: _simplify_time(
                    time=x,
                    max_simplyfi_error=timestamp_shift_tolerance,
                )
            )
        logger.debug(f"Simplified target origins: {simple_origins}")

        # 3. find common end timestamps
        if fixed_enddt is not None:
            simple_ends = pd.Series(
                data=pd.Timedelta(fixed_enddt), index=self.metadf.index
            )
        else:
            simple_ends = self.metadf["dt_end"].apply(
                lambda x: _simplify_time(
                    time=x,
                    max_simplyfi_error=timestamp_shift_tolerance,
                )
            )
        logger.debug(f"Simplified target end timestamps: {simple_ends}")
        # 4. Update the metadata, and restructure the records to these targets
        self.metadf["dataset_resolution"] = simple_freqs
        self.metadf["dt_start"] = simple_origins
        self.metadf["dt_end"] = simple_ends

        # Convert the records to clean equidistanced records for both the df and outliersdf
        # note: The outliers are taken care of as well
        self._construct_equi_spaced_records(
            timestamp_mapping_tolerance=timestamp_shift_tolerance, direction=direction
        )

        # better save than sorry: update metadf automatically
        self._get_timestamps_info(
            freq_estimation_method="highest",  # does not matter on perfect frequency
            freq_simplify_tolerance="0min",  # no simplification
            origin_simplify_tolerance="0min",  # no simplification
        )

        # # Find gaps on Import resolution
        gaps = find_gaps(
            df=self.df,
            metadf=self.metadf,
            outliersdf=self.outliersdf,
            obstypes=self.obstypes,
        )
        self._set_gaps(gaps)

    def import_only_metadata_from_file(
        self,
        input_metadata_file=None,
        template_file=None,
        kwargs_metadata_read={},
        templatefile_is_url=False,
    ):
        """Import metadata from a CSV file.

        This method will read in metadata (typical station coordinates), from
        a CSV file. The Dataset is updated with this metadata.

        This method is used when you have no observational data but want to
        work with the metadata.


        Parameters
        ----------
        input_metadata_file : string, optional
            Path to the input metadata file. If None, the input metadata path
            in the settings is used. The default is None
        template_file : string, optional
            Path to the template (JSON) file to be used on the observations
            and metadata. If None, the template path in the settings is used.
            If provided, the input_data_file must be provided as well. The
            default is None.
         kwargs_metadata_read : dict, optional
             Keyword arguments are collected in a dictionary to pass to the
             pandas.read_csv() function on the metadata file. The default is {}.
         templatefile_is_url : bool, optional
             If the path to the template file, is a url to an online template file,
             set templatefile_is_url to True. If False, the template_file is
             interpreted as a path.

        Returns
        -------
        None.

        Warning
        ----------
        Since no observation data is set, a lot of `Dataset` methods will not
        be applicable to a "metadata-only" `Dataset`. If you call one of these
        methods, an Exception will be raised.


        Warning
        ---------
        All CSV data files must be in *UTF-8 encoding*. For most CSV files,
        this condition has already been met. To make sure, in Microsoft Excel (or
        similar), you can specify to export as **`CSV UTF-8`**. If you
        encounter an error, mentioning a `"/ueff..."` tag in a CSV file, it is
        often solved by converting the CSV to UTF-8.

        Note
        -----
        If a template file is provided the mapping information of
        observational data is ignored. Thus you can use the same template
        for "metadata-only" analysis.

        See Also
        --------
        update_settings: Update the (file paths) settings of a Dataset.
        import_data_from_file: Import observational data from a CSV file.
        import_dataset: Import a dataset from a pkl file.

        Examples
        --------

        .. plot::
            :context: close-figs

            >>> import metobs_toolkit
            >>>
            >>> # Import data into a Dataset
            >>> dataset = metobs_toolkit.Dataset()
            >>> dataset.update_file_paths(
            ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
            ...                         template_file=metobs_toolkit.demo_template,
            ...                         )
            >>> dataset.import_only_metadata_from_file()

            Now you have imported only the metadata of the demo dataset.

            >>> print(dataset)
            Instance of a Dataset (metadata-only).
                 *28 stations in the metadata
                 *The following columns are present in the metadf: ['geometry', 'lat', 'lon', 'school']
                 *Coordinates are available for all stations.
                 *Known GEE datasets for: ['lcz', 'altitude', 'worldcover', 'ERA5-land']

            With no observational data, outliers or gaps.

            >>> dataset.df.empty
            True
            >>> dataset.outliersdf.empty
            True
            >>> dataset.gaps

            Although not all `Dataset` methods are applicable on a "metadata-only",
            dataset, extracting (timeseries) data from GEE is possible (if coordinates
            are provided).

            For example, we extract the LCZ at the stations and ERA5 temperatures.

            >>> lcz_at_stations = dataset.get_lcz()
            >>> lcz_at_stations
                                       lcz
            name
            vlinder01   Low plants (LCZ D)
            vlinder02        Large lowrise
            vlinder03         Open midrise
            vlinder04       Sparsely built
            vlinder05        Water (LCZ G)
            ...                        ...
            vlinder24  Dense Trees (LCZ A)
            vlinder25        Water (LCZ G)
            vlinder26         Open midrise
            vlinder27      Compact midrise
            vlinder28         Open lowrise
            <BLANKLINE>
            [28 rows x 1 columns]

            For extracting ERA5 timeseries we must specify a start and end timestamp.

            >>> import pandas as pd
            >>> start = pd.Timestamp('2017-01-01T12')
            >>> end = pd.Timestamp('2017-01-02T16')
            >>>
            >>> era5 = dataset.get_modeldata(Model=dataset.gee_datasets['ERA5-land'],
            ...                              obstypes=['temp'],
            ...                              startdt=start,
            ...                              enddt=end)

            >>> era5.make_plot()
            <Axes: title={'center': 'ERA5-land'}, ylabel='temp (Celsius)...

        """

        # get template file
        if template_file is not None:
            trg_templ_file = str(template_file)
        elif self.settings.templatefile is not None:
            trg_templ_file = str(self.settings.templatefile)
        else:
            raise MetobsDatasetError(f"No templatefile is specified.")

        # get metadatafile
        if input_metadata_file is not None:
            trg_meta_file = str(template_file)
        elif self.settings.IO["input_metadata_file"] is not None:
            trg_meta_file = str(self.settings.IO["input_metadata_file"])
        else:
            raise MetobsDatasetError(f"No metadata file is specified.")

        # Read template
        logger.info(f"Reading the templatefile at {trg_templ_file}")
        self.template.read_template_from_file(
            jsonpath=trg_templ_file, templatefile_is_url=templatefile_is_url
        )

        meta_df = import_metadata_from_csv(
            input_file=trg_meta_file,
            template=self.template,
            kwargs_metadata_read=kwargs_metadata_read,
        )

        # handle the name

        if "name" in meta_df.columns:
            pass  # normal case
        elif meta_df.shape[0] > 1:
            # no names and multiple stations (= rows) are found --> create dummies
            logger.warning(
                "No column specified that represents the name of the stations, dummy names are created."
            )
            meta_df["name"] = [f"dummy_{i}" for i in range(meta_df.shape[0])]
        else:
            # No name found, and only one station found --> use the single station name
            meta_df["name"] = self.template._get_single_station_default_name()

        meta_df = meta_df.set_index("name")

        # Construct columns that are required
        if "lat" not in meta_df.columns:
            meta_df["lat"] = np.nan
        if "lon" not in meta_df.columns:
            meta_df["lon"] = np.nan

        # Convert to geopandas dataframe
        meta_df = metadf_to_gdf(meta_df)

        # set the attribute
        self._set_metadf(metadf=meta_df)

    def import_data_from_file(
        self,
        input_data_file=None,
        input_metadata_file=None,
        template_file=None,
        freq_estimation_method="highest",
        freq_estimation_simplify_tolerance="2min",
        origin_simplify_tolerance="5min",
        timestamp_tolerance="4min",
        kwargs_data_read={},
        kwargs_metadata_read={},
        templatefile_is_url=False,
    ):
        """Read observations from a csv file and fill the Dataset.


        The input data (and metadata) are interpreted by using a template
        (JSON file).

        In order to locate gaps, an ideal set of timestamps is expected. This
        set of timestamps is computed for each station separatly by:

         * Assuming a constant frequency. This frequency is estimated by using
           a freq_estimation_method. If multiple observationtypes are present,
           the assumed frequency is the highest estimated frequency among
           the differnt observationtypes. To simplify the estimated frequency a
           freq_estimation_simplify_error can be specified.
         * A start timestamp (origin) is found for each station. If multiple
           observationtypes are present,
           the start timestamp is the first timestamp among
           the different observationtypes. The start
           timestamp can be simplified by specifying an origin_simplify_tolerance.
         * The last timestamp is found for each station by taking the timestamp
           which is closest and smaller than the latest timestamp found of a station,
           and is an element of the ideal set of timestamps.

        Each present observation record is linked to a timestamp of this ideal set,
        by using a 'nearest' merge. If the timedifference is smaller than the
        timestamp_tolerance, the ideal timestamp is used. Otherwise, the timestamp
        will be interpreted as a (part of a) gap.


        The Dataset attributes are set and the following checks are executed:
                * Duplicate check
                * Invalid input check
                * Find gaps


        Parameters
        ----------
        input_data_file : string, optional
            Path to the input data file with observations. If None, the input
            data path in the settings is used. If provided, the
            template_file must be provided as well. The default is None.
        input_metadata_file : string, optional
            Path to the input metadata file. If None, the input metadata path
            in the settings is used. The default is None
        template_file : string, optional
            Path to the template (JSON) file to be used on the observations
            and metadata. If None, the template path in the settings is used.
            If provided, the input_data_file must be provided as well. The
            default is None.
        freq_estimation_method : 'highest' or 'median', optional
            Select which method to use for the frequency estimation. If
            'highest', the highest appearing frequency is used. If 'median', the
            median of the appearing frequencies is used. The default is 'highest'.
        freq_estimation_simplify_tolerance : Timedelta or str, optional
            The tolerance string or object represents the maximum translation
            in time to form a simplified frequency estimation.
            Ex: '5min' is 5 minutes, '1H', is one hour. A zero-tolerance (thus no
            simplification) can be set by '0min'. The default is '2min' (2 minutes).
        origin_simplify_tolerance : Timedelta or str, optional
            The tolerance string or object representing the maximum translation
            in time to apply on the start timestamp to create a simplified timestamp.
            Ex: '5min' is 5 minutes, '1H', is one hour. A zero-tolerance (thus no
            simplification) can be set by '0min'. The default is '5min' (5 minutes).
        timestamp_tolerance : Timedelta or str, optional
            The tolerance string or object represents the maximum translation
            in time to apply on a timestamp for conversion to an ideal set of timestamps.
            Ex: '5min' is 5 minutes, '1H', is one hour. A zero-tolerance (thus no
            simplification) can be set by '0min'. The default is '4min' (4 minutes).
        kwargs_data_read : dict, optional
            Keyword arguments are collected in a dictionary to pass to the
            pandas.read_csv() function on the data file. The default is {}.
        kwargs_metadata_read : dict, optional
            Keyword arguments are collected in a dictionary to pass to the
            pandas.read_csv() function on the metadata file. The default is {}.
        templatefile_is_url : bool, optional
            If the path to the template file, is a url to an online template file,
            set templatefile_is_url to True. If False, the template_file is
            interpreted as a path.

        Returns
        -------
        None.

        Note
        --------
        In practice, the default arguments will be sufficient for most applications.

        Note
        --------
        If options are present in the template, these will have priority over the arguments of this function.

        Warning
        ---------
        All CSV data files must be in *UTF-8 encoding*. For most CSV files,
        this condition has already been met. To make sure, in Microsoft Excel (or
        similar), you can specify to export as **`CSV UTF-8`**. If you
        encounter an error, mentioning a `"/ueff..."` tag in a CSV file, it is
        often solved by converting the CSV to UTF-8.

        See Also
        --------
        update_settings: Update the (file paths) settings of a Dataset.
        import_only_metadata_from_file: Import metadata without observational data.
        import_dataset: Import a dataset from a pkl file.

        Examples
        --------

        >>> import metobs_toolkit
        >>>
        >>> # Import data into a Dataset
        >>> dataset = metobs_toolkit.Dataset()
        >>> dataset.update_file_paths(
        ...                         input_data_file=metobs_toolkit.demo_datafile,
        ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                         template_file=metobs_toolkit.demo_template,
        ...                         )
        >>> dataset.import_data_from_file()

        """
        # Special argschecks
        freq_estimation_simplify_tolerance = self._timedelta_arg_check(
            freq_estimation_simplify_tolerance
        )
        origin_simplify_tolerance = self._timedelta_arg_check(origin_simplify_tolerance)
        timestamp_tolerance = self._timedelta_arg_check(timestamp_tolerance)

        # check input file args
        if (input_data_file is not None) & (template_file is None):
            raise MetobsDatasetError(
                f"A input_data_file is provided, but the template_file is missing."
            )
        if (input_data_file is None) & (template_file is not None):
            raise MetobsDatasetError(
                f"A template_file is provided, but the input_data_file is missing."
            )

        # Update paths to the input files, if given.
        if input_data_file is not None:
            self.update_file_paths(
                input_data_file=input_data_file,
                template_file=template_file,
                input_metadata_file=input_metadata_file,
            )

        logger.info(f'Importing data from file: {self.settings.IO["input_data_file"]}')

        assert self.settings.templatefile is not None, "No templatefile is specified."

        # Read template
        logger.info(f"Reading the templatefile")
        self.template.read_template_from_file(
            jsonpath=self.settings.templatefile, templatefile_is_url=templatefile_is_url
        )

        # Read observations into pandas dataframe
        logger.info(f"Reading the observations from file")
        df = import_data_from_csv(
            input_file=self.settings.IO["input_data_file"],
            template=self.template,
            known_obstypes=list(self.obstypes.keys()),
            kwargs_data_read=kwargs_data_read,
        )

        # drop Nat datetimes if present
        df = df.loc[pd.notnull(df.index)]

        logger.debug(
            f'Data from {self.settings.IO["input_data_file"]} \
                      imported to dataframe {df.head()}.'
        )

        if self.settings.IO["input_metadata_file"] is None:
            logger.warning(
                "No metadata file is defined,\
                    no metadata attributes can be set!"
            )

            use_metadata = False

        else:
            logger.info(
                f'Importing metadata from file: {self.settings.IO["input_metadata_file"]}'
            )
            use_metadata = True
            meta_df = import_metadata_from_csv(
                input_file=self.settings.IO["input_metadata_file"],
                template=self.template,
                kwargs_metadata_read=kwargs_metadata_read,
            )
            # in dataset of one station
            if self.template._is_data_single_station():
                # logger.warning("No station names find in the observations!")

                # If there is ONE name in the metadf, than we use that name for
                # the df, else we use the default name
                if ("name" in meta_df.columns) & (meta_df.shape[0] == 1):
                    name = meta_df["name"].iloc[0]
                    df["name"] = name
                    logger.warning(
                        f"One stationname found in the metadata: {name}, this name is used for the data."
                    )
                else:
                    df["name"] = str(self.settings.app["default_name"])
                    # for later merging, we add the name column with the default
                    # also in the metadf
                    meta_df["name"] = str(self.settings.app["default_name"])
                    logger.warning(
                        f'Assume the dataset is for ONE station with the \
                        default name: {self.settings.app["default_name"]}.'
                    )

            # merge additional metadata to observations
            logger.debug(f"Head of data file, before merge: {df.head()}")
            logger.debug(f"Head of metadata file, before merge: {meta_df.head()}")

            meta_cols = [
                colname for colname in meta_df.columns if not colname.startswith("_")
            ]
            additional_meta_cols = list(set(meta_cols).difference(df.columns))

            if bool(additional_meta_cols):
                logger.debug(
                    f"Merging metadata ({additional_meta_cols}) to dataset data by name."
                )
                additional_meta_cols.append("name")  # merging on name
                # merge deletes datetime index somehow? so add it back.
                df_index = df.index
                df = df.merge(
                    right=meta_df[additional_meta_cols], how="left", on="name"
                )
                df.index = df_index

        # update dataset object

        # Remove stations whith only one observation (no freq estimation)
        station_counts = df["name"].value_counts()
        issue_station = station_counts[station_counts < 2].index.to_list()
        if bool(issue_station):

            logger.warning(
                f"These stations will be removed because of only having one record: {issue_station}"
            )
            df = df[~df["name"].isin(issue_station)]

        # convert dataframe to multiindex (datetime - name)
        df = df.set_index(["name", df.index])

        # Sort by name and then by datetime (to avoid negative freq)
        df = df.sort_index(level=["name", "datetime"])

        self._construct_dataset(
            df=df,
            freq_estimation_method=freq_estimation_method,
            freq_estimation_simplify_tolerance=freq_estimation_simplify_tolerance,
            origin_simplify_tolerance=origin_simplify_tolerance,
            timestamp_tolerance=timestamp_tolerance,
            use_metadata=use_metadata,
        )


def import_dataset(folder_path, filename="saved_dataset.pkl"):
    """Import a Dataset instance from a (pickle) file.

    Parameters
    ----------
    folder_path : str
        The path to the directory where the pickle file is stored.
    filename : str, optional
        The name of the output file. The default is 'saved_dataset.pkl'.

    Returns
    -------
    metobs_toolkit.Dataset
        The Dataset instance.

    See Also
    --------
    Dataset.save_dataset: Save a Dataset as a pickle file.

    Examples
    --------

    As an example, we will import a Dataset that is stored as a pkl file in the
    current working directory (`os.getcwd()`).

    >>> import os
    >>> import metobs_toolkit
    >>>
    >>> dataset=metobs_toolkit.import_dataset(folder_path=os.getcwd(),
    ...                                      filename='your_saved_dataset.pkl')  # doctest: +SKIP
    >>> dataset # doctest: +SKIP
    Instance of Dataset at ...
    """
    # check if folder_path is known and exists
    assert os.path.isdir(folder_path), f"{folder_path} is not a directory!"

    full_path = os.path.join(folder_path, filename)

    # check if file exists
    assert os.path.isfile(full_path), f"{full_path} does not exist."

    with open(full_path, "rb") as inp:
        dataset = pickle.load(inp)

    # convert metadf to a geodataframe (if coordinates are available)
    dataset.metadf = metadf_to_gdf(dataset.metadf)

    return dataset


# =============================================================================
# Exceptions
# =============================================================================


class MetobsDatasetError(Exception):
    """Exception raised for errors in the template."""

    pass


# =============================================================================
# Docstring test
# =============================================================================
if __name__ == "__main__":
    from metobs_toolkit.doctest_fmt import setup_and_run_doctest

    setup_and_run_doctest()
