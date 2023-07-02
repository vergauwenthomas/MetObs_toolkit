#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module contains the Modeldata class and all its methods.

A Modeldata holds all timeseries coming from a model and methods to use them.
"""

import pandas as pd
import logging

from metobs_toolkit.df_helpers import init_multiindexdf, conv_tz_multiidxdf, xs_save, multiindexdf_datetime_subsetting

from metobs_toolkit.landcover_functions import connect_to_gee, gee_extract_timeseries

from metobs_toolkit.plotting_functions import model_timeseries_plot, timeseries_plot

from metobs_toolkit.convertors import convert_to_toolkit_units, standard_tlk_units

from metobs_toolkit.settings import Settings

logger = logging.getLogger(__name__)

# =============================================================================
# Class Model data (collection of external model data)
# =============================================================================


class Modeldata:
    def __init__(self, modelname):
        self.df = init_multiindexdf
        self.modelname = modelname

        self._settings = Settings()
        self.mapinfo = self._settings.gee["gee_dataset_info"]
        self.mapinfo.update(self._settings.alaro['info'])

        self._df_units = {} #the units of the data stored in the df
        self.df_tz = 'UTC'# the timezone of the datetimes stored in the df

        self._is_alaro25=False


    def __str__(self):
        if self.df.empty:
            return f'Empty Modeldata instance.'
        n_stations = self.df.index.get_level_values('name').unique().shape[0]
        obstypes = self.df.columns.to_list()
        startdt = self.df.index.get_level_values('datetime').min()
        enddt = self.df.index.get_level_values('datetime').max()

        return (f"Modeldata instance containing: \n \
    * Modelname: {self.modelname} \n \
    * {n_stations} timeseries \n \
    * The following obstypes are available: {obstypes} \n \
    * Data has these units: {self._df_units} \n \
    * From {startdt} --> {enddt} (with tz={self.df_tz}) \n \n (Data is stored in the .df attribute)")

    def __repr__(self):
        return self.__str__()


    def add_gee_dataset(self, mapname, gee_location, obstype, bandname, units,
                        scale, time_res='1H', is_image=False, is_numeric=True, credentials=''):
        """
        Method to add a new gee dataset to the available gee datasets.


        Parameters
        ----------
        mapname : str
            Mapname of choice for the GEE dataset to add.
        gee_location : str
            Location of the gee dataset (like "ECMWF/ERA5_LAND/HOURLY" for ERA5).
        obstype : str
            The observation type the band corresponds to.
        bandname : str
            Name of the dataset band as stored on the GEE.
        units : str
            The units of the band.
        scale : int
            The scale to represent the dataset in. (This is a GEE concept that
            is similar to the resolution in meters).
        time_res : timedelta string, optional
            Time reoslution of the dataset, if is_image == False. The default is '1H'.
        is_image : bool, optional
            If True, the dataset is a ee.Image, else it is assumed to be an
            ee.ImageCollection. The default is False.
        is_numeric : bool, optional
            If True, the bandvalues are interpreted as numerical values rather
            than categorical.. The default is True.
        credentials : str, optional
            Extra credentials of the dataset. The default is ''.

        Returns
        -------
        None.

        Note
        -------
        To list all available gee dataset, use the .list_gee_dataset() method.

        Note
        -------
        Currently no unit conversion is perfomed automatically other than K -->
        Celcius. This will be implemented in the futur.

        """
        # check if mapname exists
        if mapname in self.mapinfo.keys():
            logger.warning(f'{mapname} is found in the list of known gee datasets: {list(self.mapinfo.keys())}, choose a different mapname.')
            return



        if is_numeric:
            val_typ='numeric'
        else:
            val_typ = 'categorical'

        if is_image:
            im_bool=True
        else:
            im_bool=False



        new_info = {
            mapname: {
                'location': f'{gee_location}',
                'usage': 'user defined addition',
                'band_of_use' :
                    {f'{obstype}':
                         {'name': f'{bandname}',
                          'units': f'{units}'}
                         },
                        'value_type': val_typ,
                        'dynamical': not bool(is_image),
                        'scale': int(scale),
                        'is_image':bool(is_image),
                        'is_imagecollection': not bool(is_image),
                        'credentials' : f'{credentials}',
                        }
                }

        if not is_image:
            new_info[mapname]['time_res'] = f'{time_res}'

        self.mapinfo.update(new_info)
        logger.info(f'{mapname} is added to the list of available gee dataset with: {new_info}')
        return




    def add_band_to_gee_dataset(self, bandname, obstype, units, overwrite=False):
        """
        A method to add a new band to the current gee dataset (by .modelname attribute).


        Parameters
        ----------

        bandname : str
            Name of the dataset band as stored on the GEE.
        obstype : str
            The observation type the band corresponds to.
        units : str
            The units of the band.
        overwrite : bool, optional
            If True, verwrite the exising bandname when the corresponding
            obstype is already mapped to a bandname. The default is False.

        Returns
        -------
        None.

        Note
        -------
        To list all available gee dataset, use the .list_gee_dataset() method.

        Note
        -------
        Currently no unit conversion is perfomed automatically other than K -->
        Celcius. This will be implemented in the futur.

        """
        mapname = self.modelname

        # check if mapname exists
        if mapname not in self.mapinfo.keys():
            logger.warning(f'{mapname} is not found in the list of known gee datasets: {list(self.mapinfo.keys())}')
            return

        if self.mapinfo[mapname]['is_image']:
            logger.warning(f'{mapname} is found as a Image. No bandnames can be added to it.')
            return


        # check if obstype is already mapped if multiple bands exist
        if not isinstance(self.mapinfo[mapname]['band_of_use'], str):
            if obstype in self.mapinfo[mapname]['band_of_use'].keys():
                if not overwrite:
                    logger.warning(f'{obstype} already mapped to a bandname for dataset: {mapname}.')
                    return


        # update the dict
        new_info = {obstype: {'name' : bandname,
                              'units' : units}}
        self.mapinfo[mapname]['band_of_use'].update(new_info)

        logger.info(f'{new_info} is added to the {mapname} bands of use.')
        return



    def list_gee_datasets(self):
        """
        Print out all the available gee datasets

        Returns
        -------
        None.

        """
        print(f'The following datasets are found: ')
        for geename, info in self.mapinfo.items():
            print('\n --------------------------------')
            print(f'{geename} : \n')
            print(f'{info}')


    def _conv_to_timezone(self, tzstr):
        """
        Convert the timezone of the datetime index of the df attribute.

        Parameters
        ----------
        tzstr : str
            TImezonstring from the pytz module.

        Returns
        -------
        None.

        """
        # get tzstr by datetimindex.tz.zone


        df=self.df
        df['datetime_utc'] = df.index.get_level_values('datetime').tz_convert(tzstr)
        df = df.reset_index()
        df = df.drop(columns=['datetime'])
        df = df.rename(columns={'datetime_utc': 'datetime'})
        df = df.set_index(['name', 'datetime'])
        self.df = df
        self.df_tz = tzstr

    def convert_units_to_tlk(self, obstype, target_unit_name='Celsius',
                      conv_expr=None):
        """
        Method to convert the model data of one observation to the standard
        units as used by the metobs_toolkit.

        If No standard unit is present, you can give a conversion expression.

        The data attributes will be updated.

        Parameters
        ----------
        obstype : str
            Observation type to convert to standard units.
        target_unit_name : str, optional
            Target unit name to convert to. The default is 'Celsius'.
        conv_expr : str, optional
            If the target_unit_name is not a default, you can add the
            conversion expression here (i.g. "x - 273.15"). The default is None.

        Returns
        -------
        None.

        Note
        -------
        All possible mathematical operations for the conv_expr are [+, -, *, /].
        x represent the value in the current units. So "x - 273.15" represents
        the conversion from Kelvin to Celcius.

        """



        # chech if data is available
        if self.df.empty:
            logger.warning('No data to set units for.')
            return
        if not obstype in self.df.columns:
            logger.warning('{obstype} not found as observationtype in the Modeldata.')
            return

        if not conv_expr is None:
            new_unit_def = {target_unit_name:{
                                self._df_units[obstype] : f'{conv_expr}'}}
        else:
            new_unit_def={}

        new_data, new_unit = convert_to_toolkit_units(data=self.df[obstype],
                                                    data_unit=self._df_units[obstype],
                                                    new_units=new_unit_def)

        logger.info(f'{obstype} are converted from {self._df_units[obstype]} --> {new_unit}.')

        self.df[obstype] = new_data
        self._df_units[obstype] = new_unit



    def get_gee_dataset_data(self, mapname, metadf,
                             startdt_utc, enddt_utc, obstype='temp',
                             target_unit_name='new unit', conv_expr=None):

        """
        Extract timeseries of a gee dataset. The extraction can only be done
        if the gee dataset bandname (and units) corresponding to the obstype
        is known.

        The units are converted to the toolkit standard units.


        Parameters
        ----------
        mapname : str
            Mapname of choice of the GEE dataset to extract data from.
        metadf : pandas.DataFrame
            A dataframe with a 'name' index and  'lat', 'lon' columns.
            Timeseries are extracted for these locations.
        startdt_utc : datetime.datetime
            Start datetime of the timeseries in UTC.
        enddt_utc : datetime.datetime
            Last datetime of the timeseries in UTC.
        obstype : str, optional
            Toolkit observation type to extract data from. There should be a
            bandname mapped to this obstype for the gee map. The default is
            'temp'.
        target_unit_name : str, optional
            If there is on standard unit for your obstype, or if you do not
            want to convert to the standard unit, you can specify the name of
            the unit you whant to convert to. This will only be used when a
            conversion expression is provided using the conv_expr argument. The
            default is 'new unit'.
        conv_expr : str, optional
            If the target_unit_name is not a default, you can add the
            conversion expression here (i.g. "x - 273.15"). The default is None.


        Returns
        -------
        None.

        Note
        ------
        When extracting large amounts of data, the timeseries data will be
        writen to a file and saved on your google drive. In this case, you need
        to provide the Modeldata with the data using the .set_model_from_csv()
        method.

        Note
        -------
        All possible mathematical operations for the conv_expr are [+, -, *, /].
        x represent the value in the current units. So "x - 273.15" represents
        the conversion from Kelvin to Celcius.

        """

        # ====================================================================
        # Test input
        # ====================================================================
        if metadf.empty:
            logger.warning(f'The metadf is empty!')
            return

        # Subset metadf to stations with coordinates
        no_coord_meta = metadf[metadf[['lat','lon']].isna().any(axis=1)]
        if not no_coord_meta.empty:
            logger.warning(f'Following stations do not have coordinates, and thus no modeldata extraction is possible: {no_coord_meta.index.to_list()}')
            metadf = metadf[~metadf[['lat','lon']].isna().any(axis=1)]

        # is mapinfo available
        if mapname not in self.mapinfo.keys():
            logger.warning(f'{mapname} is not a known gee dataset.')
            return

        geeinfo = self.mapinfo[mapname]

        # does dataset contain time evolution
        if not geeinfo['dynamical']:
            logger.warning(f'{mapname} is a static dataset, this method does not work on static datasets')
            return

        # is obstype mapped?
        if not obstype in geeinfo['band_of_use'].keys():
            logger.warning(f'{obstype} is not yet mapped to a bandname in the {mapname} dataset.')
            return

        # can observation be converted to standaard units?
        try:
            convert_to_toolkit_units(data = [10,20,30],
                                     data_unit = geeinfo['band_of_use'][obstype]['units'])
        except:
            logger.warning(f"The {geeinfo['band_of_use'][obstype]['units']} cannot be converted to standard toolkit units: ")
            # this prints more details
            convert_to_toolkit_units(data = [10,20,30],
                                     data_unit = geeinfo['band_of_use'][obstype]['units'])


        # ====================================================================
        # GEE api extraction
        # ====================================================================

        # Connect to Gee
        connect_to_gee()

        # Get data using GEE
        df = gee_extract_timeseries(
                                    metadf=metadf,
                                    mapinfo=geeinfo,
                                    startdt=startdt_utc,
                                    enddt=enddt_utc,
                                    obstype=obstype,
                                    latcolname="lat",
                                    loncolname="lon",
                                    )



        if not df.empty:
            self._df_units[obstype] = geeinfo['band_of_use'][obstype]['units']
            if conv_expr is None:
                # use standard units
                self.convert_units_to_tlk(obstype=obstype,
                                      target_unit_name=standard_tlk_units[obstype],
                                      )
            else:
                self.convert_units_to_tlk(obstype=obstype,
                                      target_unit_name=target_unit_name,
                                      conv_expr=conv_expr
                                      )

            self.df_tz='UTC'
        else:
            self._data_stored_at_drive=True

        self.df = df
        self.modelname = mapname



    def get_ERA5_data(self, metadf, startdt_utc, enddt_utc, obstype='temp'):
        """
        Extract timeseries of the ERA5_hourly dataset.

        The units are converted to the toolkit standard units.

        (This method is a specific ERA5_hourly wrapper on the
         get_gee_dataset_data() method)

        Parameters
        ----------

        metadf : pandas.DataFrame
            A dataframe with a 'name' index and  'lat', 'lon' columns.
            Timeseries are extracted for these locations.
        startdt_utc : datetime.datetime
            Start datetime of the timeseries in UTC.
        enddt_utc : datetime.datetime
            Last datetime of the timeseries in UTC.
        obstype : str, optional
            Toolkit observation type to extract data from. There should be a
            bandname mapped to this obstype for the gee map. The default is
            'temp'.


        Returns
        -------
        None.

        Note
        ------
        When extracting large amounts of data, the timeseries data will be
        writen to a file and saved on your google drive. In this case, you need
        to provide the Modeldata with the data using the .set_model_from_csv()
        method.

        """

        self.get_gee_dataset_data(mapname='ERA5_hourly',
                                  metadf=metadf,
                                  startdt_utc=startdt_utc,
                                  enddt_utc=enddt_utc,
                                  obstype=obstype)


    def set_alaro_25_model_from_csv(self, csvpath):
        """
        (This is for the participants of the Cost FAIRNESS Summerschool in Ghent.)

        This method will import the data from the ALARO model, that was send
        to you.


        Parameters
        ----------
        csvpath : str
            Path to the datafile with ALARO timeseries. (This file was send
            to you by email).

        Returns
        -------
        None.

        """

        # update name
        if self.modelname != 'ALARO_2.5':
            logger.info(f'Converting modelname: {self.modelname} --> ALARO_2.5')
            self.modelname ='ALARO_2.5'


        info = self.mapinfo['ALARO_2.5']


        # read in file
        df = pd.read_csv(csvpath, sep=",")

        # Subset to columns in the template
        keep_cols = [val['name'] for val in info['band_of_use'].values()]
        keep_cols.append(info['other_mapping']['datetime']['name'])
        keep_cols.append(info['other_mapping']['name']['name'])
        df = df[keep_cols]


        # rename columns to 'defaults'
        rename_dict = {val['name'] : key for key, val in info['band_of_use'].items()}
        rename_dict[info['other_mapping']['datetime']['name']] = 'datetime'
        rename_dict[info['other_mapping']['name']['name']] = 'name'
        df = df.rename(columns=rename_dict)

        # unit conversion
        for col in info['conversions'].keys():
            df[col] = df[col] * info['conversions'][col]


        # format datatime
        df["datetime"] = pd.to_datetime(df["datetime"],
                                        format=info['other_mapping']['datetime']['fmt'])

        df["datetime"] = df["datetime"].dt.tz_localize(info['other_mapping']['datetime']['tz'])


        # Make multiidx structure:
        df = df.set_index(['name', 'datetime'])

        # 3. update attributes
        self.df = df
        self.df_tz = info['other_mapping']['datetime']['tz']

        unit_dict = {key: val['units'] for key, val in info['band_of_use'].items() if 'units' in val}
        self._df_units.update(unit_dict)

        self._is_alaro25=True




    def set_model_from_csv(self, csvpath):
        """
        This method loads timeseries data that is stored in a csv file.
        The name of the gee dataset the timeseries are coming from must be the
        same as the .modelname attribute of the Modeldata.


        The timeseries will be formatted and converted to standard toolkit
        units.

        Parameters
        ----------
        csvpath : str
            Path of the csv file containing the modeldata timeseries.

        Returns
        -------
        None.

        """

        # tests ----
        if not self.modelname in self.mapinfo.keys():
            logger.warning(f'{self.modelname} is not found in the gee datasets.')
            return

        # 1. Read csv and set timezone
        df = pd.read_csv(csvpath, sep=",")
        # format datetime
        df["datetime"] = pd.to_datetime(df["datetime"], format="%Y%m%d%H%M%S")
        # (assume all gee dataset are in UTC)
        df["datetime"] = df["datetime"].dt.tz_localize('UTC')
        # self.df_tz='UTC'

        # 2. Format dataframe
        # format index
        df = df.set_index(["name", "datetime"])
        df = df.sort_index()

        # rename to values to toolkit space

        bandname = df.columns[0] #assume only one column
        # scan to the geeinfo to found which obstype and unit the bandname represents
        geeinfo = self.mapinfo[self.modelname]
        obstype = [obs for obs, val in geeinfo['band_of_use'].items() if val['name'] == bandname][0]
        cur_unit = [val['units'] for obs, val in geeinfo['band_of_use'].items() if val['name'] == bandname][0]

        df = df.rename(columns={bandname: obstype})


        # 3. update attributes
        self.df = df[[obstype]]
        self.df_tz = 'UTC'
        self._df_units[obstype] = cur_unit


        # 4. Convert units
        self.convert_units_to_tlk(obstype=obstype,
                                  target_unit_name=standard_tlk_units[obstype])




    def interpolate_modeldata(self, to_multiidx, obstype="temp"):
        """
        Interpolate the modeldata timeseries, of an obstype, to a
        given name-datetime multiindex.

        The modeldata will be converted to the timezone of the multiindex.

        If no interpolation can be done, Nan values are used.

        Parameters
        ----------
        to_multiidx : pandas.MultiIndex
            A name - datetime (tz-aware) multiindex to interpolate the
            modeldata timeseries to.
        obstype : str, optional
            Observation type of the timeseries. obstype must be a column in the
            Modeldata.df. The default is "temp".

        Returns
        -------
        returndf : pandas.DataFrame
            A dataframe with to_multiidx as an index and obstype as a column.
            The values are the interpolated values.

        """
        returndf = init_multiindexdf()

        recordsdf = init_multiindexdf()
        recordsdf.index = to_multiidx
        # iterate over stations check to avoid extrapolation is done per stations
        for sta in recordsdf.index.get_level_values("name").unique():
            sta_recordsdf = xs_save(recordsdf, sta, level="name", drop_level=False)
            sta_moddf = xs_save(self.df, sta, level="name", drop_level=False)

            # convert modeldata to timezone of observations
            sta_moddf = conv_tz_multiidxdf(
                df=sta_moddf,
                timezone=sta_recordsdf.index.get_level_values("datetime").tz,
            )

            # check if modeldata is will not be extrapolated !
            if min(sta_recordsdf.index.get_level_values("datetime")) < min(
                sta_moddf.index.get_level_values("datetime")
            ):
                logger.warning("Modeldata will be extrapolated")
            if max(sta_recordsdf.index.get_level_values("datetime")) > max(
                sta_moddf.index.get_level_values("datetime")
            ):
                logger.warning("Modeldata will be extrapolated")

            # combine model and records
            mergedf = sta_recordsdf.merge(
                sta_moddf, how="outer", left_index=True, right_index=True
            )

            # reset index for time interpolation
            mergedf = mergedf.reset_index().set_index("datetime").sort_index()

            # interpolate missing modeldata
            mergedf[obstype].interpolate(
                method="time", limit_area="inside", inplace=True
            )
            # convert back to multiindex
            mergedf = mergedf.reset_index().set_index(["name", "datetime"]).sort_index()
            # filter only records
            mergedf = mergedf.loc[sta_recordsdf.index]

            returndf = pd.concat([returndf, mergedf])
        return returndf

    def make_plot(self, obstype_model="temp", dataset = None,
                  obstype_dataset=None, stationnames=None,
                  starttime=None, endtime=None, title=None, show_outliers=True,
                  show_filled=True, legend=True,
                  _ax=None, #needed for GUI, not recommended use
                  ):

        """
        This function creates a timeseries plot for the Modeldata. When a
        metobs_toolkit.Dataset is provided, it is plotted in the same figure.

        The line colors represent the timesries for different locations.



        Parameters
        ----------
        obstype_model : string, optional
             Fieldname of the Modeldata to visualise. The default is 'temp'.
        dataset : metobs_toolkit.Dataset, optional
            A Dataset instance with observations plotted in the same figure.
            Observations are represented by solid line and modeldata by dashed
            lines. The default is None.
        obstype_dataset : string, optional
            Fieldname of the Dataset to visualise. Only relevent when a dataset
            is provided. If None, obsype_dataset = obstype_model. The default
            is None.
        stationnames : list, optional
            A list with stationnames to include in the timeseries. If None is
            given, all the stations are used, defaults to None.
        starttime : datetime.datetime, optional
             Specifiy the start datetime for the plot. If None is given it will
             use the start datetime of the dataset, defaults to None.
        endtime : datetime.datetime, optional
             Specifiy the end datetime for the plot. If None is given it will
             use the end datetime of the dataset, defaults to None.
        title : string, optional
             Title of the figure, if None a default title is generated. The
             default is None.
        show_outliers : bool, optional
             If true the observations labeld as outliers will be included in
             the plot. Only relevent when a dataset is provided. The default
             is True.
        show_filled : bool, optional
             If true the filled values for gaps and missing observations will
             be included in the plot. Only relevent when a dataset is provided.
             The default is True.
        legend : bool, optional
             If True, a legend is added to the plot. The default is True.


        Returns
        -------
        axis : matplotlib.pyplot.axes
             The timeseries axes of the plot is returned.

        """


        logger.info(f"Make {obstype_model}-timeseries plot of model data")

        # Basic test
        if obstype_model not in self.df.columns:
            logger.warning(f'{obstype_model} is not foud in the modeldata df.')
            return
        if self.df.empty:
            logger.warning('The modeldata is empty.')
            return
        if obstype_dataset is None:
            obstype_dataset = obstype_model

        if (not dataset is None):
            if (obstype_dataset not in dataset.df.columns):
                logger.warning(f'{obstype_dataset} is not foud in the Dataframe df.')
                return


        model_df = self.df

        # ------ filter model ------------

        # Filter on obstype
        model_df = model_df[[obstype_model]]

        # Subset on stationnames
        if not stationnames is None:
            model_df = model_df[model_df.index.get_level_values('name').isin(stationnames)]

        # Subset on start and endtime
        model_df = multiindexdf_datetime_subsetting(model_df, starttime, endtime)


        #  -------- Filter dataset (if available) -----------
        if not dataset is None:
            # combine all dataframes
            mergedf = dataset.combine_all_to_obsspace()

            # subset to obstype
            mergedf = xs_save(mergedf, obstype_dataset, level='obstype')

            # Subset on stationnames
            if not stationnames is None:
                mergedf = mergedf[mergedf.index.get_level_values('name').isin(stationnames)]

            # Subset on start and endtime
            mergedf = multiindexdf_datetime_subsetting(mergedf, starttime, endtime)


        # Generate ylabel

        try:
            model_true_field_name = self.mapinfo[self.modelname]['band_of_use'][obstype_model]['name']
        except KeyError:
            logger.info(f'No model field name found for {obstype_model} in {self}.')
            model_true_field_name = 'Unknown fieldname'

        fieldname = f'{model_true_field_name}'

        if not dataset is None:
            dataset_obs_orig_name = dataset.data_template[obstype_dataset]['orig_name']
            units = dataset.data_template[obstype_dataset]['units']
            y_label = f'{fieldname} \n {dataset_obs_orig_name} ({units})'

        else:

            y_label = f'{fieldname} \n ({self._df_units[obstype_model]})'


        #Generate title
        title = f'{self.modelname} : {model_true_field_name}'
        if not dataset is None:
            title = f'{title} and {dataset_obs_orig_name} observations.'


        # make plot of the observations
        if not dataset is None:
            # make plot of the observations
            ax, col_map = timeseries_plot(mergedf=mergedf,
                                 title=title,
                                 ylabel=y_label,
                                 colorby='name',
                                 show_legend=legend,
                                 show_outliers=show_outliers,
                                 show_filled=show_filled ,
                                 settings = dataset.settings)


            # Make plot of the model on the previous axes
            ax, col_map = model_timeseries_plot(
                                    df=model_df,
                                    obstype=obstype_model,
                                    title=title,
                                    ylabel=y_label,
                                    settings = self._settings,
                                    show_primary_legend=False,
                                    add_second_legend=True,
                                    _ax = ax,
                                    colorby_name_colordict=col_map)


        else:

            # Make plot of model on empty axes
            ax, _colmap = model_timeseries_plot(
                    df=model_df,
                    obstype=obstype_model,
                    title=title,
                    ylabel=y_label,
                    settings = self._settings,
                    show_primary_legend=legend,
                    add_second_legend=False,
                    _ax = None
                    )

        return ax