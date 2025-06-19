import logging

import pandas as pd
import geopandas as gpd
import xarray as xr
import regionmask

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy
import numpy as np
from sklearn.neighbors import BallTree

from metobs_toolkit.backend_collection.df_helpers import to_timedelta

from metobs_toolkit.nwp_collection.error_classes import MetObsFieldNotFound
from metobs_toolkit.nwp_collection.xr_helpers import variable_names, coord_names
from metobs_toolkit.nwp_collection.field_defenitions import default_SFX_fields
from metobs_toolkit.modeltimeseries import ModelTimeSeries



# Use logger with name "<metobs_toolkit>"
logger = logging.getLogger("<metobs_toolkit>")

class Grid:
    def __init__(self, gridds):
        self.gridds = gridds.compute()
        #setup gridpoint index
        self.gp_index, self.gp_1d = self._construct_gp_index()


    def _construct_gp_index(self):

        #Create a array of [lat,lon] pairs for each gridpoint
        latdata = self.gridds.lat.data
        londata = self.gridds.lon.data

        d3comb = np.stack((latdata, londata), axis=-1) #as a 3D array
        d1comb = d3comb.reshape(-1, d3comb.shape[-1]) #as a 1D array of pairs

        #Create an indexer, with 'haversine' metric
        index = BallTree(np.deg2rad(d1comb), 
                    leaf_size=40, #default
                    metric='haversine', #haversine approx metric for latlon space
                    #other arguments are passed to the metric
                    )
        return index, d1comb
    
    def _get_nearest_1d_index(self, trg_lat, trg_lon):

        #Find the nearest gridpoint to the target point
        trgpoint = [[trg_lat, trg_lon]]
        

        #find nearset index to target locations
        dist, trgidx1D = self.gp_index.query(
                                np.deg2rad(trgpoint), 
                                k=1, #only the NN
                                return_distance=True)
        
        # metric_dist_approx= dist*63781370000 #earth radius in meters
        return trgidx1D.ravel()[0]
        #get lat,lon pair of nearest gp
        return nn_gp[0], nn_gp[1]
    
    def _get_xy_of_1d_index(self, indexint):
        latdata = self.gridds.lat.data
        londata = self.gridds.lon.data
        d3comb = np.stack((latdata, londata), axis=-1)
        i,j = np.unravel_index(indexint, d3comb.shape[:-1])

        x = self.gridds['x'].data[i]
        y =  self.gridds['y'].data[j]
        return x, y
    def get_nearest_gp(self, trg_lat, trg_lon):
        #Get nearest gridpoint in 1D space (as index)
        idx_1d= self._get_nearest_1d_index(
                                        trg_lat=trg_lat,
                                        trg_lon = trg_lon)
        
        #Get the nearest gridpoint in lat-lon space
        # nn_gp = self.gp_1d[idx_1d]

        #Get the nearest gridpoint in model-space (x, y)
        nn_x, nn_y = self._get_xy_of_1d_index(idx_1d)

        #Safty, recompute the nnlat and nnlon from the gridpoint x,y.
        #so the the xy and latlon coordinates are surely linked 
        nnlat = float(self.gridds['lat'].sel(x=nn_x, y=nn_y).data)
        nnlon = float(self.gridds['lon'].sel(x=nn_x, y=nn_y).data)
        
        return {"nnlat": nnlat,
                "nnlon": nnlon,
                "gridpoint_x": int(nn_x),
                "gridpoint_y": int(nn_y)}
        
    


class ModelDataset:
    """
    Wrapper for an xarray.Dataset representing a model dataset.

    Ensures that the dataset contains 'latitude', 'longitude', and 'validtime' coordinates.
    """

    REQUIRED_COORDS = ['lat', 'lon', 'validtime', 'reference_time']

    def __init__(self, modelID:str, dataset: xr.Dataset, field_defenitions: list):
        """
        Initialize the ModelDataset.

        Parameters
        ----------
        dataset : xarray.Dataset
            The xarray Dataset to wrap.

        Raises
        ------
        ValueError
            If any of the required coordinates are missing.
        """

        self.modelID = modelID
        self.dataset = dataset
        self.field_defs = field_defenitions

        #Test and format the datset
        self._check_required_coords()
        self._subset_to_mapped_fields()
        self._too_standard_units()

        #Set grid
        self.grid = Grid(self.dataset[['lat', 'lon']])

        self._tails_are_removed = False

    def __repr__(self):
        return f"Gridded output of {self.modelID}:\n {self.dataset}"
    
    @property
    def variable_names(self):
        return variable_names(self.dataset)
   
    def _subset_to_mapped_fields(self):

        in_dataset = self.variable_names
        
        in_defs = [modlobs.model_band for modlobs in self.field_defs]

        #log unmapped
        droped_vars = list(set(in_dataset) - set(in_defs))
        if bool(droped_vars):
            logger.warning(f'The following variables are unmapped and are removed from the modeldataset: {droped_vars}')

        #Do not drop the coordinates
        target_subset = set(in_dataset).intersection(set(in_defs))
        self.dataset = self.dataset[list(target_subset) + list(coord_names(self.dataset))]
        return
    
    def _too_standard_units(self):
        for field in self.variable_names:
            print(field)
            modelobs = self._modelobstype_from_variablename(field)
            self.dataset[field].data = modelobs.convert_to_standard_units(
                                        input_data = self.dataset[field].data,
                                        input_unit = modelobs.model_unit
                                         )

    def _cleanup_coordinates(self):
        ### when subsetting to the variables, only the coordinates
        # used by the variables are kept
        self.dataset = self.dataset[self.variable_names]
        return
    
    def _check_required_coords(self):
        missing = [coord for coord in self.REQUIRED_COORDS if coord not in self.dataset.coords]
        if missing:
            raise ValueError(f"Missing required coordinates: {missing}")
        
        # Ensure 'validtime' is a dimension
        if 'validtime' not in self.dataset.dims:
            self.dataset = self.dataset.expand_dims('validtime')

        # Ensure 'validtime' is an index
        # Needed for the .sel() method on validtime
        if 'validtime' not in self.dataset.indexes:
            self.dataset = self.dataset.set_index({'validtime': 'validtime'})

    def spatial_plot(self,
                    target_fieldname: str | None=None, 
                    validtime: pd.Timestamp | None = None,
                    **kwargs,
                    ):
        
        #Set fieldname
        if target_fieldname is None:
            #pick the first
            target_fieldname = self.variable_names[0]
        else:
            self._check_fieldname_is_known(target_fieldname)


        #Set validtime
        if validtime is None:
            #Pick the first present validtime
            validtime = self.dataset['validtime'].data[0]
        else: 
            #test if validtime is in the dataset
            if not isinstance(validtime, pd.Timestamp):
                raise TypeError("validtime must be a pandas.Timestamp")
            if validtime not in self.dataset['validtime'].values:
                raise ValueError(f"validtime {validtime} not found in dataset validtime dimension")

        #Construct Xr array (2D repr)
        to_plot = self.dataset[target_fieldname].sel({'validtime': validtime})

        # Set up a standard map for latlon data.
        fig, geo_axes= plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
        # Make plot
        to_plot.plot.pcolormesh(
                x='lon',
                y='lat',
                transform=ccrs.PlateCarree(),
                cbar_kwargs={"orientation": "horizontal", "shrink": 0.8, "aspect": 40},                                                     
                robust=True,
                ax=geo_axes,
                **kwargs,
                )    
        
        trgmodelobstype = self._modelobstype_from_variablename(target_fieldname)


        geo_axes.set_title(f'{trgmodelobstype} at validate {validtime}')
        #adding features
        geo_axes.coastlines()
        geo_axes.add_feature(cartopy.feature.BORDERS)

        #scaling    
        geo_axes.autoscale_view()
        return geo_axes




    def _get_nn_gridpoint(self, trg_lat, trg_lon):

        return self.grid.get_nearest_gp(trg_lat=trg_lat, trg_lon=trg_lon)
        

    def _modelobstype_from_variablename(self, variabelname):
        trgmodelobscandiates = [modelobs for modelobs in self.field_defs if modelobs.model_band == variabelname]
        if not bool(trgmodelobscandiates):
            raise MetObsFieldNotFound(f'{variabelname} is not known field.')
        
        return trgmodelobscandiates[0]
    
    def _check_fieldname_is_known(self, fieldname):
        if fieldname not in self.variable_names:
            raise MetObsFieldNotFound(f'{fieldname} not a known fieldname. These fields are present: {self.variable_names}')
        return

    def insert_modeltimeseries(
            self,
            stationlist: list,
            target_variables: list,
            compute_before_assign=False,
            force_update: bool = False) -> None: 
        

        #Update the stations sites
        stations_grid_info = {
                'name':[],
                'x':[],
                'y':[]}

        logger.info('finding nn gp for all stations')
        for sta in stationlist:
            #Get the nn gridpoint for the station
            sta.site._grid_info = self._get_nn_gridpoint(
                                        trg_lat=sta.site.lat,
                                        trg_lon=sta.site.lon)
            #Add it to the dict
            stations_grid_info['name'].append(sta.name)
            stations_grid_info['x'].append(sta.site.get_gp_x_index())
            stations_grid_info['y'].append(sta.site.get_gp_y_index())

        logger.info('clipping nn from fc')
        # create xr dataset with station locations
        targetds = xr.Dataset(pd.DataFrame(data=stations_grid_info).set_index('name'))
        
        #clip the stations from the model
        targetds = self.dataset.sel(targetds)

        def getmemsize(obj):
            mem_bytes = obj.nbytes if hasattr(obj, 'nbytes') else obj.__sizeof__()
            mem_gb = mem_bytes / (1024 ** 3)
            return 
        logger.info(f"Memory size of targetds: {getmemsize(targetds):.3f} GB")

        if compute_before_assign:
            logger.info('Computing the targetds, since compute_before_assign is true')
            targetds = targetds.compute()
            logger.info(f"Memory size of targetds after compute: {getmemsize(targetds):.3f} GB")


        logger.info(f"Memory size of self: {getmemsize(self):.3f} GB")
        #Memory is more limitting then computatio time, so favour memory efficient
        i = 1
        for sta in stationlist:
            logger.info(f'creating modeldata for {sta.name} ({i}/{len(stationlist)+1}) ')
            logger.info(f"Memory size of sta before modelfil: {getmemsize(sta):.3f} GB")
            for var in target_variables:
                 sta.add_to_modeldata(
                     ModelTimeSeries(
                            site=sta.site,
                            datarecords=targetds[var].sel(name=sta.name).to_numpy(),
                            timestamps=targetds['validtime'].to_numpy(),
                            modelobstype=self._modelobstype_from_variablename(var),
                            datadtype=np.float32,
                            timezone='UTC',
                            modelID=self.modelID,
                            ),
                            force_update=force_update)
            logger.info(f"Memory size of sta after modelfil: {getmemsize(sta):.3f} GB")
            logger.info(f"Memory size of self after modelfil of {sta.name}: {getmemsize(self):.3f} GB")                
            i+=1
        #Construct Dataframe
        # trg_index = ['name', 'validtime'] #TODO add other dimensions for cycle applications
        # targetdf = targetds[trg_index + target_variables].compute().to_dataframe()
        # targetdf = targetdf.reset_index().set_index(['name', 'validtime'])
        # targetdf = targetdf[target_variables] #drop columns representing coordinates
        
        # #add it as modeltimeseries to the stations
        # for sta in stationlist:
        #     print(f'creating modeldata for {sta.name}')
        #     stadf = targetdf.xs(sta.name, level='name')
        #     #each column represents bandvalues
        #     for colname in stadf.columns:

        #         sta.add_to_modeldata(ModelTimeSeries(
        #                     site=sta.site,
        #                     datarecords=stadf[colname].to_numpy(),
        #                     timestamps=stadf.index.to_numpy(),
        #                     modelobstype=self._modelobstype_from_variablename(colname),
        #                     datadtype=np.float32,
        #                     timezone='UTC',
        #                     modelID=self.modelID,
        #                     ),
        #                     force_update=force_update)
        #     #Save memory
        #     del stadf
        #     targetdf = targetdf.drop(sta.name, level='name')
        



    def trim_to_box(self, minlat, maxlat, minlon, maxlon, drop=False):
        self.dataset = self.dataset.where(
                    ((self.dataset['lat'] <= maxlat) &
                        (self.dataset['lat'] >= minlat) &
                        (self.dataset['lon'] <= maxlon) &
                        (self.dataset['lon'] >= minlon)),
                    drop=drop)
        

    def trim_to_country(self, country_shp_file,
                        country='Belgium',
                        drop=False):
        
        #read the shape file
        shpfile = gpd.read_file(country_shp_file)

        #subset one country
        trg_shp = shpfile.loc[shpfile['name'] == country]
        
        #convert trg_shp to  regionmask
        trg_region = regionmask.from_geopandas(trg_shp) #to region
        
        #rasterize (convert to a dataset with nan's out of the trg region)
        maskraster = trg_region.mask(self.dataset,
                                        lon_name = 'lon',
                                        lat_name='lat')
        
        self.dataset = self.dataset.where(maskraster.notnull(), drop=drop)
        
    def trim_spinup_and_tails(self, spinup_duration=pd.Timedelta('4h'),
                              cycle_freq:pd.Timestamp|None = None,
                              validtimename='validtime',
                              referecetimename='reference_time',
                              drop=True):
        #1. add leadtime
        self.add_leadtime_dimension(validtimename=validtimename,
                                    referecetimename=referecetimename)

        # 2. construct limits of leadtime
        if cycle_freq is None:
            #estimate the cycle frequency
            cycle_freq = pd.infer_freq(pd.to_datetime(self.dataset[referecetimename]))
        
            
        #to timedeltas
        ref_time_freq = to_timedelta(cycle_freq)
        spinup = to_timedelta(spinup_duration)

        min_leadtime = spinup
        max_leadtime = ref_time_freq + spinup

        # 3. subset data
        #NOTE: Make sure that always one of the bounds is included and the other
        # is not! Else there will be multiple values for one validtime (based
        # on referene time).
        self.dataset = self.dataset.where(
            (self.dataset['leadtime'] >= min_leadtime) & #leadtime included
            (self.dataset['leadtime'] < max_leadtime), #Maxleadtime excluded !!! 
            drop=drop)
        
        #4. Reconstruct the 1D time timensions

        #Note the leadtime is collapsed with another method, since the non-physical
        #elements on this coordinate are not Nan (but negative or in-spinup values)
        self.dataset['leadtime1d'] = self.dataset['leadtime'].where(self.dataset['leadtime'] >=min_leadtime).min(dim='reference_time', skipna=True)
        self.dataset = self.dataset.assign_coords(leadtime1d = self.dataset['leadtime1d'])
        self.dataset = self.dataset.drop_vars('leadtime')
        self.dataset = self.dataset.rename({'leadtime1d': 'leadtime'})

        # Remove reference_time by collapsing it, the agg function is irrelevant if it skips nan's
        self.dataset = self.dataset.ffill(dim='reference_time').isel(reference_time=-1)
        
        #recreate reference_time
        self.dataset = self.dataset.drop_vars('reference_time')
        self.dataset = self.dataset.assign_coords(reference_time=
                    (self.dataset['validtime'] - self.dataset['leadtime']))

    # def trim_spinup_period(self, spinup_duration=pd.Timedelta('4h')):
    #     if self._tails_are_removed:
    #         raise RuntimeError("trim_spinup_period must be applied before the tails are removed. Call trim_spinup_period before trim_tails_of_cycled_ds.")

    #     if isinstance(spinup_duration, str):
    #         #convert to pandas Timedelta
    #         spinup_duration = pd.Timedelta(spinup_duration)

    #     self.dataset = self.dataset.where(  
    #                 self.dataset['validtime'] - self.dataset['reference_time'] >= spinup_duration,
    #                 drop=True) #Drop must be true for the tail removal to work
        
    

    def add_leadtime_dimension(self,
                        validtimename='validtime',
                        referecetimename='reference_time'):
        #Create leadtime dimension
        leadtime = (self.dataset['validtime'] - self.dataset['reference_time']).astype('timedelta64[s]')
        self.dataset = self.dataset.assign_coords(leadtime=leadtime)


    # def trim_tails_of_cycled_ds(self,
    #                             validtimename='validtime',
    #                             referecetimename='reference_time'):
    #     """ Drop the tails of data that are captured by a newer cycle for all variables. """
    #     #TODO: remove spinup first !!! 
    #     self._tails_are_removed = True
    #     self.add_leadtime_dimension()

    #     #ASSUMPTION: the time-related dimensions are the same for all variables !!!!
    #     # this is typical true for FA files, which contains variables at the same validtime and refernce_time. 
       
    #     def get_no_tail_indices(da):
    #         # Mask out negative leadtimes
    #         masked = da.where(da['leadtime'] >= 0)
    #         # If all values are NaN, return 0 or np.nan (choose what makes sense for your use case)
    #         if masked.isnull().all():
    #             # Option 1: return np.nan (will drop this group later)
    #             # return np.nan
    #             return xr.DataArray(np.nan)
    #             # Exclude 'reference_time' from coords and dims
    #             # Option 2: return 0 (will select the first, but may not be what you want)
    #             # return xr.DataArray(0, coords=masked.coords, dims=[])
    #         # Otherwise, return the index of the minimum
    #         return masked.argmin(dim=referecetimename)
    #     #Find no-tail indices on a single field (Speedup + memory saving)
    #     dummy_field = self.variable_names[0]
    #     no_tail_idices = self.dataset[dummy_field].groupby(validtimename).map(get_no_tail_indices)

    #     #Drop nans
    #     no_tail_idices = no_tail_idices.dropna(validtimename)

    #     if  hasattr(no_tail_idices, "compute"):
    #         #issue is that isel() does not take chunked dask array as input, so they need
    #         #to be computed 
    #         no_tail_idices = no_tail_idices.compute()

    #     self.dataset = self.dataset.isel({referecetimename: no_tail_idices})
       
    #     self._tails_are_removed = True #To be checked when calling remove spinup after tail removal
    #     #Fix the leadtime  + refernce_tim coordinate
    #     # for some reason it is strangled with x and y dimension, so unravel
    #     # and make it dependant only on validtime 

    #     new_lt = self.dataset['leadtime'].isel(x=1, y=1)
    #     self.dataset = self.dataset.drop_vars(['leadtime', 'reference_time'])
    #     #Construct leadtime coord (1D depending on validtime)
    #     self.dataset = self.dataset.assign_coords(leadtime=("validtime", new_lt.data))
    #     #Construct reference time coord (1D depending on validtime)  
    #     self.dataset = self.dataset.assign_coords(reference_time=(
    #         ("validtime",
    #          (self.dataset.validtime - new_lt).data)
    #     )) 
        
        


    