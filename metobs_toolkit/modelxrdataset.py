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
            force_update: bool = False) -> None: 
        

        #Update the stations sites
        stations_grid_info = {
                'name':[],
                'x':[],
                'y':[]}


        for sta in stationlist:
            #Get the nn gridpoint for the station
            sta.site._grid_info = self._get_nn_gridpoint(
                                        trg_lat=sta.site.lat,
                                        trg_lon=sta.site.lon)
            #Add it to the dict
            stations_grid_info['name'].append(sta.name)
            stations_grid_info['x'].append(sta.site.get_gp_x_index())
            stations_grid_info['y'].append(sta.site.get_gp_y_index())


        # create xr dataset with station locations
        targetds = xr.Dataset(pd.DataFrame(data=stations_grid_info).set_index('name'))

        #clip the stations from the model
        targetds = self.dataset.sel(targetds)

        #Construct Dataframe
        trg_index = ['name', 'validtime'] #TODO add other dimensions for cycle applications
        targetdf = targetds[trg_index + target_variables].compute().to_dataframe()
        targetdf = targetdf[target_variables] #drop columns representing coordinates
        
        #add it as modeltimeseries to the stations
        for sta in stationlist:
            stadf = targetdf.xs(sta.name, level='name')
            #each column represents bandvalues
            for colname in stadf.columns:

                sta.add_to_modeldata(ModelTimeSeries(
                            site=sta.site,
                            datarecords=stadf[colname].to_numpy(),
                            timestamps=stadf.index.to_numpy(),
                            modelobstype=self._modelobstype_from_variablename(colname),
                            datadtype=np.float32,
                            timezone='UTC',
                            modelID=self.modelID,
                            ),
                            force_update=force_update)
        


   

    # package_root=Path(__file__).parent
    # #file with country shapes
    # country_shp = os.path.join(package_root, 'data', 'world-administrative-boundaries.shp')


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
        
        


    def trim_spinup_period(self, spinup_duration=pd.Timedelta('4h'), drop=False):
        if isinstance(spinup_duration, str):
            #convert to pandas Timedelta
            spinup_duration = pd.Timedelta(spinup_duration)
        
        self.dataset = self.dataset.where(  
                    self.dataset['validtime'] - self.dataset['reference_time'] >= spinup_duration,
                    drop=drop)
    


    def trim_tails_of_cycled_ds(self, keep_smallest_lt=True,
                                validtimename='validtime',
                                referecetimename='reference_time'):
        """ Drop the tails of data that are captured by a newer cycle for all variables. """
        #WARNING! make shure that the leadtime is extracted before executing this function !! 
        if keep_smallest_lt:
            select_duplicate_idx = -1 #last from duplicates
        else: 
            select_duplicate_idx = 1 #first from duplicates


        #attempt
        #reshape so that reference time is last dimension
        self.dataset = self.dataset.transpose(..., validtimename, referecetimename)  # equivalent
        #sort data by all dimensions
        self.dataset = self.dataset.sortby([validtimename, referecetimename])


        for field in self.dataset.variables: 
            if field == referecetimename:
                continue #remove after all references in the fields are removed
        
            if referecetimename in self.dataset[field].dims:
                #trgdims: the order of dimensions, referencetimename must be the last dimension !!! 
                #reshape so that reference_time is last dimension
                # trgdims = tuple(set(self.dataset[field].dims) - set([referecetimename]))+tuple([referecetimename])
                # fieldds = self.dataset[field].transpose(*trgdims)
            
                #sort data by all dimensions
                # fieldds = fieldds.sortby(variables = list(trgdims))

                #reduce the data
                # newdata = _reduce_last_dimension(arr=fieldds.data,
                #                                 keepidx=select_duplicate_idx)
                newdata = _reduce_last_dimension(arr=self.dataset[field].data,
                                                keepidx=select_duplicate_idx)

                #construct a dataarray
                newcoords = {coord: self.dataset[field][coord] for coord in self.dataset[field].coords if coord != referecetimename}
                newfieldds = xr.DataArray(
                    data=newdata,
                    coords=newcoords,
                    dims=self.dataset[field].dims[:-1],  # drop the last dimension
                    attrs=self.dataset[field].attrs
                )  # the last dimension is resolved
            
                #overwrite the field
                self.dataset[field] = newfieldds
        
        # Drop the reference_time dimension from the dataset
        self.dataset = self.dataset.drop_dims(referecetimename)

# Reduce the last dimension by taking the last non-NaN element along that dimension
def _reduce_last_dimension(arr, keepidx=-1):
    """
    Reduce the last dimension of a numpy array by selecting the keepidx-th non-NaN value
    along the last dimension for each index of the other dimensions.

    Parameters
    ----------
    arr : np.ndarray
        Input array. Reduction is performed along the last dimension.
    keepidx : int
        Index to select among the non-NaN values along the last dimension.
        -1 means the last non-NaN value (default).

    Returns
    -------
    reduced : np.ndarray
        Array with the last dimension reduced.
    """
    arr = np.asarray(arr)
    # arr shape: (..., N)
    # mask for non-NaN values
    mask = ~np.isnan(arr)
    # count valid values along last dimension
    valid_count = mask.sum(axis=-1)
    # indices of the last valid value along last dimension
    if keepidx == -1:
        # last valid
        idx = np.flip(mask, axis=-1).argmax(axis=-1)
        idx = arr.shape[-1] - 1 - idx
    else:
        # kth valid (from start)
        # get indices of all valid values
        valid_idx = np.where(mask)
        # create an array of indices for the last dimension
        last_dim_idx = np.arange(arr.shape[-1])
        # broadcast to shape of arr
        last_dim_idx = np.broadcast_to(last_dim_idx, arr.shape)
        # for each position, get the indices of valid values
        # then select the keepidx-th one
        # This is tricky to vectorize for arbitrary keepidx, so fallback to slower method
        # but only for non-default keepidx
        reduced = np.full(arr.shape[:-1], np.nan)
        it = np.nditer(valid_count, flags=['multi_index'])
        for _ in it:
            idxs = np.where(mask[it.multi_index])[0]
            if idxs.size > 0 and abs(keepidx) < idxs.size:
                reduced[it.multi_index] = arr[it.multi_index + (idxs[keepidx],)]
        return reduced


    # Use advanced indexing to select the value at idx for each position
    # Build indices for all dimensions except last
    grid = list(np.ogrid[tuple(map(slice, arr.shape[:-1]))])
    # Add the last dimension indices
    grid.append(idx)
    reduced = arr[tuple(grid)]
    # Set to nan where there are no valid values
    reduced[valid_count == 0] = np.nan
    return reduced
   


# def _reduce_last_dimension(arr, keepidx = -1):
#     # Create an array to store the reduced values
#     reduced = np.full(arr.shape[:-1], np.nan)
#     for idx in np.ndindex(arr.shape[:-1]):
#         # Extract the last non-NaN value along the last dimension
#         valid_values = arr[idx].compressed() if isinstance(arr[idx], np.ma.MaskedArray) else arr[idx][~np.isnan(arr[idx])]
#         if valid_values.size > 0:
#             reduced[idx] = valid_values[keepidx] # -1: pick last (THus newest cycle)
#     return reduced