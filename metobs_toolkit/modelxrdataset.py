import logging

import pandas as pd
import xarray as xr

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
        self.dataset = self.dataset[in_defs + coord_names(self.dataset)]
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
        


   