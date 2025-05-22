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





class ModelDataset:
    """
    Wrapper for an xarray.Dataset representing a model dataset.

    Ensures that the dataset contains 'latitude', 'longitude', and 'validtime' coordinates.
    """

    REQUIRED_COORDS = ['lat', 'lon', 'validtime', 'reference_time']

    def __init__(self, ID:str, dataset: xr.Dataset, field_defenitions: list):
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

        self.ID = ID
        self.dataset = dataset
        self.field_defs = field_defenitions

        #Test and format the datset
        self._check_required_coords()
        self._subset_to_mapped_fields()
        self._too_standard_units()

    def __repr__(self):
        return str(self.dataset)
    
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
            validtime = dataset['validtime'].data[0]
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

        #IDEA:  The creattion of the index is recalculated for each nn computation,
        #this is not needed, it is overkill

        #Create a array of [lat,lon] pairs for each gridpoint
        latdata = self.dataset.lat.data.compute()
        londata = self.dataset.lon.data.compute()

        d3comb = np.stack((latdata, londata), axis=-1) #as a 3D array
        d1comb = d3comb.reshape(-1, d3comb.shape[-1]) #as a 1D array of pairs

        #Create an indexer, with 'haversine' metric
        index = BallTree(np.deg2rad(d1comb), #RADIANS for haversine !! 
                    leaf_size=40, #default
                    metric='haversine', #haversine approx metric for latlon space
                    #other arguments are passed to the metric
                    )
        
        #Find the nearest gridpoint to the target point
        trgpoint = [[trg_lat, trg_lon]]
        #find nearset index to target locations
        dist, trgidx1D = index.query(
                                np.deg2rad(trgpoint), #RADIANS for haversine
                                k=1, #only the NN
                                return_distance=True)
        
        metric_dist_approx= dist*63781370000 #earth radius in meters

        nnlat = d1comb[trgidx1D.ravel()[0]][0]
        nnlon = d1comb[trgidx1D.ravel()[0]][1]

        # Get the x, y of the gridpoint
        gridpoint_x, gridpoint_y = np.argwhere(latdata == nnlat)[0]
        

        return {
            "nnlat": float(nnlat),
            "nnlon": float(nnlon),
            "gridpoint_x": int(gridpoint_x),
            "gridpoint_y": int(gridpoint_y),
            "distance_m": float(metric_dist_approx.ravel()[0])
        }

    def _modelobstype_from_variablename(self, variabelname):
        trgmodelobscandiates = [modelobs for modelobs in self.field_defs if modelobs.model_band == variabelname]
        if not bool(trgmodelobscandiates):
            raise MetObsFieldNotFound(f'{variabelname} is not known field.')
        
        return trgmodelobscandiates[0]
    
    def _check_fieldname_is_known(self, fieldname):
        if fieldname not in self.variable_names:
            raise MetObsFieldNotFound(f'{fieldname} not a known fieldname. These fields are present: {self.variable_names}')
        return





    def extract_modeltimeseries(self, station, trg_variable) -> None:
        self._check_fieldname_is_known(trg_variable)
        #subset the grid to single gp
        pointds = self.dataset.sel(
                        {'x': station.site.get_gp_x_index(),
                         'y': station.site.get_gp_y_index()}).compute()
        
        #Extract timeseries
        return ModelTimeSeries(
                    site=station.site,
                    datarecords=pointds[trg_variable].data,
                    timestamps=pointds.validtime.data, 
                    obstype=self._modelobstype_from_variablename(trg_variable),
                    datadtype=np.float32,
                    timezone='UTC',
                    modelname=self.ID,
                    modelvariable=trg_variable)
            