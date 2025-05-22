import xarray as xr

def variable_names(ds):
    return list(set(ds.variables.keys()) - set(ds.coords.keys()) )

def coord_names(ds):
    return list(ds.coords.keys())

