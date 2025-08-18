import numpy as np
from sklearn.metrics import mean_squared_error, mean_absolute_error


def rmse_func(obs, model) -> float:
    """Calculate the Root Mean Square Error (RMSE) between observed and model values."""
    return np.sqrt(mean_squared_error(obs, model))

def mae_func(obs, model) -> float:
    return mean_absolute_error(obs, model)

def bias_func(obs, model) -> float:
    """Calculate the bias between observed and model values."""
    return np.mean(model - obs)

def N_samples(obs, model) -> int:
    """Calculate the number of samples."""
    return np.count_nonzero(~np.isnan(obs) & ~np.isnan(model))