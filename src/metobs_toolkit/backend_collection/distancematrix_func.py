"""
Distance matrix generation using BallTree with haversine distance metric.
"""

import pandas as pd
import numpy as np
from sklearn.neighbors import BallTree

from metobs_toolkit.backend_collection.decorators import log_entry


@log_entry
def generate_distance_matrix(
    metadf: pd.DataFrame, lat_col: str = "lat", lon_col: str = "lon"
) -> pd.DataFrame:
    """
    Generate a distance matrix using BallTree with haversine distance metric.

    This function calculates great-circle distances between stations using the
    haversine formula, which is appropriate for calculating distances on the
    Earth's surface given latitude and longitude coordinates.

    Parameters
    ----------
    metadf : pd.DataFrame
        DataFrame with station names as index and latitude/longitude as columns.
        Coordinates should be in decimal degrees.
    lat_col : str, optional
        Name of the latitude column. Default is 'lat'.
    lon_col : str, optional
        Name of the longitude column. Default is 'lon'.

    Returns
    -------
    pd.DataFrame
        Symmetric distance matrix with station names as both index and columns.
        Distances are in meters. The diagonal contains zeros (distance from
        each station to itself).

    Raises
    ------
    ValueError
        If required columns are missing, DataFrame is empty, or coordinates
        contain missing values.

    Notes
    -----
    The haversine distance metric calculates great-circle distances between points
    on a sphere given their latitude and longitude in radians. The Earth's radius
    is assumed to be 6371 km (637100 m). The BallTree algorithm provides efficient
    nearest neighbor queries for geographic data.

    Examples
    --------
    >>> import pandas as pd
    >>> metadf = pd.DataFrame({
    ...     'lat': [50.8503, 51.2194, 50.4501],
    ...     'lon': [4.3517, 4.4025, 3.9538]
    ... }, index=['Brussels', 'Antwerp', 'Ghent'])
    >>> distance_matrix = generate_distance_matrix(metadf)
    >>> print(distance_matrix.round(0))
              Brussels  Antwerp   Ghent
    Brussels       0.0  41200.0  52600.0
    Antwerp    41200.0      0.0  91200.0
    Ghent      52600.0  91200.0      0.0
    """
    # Validate input
    if lat_col not in metadf.columns or lon_col not in metadf.columns:
        raise ValueError(
            f"Columns '{lat_col}' and/or '{lon_col}' not found in DataFrame"
        )

    if metadf.empty:
        raise ValueError("Input DataFrame is empty")

    # Check for missing values
    if metadf[[lat_col, lon_col]].isnull().any().any():
        raise ValueError("Missing values found in latitude or longitude columns")

    # Extract coordinates and convert to radians for haversine
    coords = metadf[[lat_col, lon_col]].values
    coords_rad = np.radians(coords)

    # Create BallTree with haversine metric
    tree = BallTree(coords_rad, metric="haversine")

    # Calculate pairwise distances between all points
    n_stations = len(metadf)

    # Earth's radius in meters
    earth_radius_m = 6371000.0

    # Use tree to calculate distances from all points to all points
    distances_rad, indices = tree.query(coords_rad, k=n_stations, return_distance=True)

    # Convert distances from radians to meters
    distances_km = distances_rad * earth_radius_m

    # Reorder the distances matrix to match the original order
    # The BallTree returns results sorted by distance, so we need to reorder
    ordered_distances = np.zeros((n_stations, n_stations))
    for i in range(n_stations):
        for j in range(n_stations):
            # Find where station j appears in the results for station i
            pos = np.where(indices[i] == j)[0][0]
            ordered_distances[i, j] = distances_km[i, pos]

    # Create distance matrix DataFrame
    distance_matrix = pd.DataFrame(
        ordered_distances, index=metadf.index, columns=metadf.index
    )

    return distance_matrix
