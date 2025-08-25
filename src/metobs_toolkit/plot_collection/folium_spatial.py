import logging
import sys
from typing import List, Dict, Optional

import pandas as pd
import matplotlib
import branca.colormap as brcm
import geemap.foliumap as geemap
import folium
from folium import plugins as folium_plugins

# Configure logging
from metobs_toolkit.backend_collection.loggingmodule import log_entry

logger = logging.getLogger("<metobs_toolkit>")


def _get_init_mapcenter(gdf: pd.DataFrame) -> List[float]:
    """
    Calculate the initial map center based on the GeoDataFrame.

    Parameters
    ----------
    gdf : pd.DataFrame
        GeoDataFrame containing the data.

    Returns
    -------
    List[float]
        A list containing the latitude and longitude of the map center.
    """
    logger.debug("Calculating initial map center.")
    centroid = gdf["geometry"].unary_union.centroid
    return [centroid.y, centroid.x]


@log_entry
def folium_map() -> geemap.Map:
    """
    Create a folium map using geemap.

    Returns
    -------
    geemap.Map
        A folium map object.
    """
    Map = geemap.Map(add_google_map=False)
    return Map


@log_entry
def add_title_to_folium_map(title: str, Map: folium.Map) -> folium.Map:
    """
    Add a title to a folium map.

    Parameters
    ----------
    title : str
        The title to add to the map.
    Map : folium.Map
        The folium map object.

    Returns
    -------
    folium.Map
        The updated folium map with the title added.
    """
    if not isinstance(title, str):
        raise TypeError("Argument 'title' must be of type str.")
    if not isinstance(Map, folium.Map):
        raise TypeError("Argument 'Map' must be of type folium.Map.")

    title_html = """
                 <h3 align="center" style="font-size:20px"><b>{}</b></h3>
                 """.format(
        title
    )

    Map.get_root().html.add_child(folium.Element(title_html))
    return Map


@log_entry
def add_stations_to_folium_map(
    Map: folium.Map, metadf: pd.DataFrame, display_cols: List[str] = ["name"]
) -> folium.Map:
    """
    Add station markers to a folium map.

    Parameters
    ----------
    Map : folium.Map
        The folium map object.
    metadf : pd.DataFrame
        A DataFrame containing station metadata.
    display_cols : List[str], optional
        Columns to display in the popup, by default ["name"].

    Returns
    -------
    folium.Map
        The updated folium map with station markers added.
    """

    metadf = metadf.reset_index()
    metadf["geometry"] = metadf["geometry"].to_crs("epsg:4326")
    for _, row in metadf.iterrows():
        point = row["geometry"]
        popuptext = ""
        for disp in display_cols:
            popuptext = f"{popuptext}{row[disp]}\n"
        folium.Marker(
            location=[point.y, point.x], fill_color="#43d9de", popup=popuptext, radius=8
        ).add_to(Map)

    return Map


@log_entry
def make_folium_html_plot(
    gdf: pd.DataFrame,
    variable_column: str,
    var_display_name: str,
    var_unit: str,
    label_column: str,
    label_col_map: Dict[str, str],
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    radius: int = 13,
    fill_alpha: float = 0.6,
    mpl_cmap_name: str = "viridis",
    max_fps: int = 4,
    dt_disp_fmt: str = "%Y-%m-%d %H:%M",
) -> folium.Map:
    """
    Create an interactive folium map with time-stamped data.

    Parameters
    ----------
    gdf : pd.DataFrame
        GeoDataFrame containing the data to plot.
    variable_column : str
        Column name for the variable to visualize.
    var_display_name : str
        Display name for the variable.
    var_unit : str
        Unit of the variable.
    label_column : str
        Column name for the labels.
    label_col_map : Dict[str, str]
        Mapping of labels to colors.
    vmin : Optional[float], optional
        Minimum value for the colormap, by default None.
    vmax : Optional[float], optional
        Maximum value for the colormap, by default None.
    radius : int, optional
        Radius of the markers, by default 13.
    fill_alpha : float, optional
        Opacity of the marker fill, by default 0.6.
    mpl_cmap_name : str, optional
        Name of the matplotlib colormap, by default "viridis".
    max_fps : int, optional
        Maximum frames per second for the animation, by default 4.
    dt_disp_fmt : str, optional
        Date-time display format, by default "%Y-%m-%d %H:%M".

    Returns
    -------
    folium.Map
        The generated folium map.
    """

    # Create a map
    m = folium.Map(
        location=_get_init_mapcenter(gdf),
        tiles="cartodbpositron",
        zoom_start=10,
        attr="<a href=https://github.com/vergauwenthomas/MetObs_toolkit </a>",
    )

    # Add extra tiles
    folium.TileLayer("OpenStreetMap", overlay=False, name="OSM").add_to(m)

    # Determine colormap bounds
    if vmin is None:
        vmin = gdf[variable_column].min()
    if vmax is None:
        vmax = gdf[variable_column].max()

    # Create colormap
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
    mapper = matplotlib.cm.ScalarMappable(
        norm=norm, cmap=matplotlib.colormaps[mpl_cmap_name]
    )
    colormap = brcm.LinearColormap(
        colors=mapper.cmap.colors,
        index=None,
        vmin=vmin,
        vmax=vmax,
        caption=f"{var_display_name} ({var_unit}) colorbar",
    )

    # Map values to colors
    @log_entry
    def map_value_to_hex(series, vmin, vmax, cmapname="viridis"):
        norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
        mapper = matplotlib.cm.ScalarMappable(
            norm=norm, cmap=matplotlib.colormaps[cmapname]
        )
        return series.apply(lambda x: str(matplotlib.colors.to_hex(mapper.to_rgba(x))))

    gdf["value_color"] = map_value_to_hex(
        gdf[variable_column], vmin, vmax, cmapname=mpl_cmap_name
    )

    # Check for unmapped labels
    unmapped_labels = [
        lab for lab in gdf[label_column].unique() if lab not in label_col_map.keys()
    ]
    if len(unmapped_labels) > 0:
        logger.error(f"Unmapped labels found: {unmapped_labels}")
        sys.exit(f"Unmapped labels found: {unmapped_labels}")

    gdf["label_color"] = gdf[label_column].map(label_col_map)

    # Serialize data to features
    @log_entry
    def make_scatter_feature(row):
        dtstring = pd.to_datetime([row["datetime"]]).strftime(dt_disp_fmt)[0]
        coords = [[row["geometry"].x, row["geometry"].y]]
        popup_str = f" <b>{row['name']}</b>  <br> {'{:.1f}'.format(row[variable_column])} {var_unit} <br> {row[label_column]}"

        features_instance = {
            "type": "Feature",
            "geometry": {
                "type": "MultiPoint",
                "coordinates": coords,
            },
            "properties": {
                "times": [dtstring],
                "popup": popup_str,
                "tooltip": f'{row["name"]}',
                "id": "geenidee",  # TYPO
                "icon": "circle",
                "iconstyle": {
                    "fillColor": row["value_color"],
                    "fillOpacity": fill_alpha,
                    "stroke": "false",
                    "radius": radius,
                    "color": row["label_color"],
                },
            },
        }
        return features_instance

    features = gdf.apply(make_scatter_feature, axis=1).to_list()

    # Add data to the map
    folium_plugins.TimestampedGeoJson(
        {
            "type": "FeatureCollection",
            "features": features,
        },
        period="PT1H",
        duration="PT1H",
        add_last_point=False,
        auto_play=False,
        loop=False,
        max_speed=max_fps,
        loop_button=True,
        date_options="YYYY/MM/DD HH:mm:ss",
        time_slider_drag_update=True,
    ).add_to(m)

    m.add_child(colormap)
    folium.LayerControl().add_to(m)

    return m
