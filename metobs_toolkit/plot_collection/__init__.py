# flake8: noqa: F401
from .default_style import default_plot_settings
from .general_functions import (
    create_axes,
    set_title,
    set_ylabel,
    set_xlabel,
    set_legend,
    create_categorical_color_map,
    format_datetime_axes,
)

from .qc_info_pies import qc_overview_pies

from .folium_spatial import (
    folium_map,
    add_stations_to_folium_map,
    add_title_to_folium_map,
)

from .timeseries_plotting import (
    plot_timeseries_color_by_label,
    plot_timeseries_color_by_station,
)

from .cycle_plotting import make_diurnal_plot
