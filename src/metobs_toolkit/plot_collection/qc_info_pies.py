import logging
import math
from typing import Tuple

import pandas as pd
import matplotlib.pyplot as plt

from metobs_toolkit.settings_collection.settings import Settings
from metobs_toolkit.plot_collection import create_axes

# Configure logging
from metobs_toolkit.backend_collection.decorators import log_entry

logger = logging.getLogger("<metobs_toolkit>")

def autopct_format(pct):
    return f'{pct:.1f}%' if pct > 0 else ''


@log_entry
def qc_overview_pies(
    end_labels_from_df: pd.Series,
    end_labels_from_outliers: pd.Series,
    per_check_labels: pd.Series,
    fig_title: str = ""
    
) -> plt.Figure:
    """
    Generate a QC overview figure with pie charts for label frequencies and per-check outcomes.

    Parameters
    ----------
    end_labels_from_df : pandas.Series
        Counts of final labels for all records (index as labels), typically from
        ``SensorData.df['label'].value_counts()``.
    end_labels_from_outliers : pandas.Series
        Counts limited to outlier records (index as labels), e.g. from
        ``SensorData.outliersdf['label'].value_counts()``. If empty, a single slice
        "No QC outliers" is drawn.
    per_check_labels : pandas.Series
        MultiIndex Series with index levels ``['checkname', 'flag']`` containing counts per
        QC check outcome (flags such as ``flagged_cond``, ``pass_cond``, ``unmet_cond``, etc.),
        as returned by ``SensorData.get_qc_freq_statistics()``.

    Returns
    -------
    matplotlib.figure.Figure
        Figure containing two large pies (all labels and outlier labels) and one small pie per
        QC check showing the distribution of its outcomes.
    """

    # Define layout
    ax = create_axes(**Settings.get("plotting_settings.pie_charts.figkwargs"))
    ax.set_axis_off()
    fig = ax.get_figure()

    ncol = Settings.get("plotting_settings.pie_charts.ncols")
    spec = fig.add_gridspec(4, ncol)
    ax_thl = fig.add_subplot(spec[0, :2])  # top half left
    ax_thr = fig.add_subplot(spec[0, 2:])  # top half right

    # Frequency with all
    colors = [Settings._get_color_from_label(label) for label in end_labels_from_df.index]
    end_labels_from_df.plot(
        ax=ax_thl,
        kind="pie",
        y="",
        autopct=autopct_format,
        legend=False,
        colors=colors,
        radius=Settings.get("plotting_settings.pie_charts.radius_big"),
        fontsize=Settings.get("plotting_settings.pie_charts.txt_label_size_big_pies"),
    )
    ax_thl.set_title("Label frequencies", **Settings.get("plotting_settings.pie_charts.big_pie_title_kwargs"))
    ax_thl.set_ylabel("")

    # Only outliers
   
    colors = [Settings._get_color_from_label(label) for label in end_labels_from_outliers.index]

    if end_labels_from_outliers.empty:
        # No outliers --> full pie with "No QC outliers" in the color of 'ok'
        end_labels_from_outliers = pd.Series([100], index=["No QC outliers"])
        colors = [Settings.get("label_def.goodrecord.color")]

    end_labels_from_outliers.plot(
        ax=ax_thr,
        kind="pie",
        y="",
        autopct=autopct_format,
        legend=False,
        colors=colors,
        radius=Settings.get("plotting_settings.pie_charts.radius_big"),
        fontsize=Settings.get("plotting_settings.pie_charts.txt_label_size_big_pies"),
    )
    ax_thr.set_title("Outlier specific frequencies", **Settings.get("plotting_settings.pie_charts.big_pie_title_kwargs"))
    ax_thr.set_ylabel("")

    # Performance per check
    per_qc_colmap = {val['label']: val['plotkwargs']['color'] for val in Settings.get('qc_status_labels_per_check').values()}

    
    i = 0
    for checkname in per_check_labels.index.get_level_values('checkname').unique():
        subax = fig.add_subplot(spec[math.floor(i / ncol) + 1, i % ncol])
        
        checkname_subset = per_check_labels.loc[checkname]
        colors = [per_qc_colmap.get(label, 'gray') for label in checkname_subset.index]
        
        checkname_subset.plot(
            ax=subax,
            kind="pie",
            autopct=autopct_format,
            legend=False,
            colors=colors,
            radius=Settings.get("plotting_settings.pie_charts.radius_small"),
            fontsize=Settings.get("plotting_settings.pie_charts.txt_label_size_small_pies"),
        )

        subax.set_title(f"{checkname}", **Settings.get("plotting_settings.pie_charts.small_pie_title_kwargs"))
        subax.set_ylabel("")

        i += 1
        
    
    fig.suptitle(fig_title, **Settings.get("plotting_settings.pie_charts.fig_title_kwargs"))
    return fig


# @log_entry
# def qc_overview_pies(
#     df: pd.DataFrame,
# ) -> plt.Figure:
#     """
#     Generate a quality control (QC) overview using pie charts.

#     Parameters
#     ----------
#     df : pandas.DataFrame
#         DataFrame containing QC data. Must include columns 'N_labeled', 'N_all', and 'N_checked'.

#     Returns
#     -------
#     matplotlib.figure.Figure
#         The generated figure containing the QC overview pie charts.

#     Raises
#     ------
#     TypeError
#         If any of the arguments are not of the expected type.
#     """

#     # Define layout
#     ax = create_axes(**Settings.get("plotting_settings.pie_charts.figkwargs"))
#     ax.set_axis_off()
#     fig = ax.get_figure()

#     ncol = Settings.get("plotting_settings.pie_charts.ncols")
#     spec = fig.add_gridspec(4, ncol)
#     ax_thl = fig.add_subplot(spec[0, :2])  # top half left
#     ax_thr = fig.add_subplot(spec[0, 2:])  # top half right

#     # Frequency with all
#     plotdf = df
#     colors = [Settings._get_color_from_label(label) for label in plotdf.index]
#     plotdf.plot(
#         ax=ax_thl,
#         kind="pie",
#         y="N_labeled",
#         autopct="%1.1f%%",
#         legend=False,
#         colors=colors,
#         radius=Settings.get("plotting_settings.pie_charts.radius_big"),
#         fontsize=Settings.get("plotting_settings.pie_charts.txt_size_big_pies"),
#     )
#     ax_thl.set_title("Label frequencies")
#     ax_thl.set_ylabel("")

#     # Outliers comparison
#     plotdf = df[
#         ~df.index.isin(
#             [
#                 Settings.get("label_def.goodrecord.label"),
#                 Settings.get("label_def.regular_gap.label"),
#             ]
#         )
#     ]

#     colors = [Settings._get_color_from_label(label) for label in plotdf.index]

#     if plotdf.empty:
#         # No outliers --> full pie with "No QC outliers" in the color of 'ok'
#         plotdf = pd.DataFrame(
#             data={"N_labeled": [100]}, index=pd.Index(data=["No QC outliers"])
#         )
#         colors = [Settings.get("label_def.goodrecord.color")]

#     plotdf.plot(
#         ax=ax_thr,
#         kind="pie",
#         y="N_labeled",
#         autopct="%1.1f%%",
#         legend=False,
#         colors=colors,
#         radius=Settings.get("plotting_settings.pie_charts.radius_big"),
#         fontsize=Settings.get("plotting_settings.pie_charts.txt_size_big_pies"),
#     )
#     ax_thr.set_title("Outlier specific frequencies")
#     ax_thr.set_ylabel("")

#     # Performance per check
#     plotdf = df[
#         ~df.index.isin(
#             [
#                 Settings.get("label_def.goodrecord.label"),
#                 Settings.get("label_def.regular_gap.label"),
#             ]
#         )
#     ]

#     # Label to QC check name map
#     label_too_qcname_map = Settings._label_to_qccheckmap()

#     i = 0
#     for idx, row in plotdf.iterrows():
#         # Target a specific axes
#         subax = fig.add_subplot(spec[math.floor(i / ncol) + 1, i % ncol])

#         # Construct a plot Series
#         plotseries = pd.Series(
#             {
#                 Settings.get("label_def.uncheckedrecord.label"): row["N_all"]
#                 - row["N_checked"],
#                 Settings.get("label_def.goodrecord.label"): row["N_checked"]
#                 - row["N_labeled"],
#                 Settings.get("label_def.outlier.label"): row["N_labeled"],
#             }
#         )
#         # Define colors
#         colors = [Settings._get_color_from_label(label) for label in plotseries.index]
#         plotseries.plot(
#             ax=subax,
#             kind="pie",
#             autopct="%1.1f%%",
#             legend=False,
#             colors=colors,
#             radius=Settings.get("plotting_settings.pie_charts.radius_small"),
#             fontsize=Settings.get("plotting_settings.pie_charts.txt_size_small_pies"),
#         )

#         subax.set_title(f"Effectiveness of {label_too_qcname_map[idx]}")
#         subax.set_ylabel("")

#         i += 1

#     logger.debug("Exiting qc_overview_pies function.")
#     return fig
