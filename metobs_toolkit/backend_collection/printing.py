#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Printing Functions

@author: thoverga
"""


def dataset_string_repr(Dataset):

    return_str = ""

    def _data_details_as_str(Dataset_with_data):
        """Info on the dataset attributes"""
        n_stations = (
            Dataset_with_data.df.index.get_level_values("name").unique().shape[0]
        )
        obstypes = (
            Dataset_with_data.df.index.get_level_values("obstype").unique().to_list()
        )
        n_obs_tot = Dataset_with_data.df["value"].count()
        n_outl = Dataset_with_data.outliersdf.shape[0]
        startdt = Dataset_with_data.df.index.get_level_values("datetime").min()
        enddt = Dataset_with_data.df.index.get_level_values("datetime").max()

        detailstr = f"\
\n     *{n_stations} stations\
\n     *{obstypes} observation types present\
\n     *{n_obs_tot} observation records (not Nan's)\
\n     *{n_outl} records labeled as outliers\
\n     *{len(Dataset_with_data.gaps)} gaps\
\n     *records range: {startdt} --> {enddt} (total duration:  {enddt - startdt})\
\n     *time zone of the records: {str(Dataset_with_data._get_tz())}"

        if (not Dataset_with_data.metadf["lat"].isnull().all()) & (
            not Dataset_with_data.metadf["lon"].isnull().all()
        ):
            detailstr += "\n     *Coordinates are available for all stations."

        return detailstr

    def _metadata_details_as_str(Dataset_without_data):
        """Info on the metadata attributes."""
        n_stations = Dataset_without_data.metadf.shape[0]
        present_cols = Dataset_without_data.metadf.columns.to_list()

        detailstr = f"\
\n     *{n_stations} stations in the metadata\
\n     *The following columns are present in the metadf: {sorted(present_cols)}"

        if (not Dataset_without_data.metadf["lat"].isnull().all()) & (
            not Dataset_without_data.metadf["lon"].isnull().all()
        ):
            detailstr += "\n     *Coordinates are available for all stations."

        return detailstr

    def _geedata_details_as_str(Dataset_with_Gee):
        return f"\n     *Known GEE datasets for: {list(Dataset_with_Gee.gee_datasets.keys())}"

    if Dataset.df.empty:
        if Dataset.metadf.empty:
            # no data and no metadata
            return_str += f"Empty instance of {type(Dataset).__name__}"
            return return_str

        else:
            # no data but with metadata
            return_str += f"Instance of a {type(Dataset).__name__} (metadata-only)."
            return_str += _metadata_details_as_str(Dataset)
    else:
        # with data and metadata
        return_str += f"{type(Dataset).__name__} instance containing:"
        return_str += _data_details_as_str(Dataset)

    if (type(Dataset).__name__ == "Dataset") | (type(Dataset).__name__ == "Station"):
        return_str += _geedata_details_as_str(Dataset)

    return return_str


def print_dataset_info(dataset, show_all_settings=False, max_disp_n_gaps=5):
    """Print out settings of a dataset.

    Parameters
    ----------
    dataset : metobs_toolkit.Dataset
        The dataset to print the settings of.
    show_all_settings : bool, optional
        If True all settings are printed else a selection of the settings is
        printed. The default is False.
    max_disp_n_gaps : int, optional
        The maximum number of gaps to print detailed information of. The
        default is 5.

    Returns
    -------
    None.

    """
    print("--------  General ---------", "\n")
    print(dataset)

    print("\n", "--------  Observation types ---------", "\n")
    for obstype in dataset.obstypes.values():
        obstype.get_info()

    print("\n", "--------  Settings ---------", "\n")
    if show_all_settings:
        dataset.show_settings()
    else:
        print(
            "(to show all settings use the .show_settings() method, or set show_all_settings = True)"
        )

    print("\n", "--------  Outliers ---------", "\n")
    if dataset.outliersdf.empty:
        print("No outliers.")
    else:
        print(
            f"A total of {dataset.outliersdf.shape[0]} found with these occurrences: \n"
        )
        print(f'{dataset.outliersdf["label"].value_counts().to_dict()}')
    print("\n", "--------  Meta data ---------", "\n")
    if dataset.metadf.empty:
        print("No metadata is found.")
    else:
        relev_columns = []
        for col in dataset.metadf.columns:
            if not dataset.metadf[col].isna().all():
                relev_columns.append(col)

        print(f"The following metadata is found: {relev_columns}")
        print("\n The first rows of the metadf looks like:")
        print(f"{dataset.metadf[relev_columns].head()}")

    # "-------- Gaps  ---------")

    if dataset.gaps is not None:
        print("\n", "--------  Gaps ---------", "\n")
        if len(dataset.gaps) <= max_disp_n_gaps:
            i = 0
            while i < len(dataset.gaps):
                print(f"\n ** GAP at index [{i}] ** ")
                dataset.gaps[i].get_info()
                i += 1
        else:
            print(
                f"The info on {len(dataset.gaps)} is to long to print. Use the .get_gaps_fill_df() to get an overview DataFrame.."
            )

    print("\n", "--------  Known GEE Modeldata---------", "\n")
    if not bool(dataset.gee_datasets):
        print("No known GEE Modeldata.")
    else:
        for geemoddat in dataset.gee_datasets.values():
            print(f" * {geemoddat}")

        print(
            f'(For more details, use the .get_info method. Ex: dataset.gee_datasets["{next(iter(dataset.gee_datasets))}"].get_info() )'
        )


# =============================================================================
# Docstring test
# =============================================================================
if __name__ == "__main__":
    from metobs_toolkit.doctest_fmt import setup_and_run_doctest

    setup_and_run_doctest()
