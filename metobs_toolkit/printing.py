#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Printing Functions

@author: thoverga
"""


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
