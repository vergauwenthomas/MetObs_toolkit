#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 14:56:26 2023

@author: thoverga
"""





def print_dataset_info(dataset, show_all_settings=False, max_disp_n_gaps = 5):

    print("\n", "--------  General ---------", "\n")
    print(dataset)

    print("\n", "--------  Settings ---------", "\n")
    if show_all_settings:
        dataset.show_settings()
    else:
        print('(to show all settings use the .show_settings() method, or set show_all_settings = True)')

    print("\n", "--------  Meta data ---------", "\n")
    if dataset.metadf.empty:
        print('No metadata is found.')
    else:
        relev_columns = []
        for col in dataset.metadf.columns:
            if not dataset.metadf[col].isna().all():
                relev_columns.append(col)

        print(f'The following metadata is found: {relev_columns}')
        print('\n The first rows of the metadf looks like:')
        print(f'{dataset.metadf[relev_columns].head()}')

    # "--------  Missing observations ---------")
    if not dataset.missing_obs is None:
        print(dataset.get_missing_obs_info())

    if not dataset.gaps is None:
        print("\n", "--------  Gaps ---------", "\n")
        if len(dataset.gaps) <=  max_disp_n_gaps:
            print(dataset.get_gaps_info())
        else:
            print(f'The info on {len(dataset.gaps)} is to long to print. Use the .get_gaps_info() to print out the details of all gaps.')