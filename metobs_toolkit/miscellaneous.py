#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A collection of miscellanious functions, for (automatic) documentation purposes.
This module is not intended for users.

Created on Mon Jul 15 11:34:07 2024

@author: thoverga
"""

from metobs_toolkit.obstypes import Obstype, tlk_obstypes


# =============================================================================
# A function to print out all standard toolkit obstypes info
# =============================================================================


def _tlk_print_standard_obstypes():
    print("The standard observations present in the Metobs toolkit")
    print(" ----------------------------------------------------- \n")
    ncol = 24
    for std_obs in tlk_obstypes.values():
        print(
            f"{std_obs.name.ljust(ncol)} | {std_obs.get_description().ljust(ncol)} | {std_obs.get_standard_unit().ljust(ncol)} "
        )


# =============================================================================
# Docstring test
# =============================================================================
if __name__ == "__main__":
    from metobs_toolkit.doctest_fmt import setup_and_run_doctest

    setup_and_run_doctest()
