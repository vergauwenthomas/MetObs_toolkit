from __future__ import annotations
import os
import logging
import concurrent.futures
from typing import Union, List, Dict, Tuple, TYPE_CHECKING, Optional

from metobs_toolkit.qcresult import (
    unchecked_cond,
    unmet_cond,
    pass_cond,
    flagged_cond,
    saved_cond,
)

import numpy as np
import pandas as pd
# Constants for buddy check status labels
BC_NOT_TESTED = "not_tested"  # Value was NaN, not tested
BC_NO_BUDDIES = "no_buddies"  # Not enough buddies to test
BC_PASSED = "passed"  # Tested and passed
BC_FLAGGED = "flagged"  # Tested and flagged as outlier
BC_SAFETYNET_SAVED = "safetynet_saved"  # Flagged but saved by safetynet
BC_SAFETYNET_OUTLIER = "safetynet_outlier"  # Flagged but not saved by safetynet
BC_WHITELIST_SAVED = "whitelist_saved"  # Flagged but saved by whitelist
BC_WHITELIST_NOT_SAVED = "whitelist_not_saved"  # Flagged but not saved by whitelist
BC_CHECK_SKIPPED = "skipped" #This check was skipped, e.g. due to arugments of the user (not whitelist, not safetynets etc)

if TYPE_CHECKING:
    from metobs_toolkit.station import Station

logger = logging.getLogger("<metobs_toolkit>")

class BuddyCheckStation:
    """Wrapper for a Station with buddy check-specific details.
    
    This class wraps a Station object and adds information about how it is
    handled during the buddy check process, including buddy assignment,
    filtering steps, and participation in buddy groups.
    
    Attributes
    ----------
    station : Station
        The wrapped Station object.
    _buddy_groups : dict
        Dictionary mapping group names to lists of buddy station names.
    flag_lapsrate_corrections : bool
        Whether lapse rate corrections have been applied.
    cor_term : float
        The correction term applied for lapse rate.
    flags : pandas.DataFrame
        DataFrame with MultiIndex (datetime, iteration) containing flag values.
        Columns are added via `add_flags` method for different check types.
    details : dict
        Dictionary storing iteration-wise detail information. Structure:
        {
            'spatial_check': {
                iteration_int: Series(index=DatetimeIndex, data=detail_strings),
                ...
            },
            'safetynet_check': {
                groupname_str: {
                    iteration_int: Series(index=DatetimeIndex, data=detail_strings),
                    ...
                },
                ...
            },
            'whitelist_check': {
                iteration_int: Series(index=DatetimeIndex, data=detail_strings),
                ...
            },
        }
    """
    station: Station
    
    def __init__(self, station: Station):
        self.station = station
        # Initialize instance-specific attributes (NOT class attributes!)
        self._buddy_groups: Dict[str, List[str]] = {
            'spatial': [],
        }
        
        # Value corrections
        self.flag_lapsrate_corrections: bool = False
        self.cor_term: float = 0.
        
        # Flags DataFrame with MultiIndex (datetime, iteration)
        self._flags: pd.DataFrame = pd.DataFrame()
        
        # Details dictionary structure
        self.details: Dict[str, Union[Dict[int, pd.Series], Dict[str, Dict[int, pd.Series]]]] = {
            'spatial_check': {},
            'safetynet_check': {},  # Dict of groupname -> Dict of iteration -> Series
            'whitelist_check': {},
        }
        
    @property
    def name(self) -> str:
        """Get the station name."""
        return self.station.name
    
    @property
    def flags(self) -> pd.DataFrame:
        """Get the flags DataFrame."""
        if self._flags.empty:
            return pd.DataFrame(index=pd.MultiIndex(levels=[[], []], codes=[[], []], names=['datetime', 'iteration']))
        
        return self._flags
    
    @flags.setter
    def flags(self, flags: pd.DataFrame) -> None:
        
        if not isinstance(flags, pd.DataFrame):
            raise ValueError("flags must be a pandas DataFrame")
        
        if not flags.empty:
            if not isinstance(flags.index, pd.MultiIndex):
                raise ValueError("flags DataFrame must have a MultiIndex")
            if flags.index.names != ['datetime', 'iteration']:
                raise ValueError("flags DataFrame MultiIndex must have levels ['datetime', 'iteration']")
        
        # Preserve column order: existing columns first, new columns at the end
        if not self._flags.empty and not flags.empty:
            existing_cols = [col for col in self._flags.columns if col in flags.columns]
            new_cols = [col for col in flags.columns if col not in self._flags.columns]
            ordered_cols = existing_cols + new_cols
            flags = flags[ordered_cols]
        
        self._flags = flags
    
   
    
    def add_flags(self, iteration: int, flag_series: pd.Series, column_name: str) -> None:
        """Add flags to the flags DataFrame for a specific iteration.
        
        Parameters
        ----------
        iteration : int
            The iteration number.
        flag_series : pd.Series
            Series with DatetimeIndex containing flag values.
        column_name : str
            The name of the column to add/update (e.g., 'spatial_check', 
            'safetynet_check:groupname', 'whitelist_check').
        """
        if flag_series.empty:
            return
            
        # Remove duplicates (keep first occurrence)
        flag_series = flag_series[~flag_series.index.duplicated(keep='first')]
        
        # Create a DataFrame with MultiIndex for the new flags
        new_flags = pd.DataFrame({
            column_name: flag_series.values
        }, index=pd.MultiIndex.from_arrays(
            [flag_series.index, [iteration] * len(flag_series)],
            names=['datetime', 'iteration']
        ))
        
        if self.flags.empty:
            self.flags = new_flags
        else:
            # Merge new flags with existing flags
            # Use combine_first to keep existing values and add new ones
            self.flags = self.flags.combine_first(new_flags)
            
            # If the column already exists, update with new values (not NaN)
            if column_name in self.flags.columns:
                # For indices that exist in both, update from new_flags
                common_idx = self.flags.index.intersection(new_flags.index)
                if not common_idx.empty:
                    self.flags.loc[common_idx, column_name] = new_flags.loc[common_idx, column_name]
            
        self.flags = self.flags.sort_index()
    
    def filter_buddies(self, filteredbuddies: List[str], groupname: str) -> None:
        self.set_buddies(filteredbuddies, groupname=f'{groupname}_filtered')
    
    def set_buddies(self, buddies: List[str], groupname: str) -> None:
        self._buddy_groups.update({groupname: buddies})
    
    def get_buddies(self, groupname: str) -> List[str]:
        if f'{groupname}_filtered' in self._buddy_groups:
            return self._buddy_groups[f'{groupname}_filtered']
        if groupname in self._buddy_groups:
            return self._buddy_groups[groupname]
        raise ValueError(f"Unknown buddy group: {groupname}")
        
    def has_enough_buddies(self, groupname: str, min_buddies: int) -> bool:
        """Check if the station has enough final buddies."""
        enough = len(self.get_buddies(groupname=groupname)) >= min_buddies
        return enough
    
    
    def _update_details(self, iteration: int, detail_series: pd.Series, groupname: str) -> None:
        if detail_series.empty:
            return
        # Remove duplicates (keep first occurrence)
        detail_series = detail_series[~detail_series.index.duplicated(keep='first')]
        
        # Store details in the details dictionary
        if iteration not in self.details[groupname]:
            self.details[groupname][iteration] = detail_series
        else:
            #FIXME: i do not think this branch is ever used
            # Append to existing series for this iteration
            existing = self.details[groupname][iteration]
            combined = pd.concat([existing, detail_series])
            # Remove duplicates keeping first
            self.details[groupname][iteration] = combined[~combined.index.duplicated(keep='first')].sort_index()
    
    def add_spatial_details(self, iteration: int, detail_series: pd.Series) -> None:
        """Add spatial check detail information for an iteration.
        
        Parameters
        ----------
        iteration : int
            The iteration number.
        detail_series : pd.Series
            Series with DatetimeIndex containing detail messages.
        """
        self._update_details(
            iteration=iteration,
            detail_series=detail_series,
            groupname='spatial_check'
        )
    
    def update_safetynet_details(self,
                                detailseries: pd.Series,
                                iteration: int,
                                groupname: str,
                                is_saved: bool = True,) -> None:
        """Update safetynet check saved information for an iteration.
        
        Parameters
        ----------
        detailseries : pd.Series
            Series with DatetimeIndex containing detail messages for saved records.
        iteration : int
            The iteration number.
        groupname : str
            The name of the safetynet group that saved the records.
        is_saved : bool, optional
            Whether the records were saved by the safetynet (True) or failed (False).
            Default is True.
        """
        # Remove duplicates (keep first occurrence)
        detailseries = detailseries[~detailseries.index.duplicated(keep='first')]
        
        if detailseries.empty:
            return
        
        if is_saved:
            flag = BC_SAFETYNET_SAVED
        else:   
            flag = BC_SAFETYNET_OUTLIER
            
        # Add flags to the flags DataFrame
        column_name = f'safetynet_check:{groupname}'
        flag_series = pd.Series(flag, index=detailseries.index)
        self.add_flags(iteration=iteration, flag_series=flag_series, column_name=column_name)
        
        # Store details in the details dictionary
        # Ensure the groupname entry exists
        if groupname not in self.details['safetynet_check']:
            self.details['safetynet_check'][groupname] = {}
        
        if iteration not in self.details['safetynet_check'][groupname]:
            self.details['safetynet_check'][groupname][iteration] = detailseries
        else:
            # Append to existing series for this iteration
            existing = self.details['safetynet_check'][groupname][iteration]
            combined = pd.concat([existing, detailseries])
            # Remove duplicates keeping first
            self.details['safetynet_check'][groupname][iteration] = combined[~combined.index.duplicated(keep='first')].sort_index()
    
    def update_whitelist_details(self, whitelistseries: pd.Series,
                                 iteration: int,
                                 is_saved: bool = True) -> None:
        """Update whitelist check saved information for an iteration.
        
        Parameters
        ----------
        whitelistseries : pd.Series
            Series with DatetimeIndex containing detail messages for whitelisted records.
        iteration : int
            The iteration number.
        is_saved : bool, optional
            Whether the records were saved by the whitelist (True) or not (False).
            Default is True.
        """
        # Remove duplicates (keep first occurrence)
        whitelistseries = whitelistseries[~whitelistseries.index.duplicated(keep='first')]
        
        if whitelistseries.empty:
            return
        
        if is_saved:
            flag = BC_WHITELIST_SAVED
        else:   
            flag = BC_WHITELIST_NOT_SAVED
        # Add flags to the flags DataFrame
        flag_series = pd.Series(flag, index=whitelistseries.index)
        self.add_flags(iteration=iteration, flag_series=flag_series, column_name='whitelist_check')
        
        # Store details in the details dictionary
        if iteration not in self.details['whitelist_check']:
            self.details['whitelist_check'][iteration] = whitelistseries
        else:
            # Append to existing series for this iteration
            existing = self.details['whitelist_check'][iteration]
            combined = pd.concat([existing, whitelistseries])
            # Remove duplicates keeping first
            self.details['whitelist_check'][iteration] = combined[~combined.index.duplicated(keep='first')].sort_index()
    
    
    
    def _get_iterations(self) -> List[int]:
        """Get all iterations that have been processed."""
        iterations = set()
        
        # From spatial_check
        iterations.update(self.details['spatial_check'].keys())
        
        # From safetynet_check (nested)
        for groupname, iter_dict in self.details['safetynet_check'].items():
            iterations.update(iter_dict.keys())
        
        # From whitelist_check
        iterations.update(self.details['whitelist_check'].keys())
        
        return sorted(iterations)
    
    def get_final_labels(self) -> pd.Series:
        flags = self.flags
        
        #if no whitelist check has been performed, create an empyt column
        if 'whitelist_check' not in flags.columns:
            flags['whitelist_check'] = BC_CHECK_SKIPPED #'skipped'
        
        final_labels = flags.groupby('datetime').apply(final_label_logic)
        final_labels.name = 'final_label'
        return final_labels
    
    def get_final_details(self) -> pd.Series:
        """Get detailed description strings for each timestamp based on final label logic.
        
        This method returns a Series with detailed description strings that illustrate
        how the final label was determined. The details are extracted from the details
        attribute and combined based on the check pipeline.
        
        Returns
        -------
        pd.Series
            Series with DatetimeIndex containing detailed description strings.
            The series name is 'final_details'.
        """
        final_details = self.flags.groupby('datetime').apply(
            lambda subset: final_detail_logic(subset, self.details)
        )
        final_details.name = 'final_details'
        return final_details
    
    
    def get_iteration_summary(self, iteration: int) -> Dict[str, int]:
        """Get a summary of record counts for a specific iteration.
        
        Parameters
        ----------
        iteration : int
            The iteration number to summarize.
            
        Returns
        -------
        dict
            Dictionary with check types as keys and record counts as values.
            For safetynet_check, includes both total and per-group counts.
        """
        # Count spatial_check
        spatial_count = 0
        if iteration in self.details['spatial_check']:
            spatial_count = len(self.details['spatial_check'][iteration])
        
        # Count safetynet_check per group
        safetynet_per_group = {}
        for groupname, iter_dict in self.details['safetynet_check'].items():
            if iteration in iter_dict:
                safetynet_per_group[groupname] = len(iter_dict[iteration])
        safetynet_total = sum(safetynet_per_group.values())
        
        # Count whitelist_check
        whitelist_count = 0
        if iteration in self.details['whitelist_check']:
            whitelist_count = len(self.details['whitelist_check'][iteration])
        
        return {
            'spatial_check': spatial_count,
            'safetynet_check': safetynet_per_group,
            'safetynet_check_total': safetynet_total,
            'whitelist_check': whitelist_count
        }
    
    def get_info(self) -> str:
        """Get a summary of the BuddyCheckStation status and attributes.
        
        Returns
        -------
        str
            Formatted string with overview of the station's buddy check status.
        """
        lines = []
        lines.append("=" * 60)
        lines.append(f"BuddyCheckStation: {self.name}")
        lines.append("=" * 60)
        
        # Buddy groups
        lines.append("\n--- Buddy Groups ---")
        if self._buddy_groups:
            for groupname, buddies in self._buddy_groups.items():
                n_buddies = len(buddies)
                buddies_str = ", ".join(buddies[:5])
                if n_buddies > 5:
                    buddies_str += f", ... (+{n_buddies - 5} more)"
                lines.append(f"  {groupname}: {n_buddies} buddies")
                if buddies:
                    lines.append(f"    [{buddies_str}]")
        else:
            lines.append("  No buddy groups assigned")
        
        # Corrections
        lines.append("\n--- Value Corrections ---")
        lines.append(f"  Lapse rate correction: {self.flag_lapsrate_corrections}")
        if self.flag_lapsrate_corrections:
            lines.append(f"  Correction term: {self.cor_term:.4f}")
        
        # Flags summary
        lines.append("\n--- Flags ---")
        if not self.flags.empty:
            lines.append(f"  Total flag entries: {len(self.flags)}")
            lines.append(f"  Flag columns: {list(self.flags.columns)}")
        else:
            lines.append("  No flags recorded")
        
        # Iteration status
        iterations = self._get_iterations()
        lines.append("\n--- Iteration Status ---")
        if iterations:
            lines.append(f"  Iterations processed: {len(iterations)}")
            
            # Totals across all iterations
            total_spatial = 0
            total_safetynet = 0
            total_whitelist = 0
            safetynet_groups_total: Dict[str, int] = {}
            
            for iteration in iterations:
                summary = self.get_iteration_summary(iteration)
                total_spatial += summary['spatial_check']
                total_safetynet += summary['safetynet_check_total']
                total_whitelist += summary['whitelist_check']
                
                # Accumulate per-group safetynet totals
                for groupname, count in summary['safetynet_check'].items():
                    safetynet_groups_total[groupname] = safetynet_groups_total.get(groupname, 0) + count
                
                lines.append(f"\n  Iteration {iteration}:")
                lines.append(f"    Spatial outliers: {summary['spatial_check']}")
                if summary['safetynet_check']:
                    lines.append(f"    Safetynet saved (total): {summary['safetynet_check_total']}")
                    for groupname, count in summary['safetynet_check'].items():
                        lines.append(f"      - {groupname}: {count}")
                else:
                    lines.append(f"    Safetynet saved: 0")
                lines.append(f"    Whitelist saved: {summary['whitelist_check']}")
            
            lines.append(f"\n  --- Totals ---")
            lines.append(f"  Total spatial outliers: {total_spatial}")
            lines.append(f"  Total safetynet saved: {total_safetynet}")
            if safetynet_groups_total:
                for groupname, count in safetynet_groups_total.items():
                    lines.append(f"    - {groupname}: {count}")
            lines.append(f"  Total whitelist saved: {total_whitelist}")
        else:
            lines.append("  No iterations processed")
        
        lines.append("=" * 60)
        
        info_str = "\n".join(lines)
        print(info_str)
        return info_str
    
    def map_timestamps(self,
                       timestamp_map: Dict[str, pd.Series]) -> None:
        """Map synchronized timestamps to original timestamps for this station.
        
        This function maps the synchronized timestamps (used during buddy check 
        processing) back to the original timestamps for this station's flags and 
        details attributes.
        
        Parameters
        ----------
        timestamp_map : dict
            Dictionary mapping station names to Series where index is synchronized
            timestamp and value is original timestamp.
            
        Returns
        -------
        None
            Modifies the station in-place.
        """
        
        # ts_map = timestamp_map[self.name]
        
        # Revert flags DataFrame timestamps (MultiIndex: datetime, iteration)
        if not self.flags.empty:
            self.flags = _map_dt_index(
                pdobj=self.flags,
                ts_map=timestamp_map,
                datetime_level='datetime'
            )
        
        # Revert details timestamps
        # spatial_check: {iteration: Series}
        for iteration, detail_series in self.details['spatial_check'].items():
            self.details['spatial_check'][iteration] =  _map_dt_index(
                pdobj=detail_series,
                ts_map=timestamp_map
            )
        
        # safetynet_check: {groupname: {iteration: Series}}
        for groupname, iter_dict in self.details['safetynet_check'].items():
            for iteration, detail_series in iter_dict.items():
                self.details['safetynet_check'][groupname][iteration] =  _map_dt_index(
                    pdobj=detail_series,
                    ts_map=timestamp_map
                )
        
        # whitelist_check: {iteration: Series}
        for iteration, detail_series in self.details['whitelist_check'].items():
            self.details['whitelist_check'][iteration] = _map_dt_index(
                pdobj=detail_series,
                ts_map=timestamp_map
            )

  
def final_label_logic(subset: pd.DataFrame) -> str:
    #the flag not tested is present on ALL iterations !
    if subset['spatial_check'].apply(lambda x: x=='not_tested').all():
        return BC_NOT_TESTED

    # --- passed condition ---- 
    #1a perfect pass (pass on last iteration of spatial check)
    if subset['spatial_check'].iloc[-1] == BC_PASSED:
        return BC_PASSED
    
    if subset['spatial_check'].iloc[-1] == BC_NO_BUDDIES:
        #Choice made: it can happen that a record passed a previous iteration,
        #but has not enough buddies in the last iteration. Applying the 'save' logic,
        #it is best to label this as no_buddies
        return BC_NO_BUDDIES

    #catched by safetynet
    #if there is at least (there can only be one), pass in the last iteration of saftynets
    saftynet_cols = [col for col in subset.columns if col.startswith('safetynet_check:')]
    if any(subset.iloc[-1][saftynet_cols] == BC_PASSED): #TODO not shure of this string
        return BC_SAFETYNET_SAVED
        

    #catched by whitelist
    if subset['whitelist_check'].iloc[-1] == BC_WHITELIST_SAVED:
        return BC_WHITELIST_SAVED
        
    # --- failed condition ----
    
    #fail in last iteration of spatial check
    if ((subset['spatial_check'].iloc[-1] == BC_FLAGGED) and
        all(subset.iloc[-1][saftynet_cols] != BC_PASSED) and #not passed is [nan, flagged, no-buddies]
        (subset['whitelist_check'].iloc[-1] != BC_WHITELIST_SAVED)): #not saved is [nan, skipped, not-saved]
        return BC_FLAGGED
    
    #fail in any previous iteration of spatial check
    if ((any(subset['spatial_check'] == BC_FLAGGED)) and
        (subset['spatial_check'].iloc[-1] == BC_NOT_TESTED)):
        return BC_FLAGGED
    
    
    raise ValueError(f"Unforeseen situartion encountered in final label logic: \n {subset}")
                


def final_detail_logic(subset: pd.DataFrame, details: Dict) -> str:
    """Extract detailed description string based on the final label logic.
    
    This function mirrors the logic of `final_label_logic` but returns
    a detailed description string extracted from the details dictionary.
    
    Parameters
    ----------
    subset : pd.DataFrame
        DataFrame subset for a single timestamp with all iterations.
    details : dict
        The details dictionary from BuddyCheckStation containing:
        - 'spatial_check': {iteration: Series}
        - 'safetynet_check': {groupname: {iteration: Series}}
        - 'whitelist_check': {iteration: Series}
        
    Returns
    -------
    str
        Detailed description string for the final label.
    """
    # Get the timestamp from the subset index
    timestamp = subset.index.get_level_values('datetime')[0]
    last_iteration = subset.index.get_level_values('iteration')[-1]
    
    # Helper to get detail from a specific check type and iteration
    def get_spatial_detail(iteration: int) -> str:
        if iteration in details['spatial_check']:
            detail_series = details['spatial_check'][iteration]
            if timestamp in detail_series.index:
                return str(detail_series.loc[timestamp])
        return ""
    
    def get_safetynet_detail(groupname: str, iteration: int) -> str:
        if groupname in details['safetynet_check']:
            if iteration in details['safetynet_check'][groupname]:
                detail_series = details['safetynet_check'][groupname][iteration]
                if timestamp in detail_series.index:
                    return str(detail_series.loc[timestamp])
        return ""
    
    def get_whitelist_detail(iteration: int) -> str:
        if iteration in details['whitelist_check']:
            detail_series = details['whitelist_check'][iteration]
            if timestamp in detail_series.index:
                return str(detail_series.loc[timestamp])
        return ""
    
    # Apply same logic as final_label_logic to determine which detail to return
    
    # Not tested condition
    if subset['spatial_check'].apply(lambda x: x == 'not_tested').all():
        return "Value was NaN, not tested."
    
    # Passed condition - pass on last iteration of spatial check
    if subset['spatial_check'].iloc[-1] == BC_PASSED:
        detail = get_spatial_detail(last_iteration)
        return detail if detail else "Passed spatial check."
    
    # No buddies condition
    if subset['spatial_check'].iloc[-1] == BC_NO_BUDDIES:
        detail = get_spatial_detail(last_iteration)
        return detail if detail else "Not enough buddies to test."
    
    # Caught by safetynet
    saftynet_cols = [col for col in subset.columns if col.startswith('safetynet_check:')]
    for col in saftynet_cols:
        if subset.iloc[-1][col] == BC_PASSED:
            groupname = col.split(':', 1)[1]
            detail = get_safetynet_detail(groupname, last_iteration)
            return detail if detail else f"Saved by safetynet ({groupname})."
    
    # Caught by whitelist
    if subset['whitelist_check'].iloc[-1] == BC_WHITELIST_SAVED:
        detail = get_whitelist_detail(last_iteration)
        return detail if detail else "Saved by whitelist."
    
    # Failed conditions
    
    # Fail in last iteration of spatial check
    if ((subset['spatial_check'].iloc[-1] == BC_FLAGGED) and
        all(subset.iloc[-1][saftynet_cols] != BC_PASSED) and
        (subset['whitelist_check'].iloc[-1] == BC_WHITELIST_NOT_SAVED)):
        # Combine details from spatial check and whitelist
        spatial_detail = get_spatial_detail(last_iteration)
        whitelist_detail = get_whitelist_detail(last_iteration)
        parts = [p for p in [spatial_detail, whitelist_detail] if p]
        return " | ".join(parts) if parts else "Flagged as outlier."
    
    # Fail in any previous iteration of spatial check
    if ((any(subset['spatial_check'] == BC_FLAGGED)) and
        (subset['spatial_check'].iloc[-1] == BC_NOT_TESTED)):
        # Find the iteration where it was flagged
        for idx, row in subset.iterrows():
            if row['spatial_check'] == BC_FLAGGED:
                flagged_iteration = idx[1]  # iteration from MultiIndex
                detail = get_spatial_detail(flagged_iteration)
                return detail if detail else f"Flagged in iteration {flagged_iteration}."
    
    return "Unknown condition."


def _map_dt_index(pdobj: pd.Series | pd.DataFrame,
    ts_map: pd.Series,
    datetime_level: str = 'datetime'
) -> pd.DataFrame:
    """Revert timestamps in a DataFrame with MultiIndex containing datetime level.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with MultiIndex containing a datetime level.
    ts_map : pd.Series
        Series mapping synchronized timestamps (index) to original timestamps (values).
    datetime_level : str
        Name of the datetime level in the MultiIndex.
        
    Returns
    -------
    pd.DataFrame
        DataFrame with reverted timestamps in the MultiIndex.
    """
    
    if isinstance(pdobj, pd.Series):
        df = pdobj.to_frame()
        returnseries = True
    else:
        df = pdobj
        returnseries = False
    
    # Get the current index
    old_index = df.index
    level_names = old_index.names
    
    df = df.reset_index()
    df['_mapped_datetime'] = df[datetime_level].map(lambda x: ts_map.get(x, x) if pd.notna(ts_map.get(x, pd.NaT)) else x)
    df = df.drop(columns=[datetime_level])
    df = df.rename(columns={'_mapped_datetime': datetime_level})
    df = df.set_index(level_names)
    df = df.sort_index()
    if returnseries:
        return df.iloc[:,0]
    return df



to_qc_labels_map = {
    BC_NOT_TESTED: unchecked_cond,  # Value was NaN, not tested
    BC_NO_BUDDIES: unmet_cond,  # Not enough buddies to test
    BC_PASSED : pass_cond, # Tested and passed
    BC_FLAGGED : flagged_cond,  # Tested and flagged as outlier
    BC_SAFETYNET_SAVED : pass_cond, # IMPORTANT !!!
    # BC_SAFETYNET_OUTLIER : flagged_cond  # Flagged but not saved by safetynet
    BC_WHITELIST_SAVED : saved_cond,  # Flagged but saved by whitelist
    # BC_WHITELIST_NOT_SAVED : flagged_cond
}