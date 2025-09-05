import logging
from typing import Literal, Union
import pandas as pd
import numpy as np

from metobs_toolkit.backend_collection.errorclasses import (
    MetObsTimeSimplifyError,
)
from metobs_toolkit.backend_collection.df_helpers import to_timedelta

from metobs_toolkit.backend_collection.loggingmodule import log_entry

logger = logging.getLogger("<metobs_toolkit>")


class TimestampMatcher:
    """
    Class for mapping and resampling timestamp records.

    Parameters
    ----------
    orig_records : pd.Series
        Original timestamp records.
    """

    def __init__(self, orig_records: pd.Series) -> None:
        """
        Initialize the TimestampMatcher with original records.

        Parameters
        ----------
        orig_records : pd.Series
            Original timestamp records.

        Raises
        ------
        TypeError
            If orig_records is not a pandas Series.
        """
        if not isinstance(orig_records, pd.Series):
            raise TypeError("orig_records must be a pandas Series")

        self.orig_records = orig_records
        self.tz = str(self.orig_records.index.tz)

        self.conv_df = pd.DataFrame()  # Holds the conversion DataFrame

        logger.debug("TimestampMatcher initialized with %d records", len(orig_records))

    @property
    def obsname(self) -> str:
        """Return the name of the original records as a string."""
        return str(self.orig_records.name)

    @property
    def target_freq(self) -> pd.Timedelta:
        """Return the inferred target frequency from the conversion DataFrame."""
        freq = pd.infer_freq(self.conv_df["datetimedummy_target"])
        if freq is None:
            raise ValueError("Frequency could not be computed.")
        return to_timedelta(freq)

    @property
    def target_records(self) -> pd.Series:
        """Return the target records after timestamp conversion."""
        if self.conv_df.empty:
            logger.error("No timestamp conversion has been applied yet.")
            raise ValueError("No timestamp conversion has been applied yet.")
        return self.conv_df.set_index("datetime")[self.obsname]

    @property
    def gap_records(self) -> pd.Series:
        """Return the records that have gaps after timestamp conversion."""
        if self.conv_df.empty:
            logger.error("No timestamp conversion has been applied yet.")
            raise ValueError("No timestamp conversion has been applied yet.")
        gapsubset = self.conv_df[self.conv_df["datetimedummy_raw"].isnull()]
        return gapsubset.set_index("datetime")[self.obsname]

    @property
    def outlier_records(self) -> pd.Series:
        """Return the records that are outliers after timestamp conversion."""
        if self.conv_df.empty:
            logger.error("No timestamp conversion has been applied yet.")
            raise ValueError("No timestamp conversion has been applied yet.")

        outliersubset = self.conv_df[
            pd.isna(self.conv_df[self.obsname])
            & pd.notna(self.conv_df["datetimedummy_raw"])
        ]
        return outliersubset.set_index("datetime")[self.obsname]

    @log_entry
    def get_outlier_map(self) -> dict:
        """
        Get a mapping of outlier records.

        Returns
        -------
        dict
            Mapping of outlier records.

        Raises
        ------
        ValueError
            If no timestamp conversion has been applied yet.
        """
        if self.conv_df.empty:
            logger.error("No timestamp conversion has been applied yet.")
            raise ValueError("No timestamp conversion has been applied yet.")

        outliersubset = self.conv_df[
            pd.isna(self.conv_df[self.obsname])
            & pd.notna(self.conv_df["datetimedummy_raw"])
        ]
        return dict(
            zip(
                outliersubset["datetimedummy_raw"],
                outliersubset["datetimedummy_target"],
            )
        )

    def _map_to_perfect_timestamps(
        self,
        target_freq: Union[pd.Timedelta, str],
        shift_tolerance: Union[pd.Timedelta, str],
        origin: Union[pd.Timestamp, None] = None,
        closing: Union[pd.Timestamp, None] = None,
        direction: Literal["backward", "forward", "nearest"] = "nearest",
    ) -> None:
        """
        Map original records to perfect timestamps.

        Parameters
        ----------
        target_freq : pd.Timedelta or str
            Target frequency for mapping.
        shift_tolerance : pd.Timedelta or str
            Tolerance for shifting timestamps.
        origin : pd.Timestamp or None, optional
            Origin timestamp for the mapping. If None, the minimum index is used.
        closing : pd.Timestamp or None, optional
            Closing timestamp for the mapping. If None, it is calculated.
        direction : str, optional
            Direction for merging, by default "nearest".

        Raises
        ------
        AssertionError
            If the closing timestamp is not a valid candidate.
        """
        target_freq = pd.to_timedelta(target_freq)
        # get origin
        if origin is None:
            origin = self.orig_records.index.min()
        else:
            origin = origin

        # get closing timestamp
        if closing is None:
            closing = pd.Timestamp(
                origin
                + (
                    int((self.orig_records.index.max() - origin) / target_freq)
                    * target_freq
                )
            )
        else:
            # test if the closing is a good candidate
            assert (closing - origin) % target_freq == pd.Timedelta(
                0
            ), f"{closing} is not a good candidate with origin: {origin} and freq {target_freq}."

        target_dtrange = pd.date_range(
            start=origin,
            end=closing,
            freq=target_freq,
            tz=self.tz,
        )

        rawdf = self.orig_records.to_frame()
        rawdf["datetimedummy"] = (
            rawdf.index
        )  # To keep track of the original raw timestamps

        targetdf = pd.DataFrame(
            data={"datetimedummy": target_dtrange}, index=target_dtrange
        )
        targetdf.index.name = "datetime"

        logger.debug(
            "Merging target and raw DataFrames with tolerance: %s", shift_tolerance
        )

        mergedf = pd.merge_asof(
            left=targetdf,
            right=rawdf.sort_index(),
            suffixes=("_target", "_raw"),
            left_on="datetime",
            right_on="datetime",
            tolerance=pd.to_timedelta(shift_tolerance),
            direction=direction,
        )

        self.conv_df = mergedf

    @log_entry
    def make_equispaced_timestamps_mapper(
        self,
        freq_estimation_method: Literal["highest", "median"],
        freq_estimation_simplify_tolerance: Union[pd.Timedelta, str],
        origin_simplify_tolerance: Union[pd.Timedelta, str],
        timestamp_tolerance: Union[pd.Timedelta, str],
        force_freq: Union[pd.Timedelta, str, None] = None,
        force_origin: Union[pd.Timestamp, str, None] = None,
        force_closing: Union[pd.Timestamp, None] = None,
    ) -> None:
        """
        Create a mapper for irregular records to perfect-frequency timestamps.

        Parameters
        ----------
        freq_estimation_method : {'highest', 'median'}
            Method to estimate frequency.
        freq_estimation_simplify_tolerance : pd.Timedelta or str
            Tolerance for simplifying frequency estimation.
        origin_simplify_tolerance : pd.Timedelta or str
            Tolerance for simplifying origin.
        timestamp_tolerance : pd.Timedelta or str
            Tolerance for timestamp matching.
        force_freq : pd.Timedelta, str, or None, optional
            Force the frequency to a specific value.
        force_origin : pd.Timestamp, str, or None, optional
            Force the origin to a specific value.
        force_closing : pd.Timestamp or None, optional
            Force the closing timestamp.

        Raises
        ------
        TypeError
            If input types are incorrect.
        """
        if not isinstance(
            freq_estimation_method, str
        ) or freq_estimation_method not in ["highest", "median"]:
            raise TypeError("freq_estimation_method must be 'highest' or 'median'")
        if not isinstance(freq_estimation_simplify_tolerance, (pd.Timedelta, str)):
            raise TypeError(
                "freq_estimation_simplify_tolerance must be a pd.Timedelta or str"
            )
        if not isinstance(origin_simplify_tolerance, (pd.Timedelta, str)):
            raise TypeError("origin_simplify_tolerance must be a pd.Timedelta or str")
        if not isinstance(timestamp_tolerance, (pd.Timedelta, str)):
            raise TypeError("timestamp_tolerance must be a pd.Timedelta or str")

        if force_freq is None:
            logger.debug(
                "Estimating target frequency using method: %s", freq_estimation_method
            )

            target_freq = get_likely_frequency(
                timestamps=self.orig_records.index,
                method=freq_estimation_method,
                max_simplify_error=freq_estimation_simplify_tolerance,
            )
        else:
            target_freq = pd.to_timedelta(force_freq)
            logger.debug(f"Force the frequency to {target_freq}")

        if force_origin is None:
            logger.debug(
                "Simplifying origin timestamp with tolerance: %s",
                origin_simplify_tolerance,
            )

            target_origin = simplify_time(
                time=self.orig_records.index.min(),
                max_simplify_error=origin_simplify_tolerance,
            )

        else:
            target_origin = pd.to_datetime(force_origin)
            logger.debug(f"Force origin to {target_origin}")

        self._map_to_perfect_timestamps(
            target_freq=target_freq,
            shift_tolerance=timestamp_tolerance,
            origin=target_origin,
            closing=force_closing,
        )


@log_entry
def simplify_time(
    time: Union[pd.Timestamp, pd.Timedelta],
    max_simplify_error: pd.Timedelta,
    zero_protection: bool = False,
) -> Union[pd.Timestamp, pd.Timedelta]:
    """
    Simplify a time (or timedelta) to a rounded value, within a tolerance.

    Parameters
    ----------
    time : pd.Timestamp or pd.Timedelta
        The time or timedelta to simplify.
    max_simplify_error : pd.Timedelta
        The maximum allowed deviation from the simplified value.
    zero_protection : bool, optional
        If True, ensures the returned simplified time is not zero. Used when applied to frequencies.

    Returns
    -------
    pd.Timestamp or pd.Timedelta
        The simplified time or timedelta.

    Raises
    ------
    MetObsTimeSimplifyError
        If no simplification is possible and zero_protection is True.
    """
    # NOTE: Make sure that the sequence goes from coarse to fine AND
    # That all elements are natural multiplicatives of each other !! (this
    # is required since the synchronization relies on it)
    simplify_sequence = ["1d", "1h", "30min", "10min", "5min", "1min", "30s"]

    for simpl_resolution in simplify_sequence:
        candidate = time.round(pd.Timedelta(simpl_resolution))
        if zero_protection:
            if candidate == pd.Timedelta(0):
                # try the ceil as candidate (is never zero)
                candidate = time.ceil(pd.Timedelta(simpl_resolution))

        # Tests if candidate meets conditions
        if abs(time - candidate) < pd.to_timedelta(max_simplify_error):
            return candidate

    # No simplification possible
    if zero_protection:
        if time == pd.Timedelta(0):
            raise MetObsTimeSimplifyError(
                f"No simplification possible for {time}, and zero_protection is set to True."
            )
    return time


@log_entry
def get_likely_frequency(
    timestamps: pd.DatetimeIndex,
    method: Literal["highest", "median"] = "highest",
    max_simplify_error: Union[str, pd.Timedelta] = "2min",
) -> pd.Timedelta:
    """
    Find the most likely observation frequency of a DatetimeIndex.

    Parameters
    ----------
    timestamps : pd.DatetimeIndex
        DatetimeIndex of the dataset.
    method : {'highest', 'median'}, optional
        Select which method to use. If 'highest', the highest appearing frequency is used.
        If 'median', the median of the appearing frequencies is used. Default is 'highest'.
    max_simplify_error : str or pd.Timedelta, optional
        The maximum deviation from the found frequency when simplifying. Default is '2min'.

    Returns
    -------
    pd.Timedelta
        The assumed (and simplified) frequency of the DatetimeIndex.

    Raises
    ------
    MetObsTimeSimplifyError
        If max_simplify_error is not a valid time indication.
    AssertionError
        If the method is not known.
    """
    logger.debug(
        f"Starting get_likely_frequency with method={method} and max_simplify_error={max_simplify_error}"
    )
    assert method in [
        "highest",
        "median",
    ], f"The method for frequency estimation ({method}) is not known. Use one of [highest, median]"

    try:
        pd.to_timedelta(max_simplify_error)
    except ValueError:
        raise MetObsTimeSimplifyError(
            f'{max_simplify_error} is not valid time indication. Example: "5min" indicates 5 minutes.'
        )

    # simplify is True if a non-zero simplify_error is provided
    if pd.to_timedelta(max_simplify_error).seconds < 1:
        simplify = False
    else:
        simplify = True

    logger.debug(f"Simplify is set to {simplify}")

    freqs_blacklist = [pd.Timedelta(0), np.nan]  # avoid a zero frequency

    freqs = timestamps.to_series().diff()
    freqs = freqs[~freqs.isin(freqs_blacklist)]

    if method == "highest":
        assume_freq = freqs.min()  # highest frequency
    elif method == "median":
        assume_freq = freqs.median()

    logger.debug(f"Assumed frequency before simplification: {assume_freq}")

    # Check if frequency estimation failed (empty timestamps or no valid differences)
    if pd.isna(assume_freq):
        raise MetObsTimeSimplifyError(
            f"Cannot estimate frequency from the provided timestamps. "
            f"This typically occurs when there are no valid timestamps, "
            f"only a single timestamp, or all timestamps are identical. "
        )

    if simplify:
        assume_freq = simplify_time(
            time=assume_freq,
            max_simplify_error=max_simplify_error,
            zero_protection=True,
        )

        logger.debug(f"Assumed frequency after simplification: {assume_freq}")

    logger.debug(f"Final assumed frequency: {assume_freq}")

    return pd.to_timedelta(assume_freq)
