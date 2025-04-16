import logging
from typing import Literal, Union
import pandas as pd
import numpy as np


from metobs_toolkit.backend_collection.errorclasses import *
from metobs_toolkit.backend_collection.df_helpers import to_timedelta


logger = logging.getLogger(__file__)


class TimestampMatcher:
    """Class for mapping and resampling timestamp records."""

    def __init__(self, orig_records: pd.Series):
        """
        Initialize the TimestampMatcher with original records.

        Parameters
        ----------
        orig_records : pd.Series
            Original timestamp records.
        """
        if not isinstance(orig_records, pd.Series):
            raise TypeError("orig_records must be a pandas Series")

        self.orig_records = orig_records
        self.tz = str(self.orig_records.index.tz)  # TYPO

        self.conv_df = pd.DataFrame()  # Holds the conversion dataframe

        logger.info("TimestampMatcher initialized with %d records", len(orig_records))

    @property
    def obsname(self) -> str:
        """Get the name of the original records.

        Returns
        -------
        str
            Name of the original records.
        """
        return str(self.orig_records.name)  # TYPO

    @property
    def target_freq(self) -> pd.Timedelta:
        """Infer the target frequency from the conversion dataframe.

        Returns
        -------
        pd.Timedelta
            Inferred target frequency.
        """
        freq = pd.infer_freq(self.conv_df["datetimedummy_target"])
        if freq is None:
            raise ValueError("Frequency could not be computed.")
        return to_timedelta(freq)

    @property
    def target_records(self) -> pd.Series:
        """
        Get the target records after timestamp conversion.

        Returns
        -------
        pd.Series
            Target records with converted timestamps.
        """
        if self.conv_df.empty:
            logger.error("No timestamp conversion has been applied yet.")
            raise ValueError("No timestamp conversion has been applied yet.")
        return self.conv_df.set_index("datetime")[self.obsname]  # TYPO

    @property
    def gap_records(self) -> pd.Series:
        """
        Get the records that have gaps after timestamp conversion.

        Returns
        -------
        pd.Series
            Records with gaps.
        """
        if self.conv_df.empty:
            logger.error("No timestamp conversion has been applied yet.")
            raise ValueError("No timestamp conversion has been applied yet.")
        gapsubset = self.conv_df[self.conv_df["datetimedummy_raw"].isnull()]  # TYPO
        return gapsubset.set_index("datetime")[self.obsname]  # TYPO

    @property
    def outlier_records(self) -> pd.Series:
        """
        Get the records that are outliers after timestamp conversion.

        Returns
        -------
        pd.Series
            Outlier records.
        """
        if self.conv_df.empty:
            logger.error("No timestamp conversion has been applied yet.")
            raise ValueError("No timestamp conversion has been applied yet.")

        outliersubset = self.conv_df[
            pd.isna(self.conv_df[self.obsname])
            & pd.notna(self.conv_df["datetimedummy_raw"])  # TYPO
        ]  # TYPO
        return outliersubset.set_index("datetime")[self.obsname]  # TYPO

    def get_outlier_map(self) -> dict:
        """
        Get a mapping of outlier records.

        Returns
        -------
        dict
            Mapping of outlier records.
        """
        if self.conv_df.empty:
            logger.error("No timestamp conversion has been applied yet.")
            raise ValueError("No timestamp conversion has been applied yet.")

        outliersubset = self.conv_df[
            pd.isna(self.conv_df[self.obsname])
            & pd.notna(self.conv_df["datetimedummy_raw"])  # TYPO
        ]  # TYPO
        return dict(
            zip(
                outliersubset["datetimedummy_raw"],
                outliersubset["datetimedummy_target"],
            )
        )  # TYPO

    def _map_to_perfect_timestamps(
        self,
        target_freq,
        shift_tolerance,
        origin=None,
        closing=None,
        direction="nearest",
    ):

        target_freq = pd.to_timedelta(target_freq)
        # # check if the records are perfect
        # freq = pd.infer_freq(self.orig_records.index)
        # if freq is None:
        #     raise ValueError(
        #         f"{self} does not hold perfec-freq-records, and thus reindexing is not possible without introducing errors"
        #     )

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
            tz=self.tz,  # TYPO
        )

        rawdf = self.orig_records.to_frame()
        rawdf["datetimedummy"] = (
            rawdf.index
        )  # To keep track of the original raw timestamps

        targetdf = pd.DataFrame(
            data={"datetimedummy": target_dtrange}, index=target_dtrange  # TYPO
        )
        targetdf.index.name = "datetime"

        logger.debug(
            "Merging target and raw dataframes with tolerance: %s", shift_tolerance
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

    def make_equispaced_timestamps_mapper(
        self,
        freq_estimation_method: Literal["highest", "median"],
        freq_estimation_simplify_tolerance: Union[pd.Timedelta, str],
        origin_simplify_tolerance: Union[pd.Timedelta, str],
        timestamp_tolerance: Union[pd.Timedelta, str],
        force_freq=None,
        force_origin=None,
        force_closing=None,
    ) -> None:
        """
        Create a mapper for irregular records to perfect-freq-timestamps.

        Parameters
        ----------
        freq_estimation_method : Literal['highest', 'median']
            Method to estimate frequency.
        freq_estimation_simplify_tolerance : pd.Timedelta or str
            Tolerance for simplifying frequency estimation.
        origin_simplify_tolerance : pd.Timedelta or str
            Tolerance for simplifying origin.
        timestamp_tolerance : pd.Timedelta or str
            Tolerance for timestamp matching.

        Returns
        -------
        Tuple[pd.Series, pd.DataFrame]
            Target series and conversion dataframe.
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
                timestamps=self.orig_records.index,  # TYPO
                method=freq_estimation_method,
                max_simplify_error=freq_estimation_simplify_tolerance,
            )  # TYPO
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


def simplify_time(
    time: Union[pd.Timestamp, pd.Timedelta],
    max_simplify_error: pd.Timedelta,
    zero_protection: bool = False,
):
    """Simplifies a time (or timedelta) to a rounded value, within a tolerance.

    time can be a timestamp of a timedelta.

    zero_protections make shure that the returned simplified time is not zero.
    This is in practice used when applied to frequencies.
    """

    # NOTE: Make sure that the sequence goes from coarce to fine AND
    # That all elements are natural multiplicatives of each other !! (this
    # is required since the syncronization relies on it)
    simplify_sequence = ["1d", "1h", "30min", "10min", "5min", "1min", "30s"]

    for simpl_resolution in simplify_sequence:
        candidate = time.round(pd.Timedelta(simpl_resolution))
        if zero_protection:
            if candidate == pd.Timedelta(0):
                # try the ceil as candidate (is never zero)
                candidate = time.ceil(pd.Timedelta(simpl_resolution))

        # Tests if candidate mets conditions
        if abs(time - candidate) < pd.to_timedelta(max_simplify_error):
            return candidate

    # No simplyfication posible
    if zero_protection:
        if time == pd.Timedelta(0):
            raise MetObsTimeSimplifyError(
                f"No simplification possible for {time}, and zero_protection is set to True."
            )
    return time


def get_likely_frequency(
    timestamps: pd.DatetimeIndex,
    method: Literal["highest", "median"] = "highest",
    max_simplify_error: Union[str, pd.Timedelta] = "2min",
) -> pd.Timedelta:
    """Find the most likely observation frequency of a datetimeindex.

    Parameters
    ----------
    timestamps : pandas.Datetimeindex()
        Datetimeindex of the dataset.df.
    method : 'highest' or 'median', optional
        Select wich method to use. If 'highest', the highest apearing frequency is used.
        If 'median', the median of the apearing frequencies is used. The default is 'highest'.
    max_simplify_error : datetimestring, optional
        The maximum deviation from the found frequency when simplifying. The default is '2min'.

    Returns
    -------
    assume_freq : datetime.timedelta
        The assumed (and simplified) frequency of the datetimeindex.

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
            f'{max_simplify_error} is not valid timeindication. Example: "5min" indicates 5 minutes.'
        )

    # simplify is true if a non-zero simplify_error is provided
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

    if simplify:
        assume_freq = simplify_time(
            time=assume_freq,
            max_simplify_error=max_simplify_error,
            zero_protection=True,
        )

        logger.debug(f"Assumed frequency after simplification: {assume_freq}")

    # if assume_freq == pd.to_timedelta(0):  # highly likely due to a duplicated record
    #     # select the second highest frequency
    #     assume_freq = abs(
    #         timestamps.to_series().diff().value_counts().index
    #     ).sort_values(ascending=True)[1]

    logger.debug(f"Final assumed frequency: {assume_freq}")

    return pd.to_timedelta(assume_freq)
