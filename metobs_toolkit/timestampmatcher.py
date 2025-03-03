import logging
from typing import Literal, Tuple
import pandas as pd
from metobs_toolkit.backend_collection.df_helpers import (
    get_likely_frequency,
    simplify_time,
)

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
        # note: sometimes 'h' is returned, and this gives issues, so add a 1 in front
        if not freq[0].isdigit():
            freq = "1" + freq

        return pd.Timedelta(freq)

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
            origin = self.orig_records.index.min().floor(
                target_freq
            )  # floor target freq!
        else:
            origin = origin.floor(target_freq)  # make sure origin is floored

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
        freq_estimation_simplify_tolerance: pd.Timedelta | str,
        origin_simplify_tolerance: pd.Timedelta | str,
        timestamp_tolerance: pd.Timedelta | str,
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
                time=self.orig_records.index.min().floor(target_freq),
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
