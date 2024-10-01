#!/usr/bin/env python3

"""
Module for handling and analyzing data differences between fresh MDSplus data
and cached SQL data.
"""

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from disruption_py.config import config
from disruption_py.core.utils.misc import safe_cast


@dataclass
class DataDifference:
    """
    Data difference between fresh MDSplus data and cached SQL data.
    """

    shot_id: int
    data_column: str

    missing_cache_data: bool
    missing_fresh_data: bool

    anomalies: np.ndarray = field(init=False)  # 1 if anomaly, 0 o.w.
    relative_difference: np.ndarray = field(init=False)
    fresh_column_data: pd.Series
    cache_column_data: pd.Series
    expect_failure: bool

    fresh_time: pd.Series
    cache_time: pd.Series

    def __post_init__(self):
        self.anomalies, self.relative_difference = self.compute_numeric_anomalies()

    @property
    def num_anomalies(self) -> int:
        """Sum the number of anomalies"""
        return np.sum(self.anomalies)

    @property
    def timebase_length(self) -> int:
        """Get the timebase length"""
        return len(self.anomalies)

    @property
    def missing_data(self) -> bool:
        """Return True if either fresh or cache is missing data."""
        return self.missing_cache_data or self.missing_fresh_data

    @property
    def failed(self) -> str:
        """Return True if missing data or if there are too many anomalies."""
        if self.missing_data:
            return True
        return (
            self.num_anomalies / self.timebase_length
            > 1 - config().testing.MATCH_FRACTION
        )

    @property
    def failure_ratio_string(self) -> str:
        """Create a string for the fraction of datapoints that are anomalies"""
        return f"{self.num_anomalies / self.timebase_length:.4f}"

    @property
    def column_mismatch_string(self) -> str:
        """Create a string showing the difference between fresh and cache data."""
        # Missing data handled here because difference_df expects data to exist
        s = f"Shot {self.shot_id} column {self.data_column}"
        if self.missing_cache_data or self.missing_fresh_data:
            fresh_str = (
                "Missing fresh data" if self.missing_fresh_data else "Have fresh data"
            )
            cache_str = (
                "Missing cache data" if self.missing_cache_data else "Have cache data"
            )
            return f"{s}: {fresh_str} and {cache_str}"
        return s + f" with arrays:\n{self.difference_df.to_string()}"

    @property
    def difference_df(self) -> pd.DataFrame:
        """
        Create a dataframe with columns for time, fresh data, cache data, the
        ratio between the two data, the absolute difference between them, the relative
        difference, and whether the point is an anomaly.
        """
        indices = (
            np.arange(self.timebase_length)
            if config().testing.VERBOSE_OUTPUT
            else self.anomalies.flatten()
        )
        anomaly = self.anomalies[indices]
        fresh_data = self.fresh_column_data.iloc[indices]
        cache_data = self.cache_column_data.iloc[indices]
        return pd.DataFrame(
            {
                "Time": self.fresh_time[indices],
                "Fresh Data": fresh_data,
                "Cache Data": cache_data,
                "Fresh/Cache": fresh_data / cache_data,
                "Absolute difference": abs(fresh_data - cache_data),
                "Relative difference": self.relative_difference[indices],
                "Anomaly": anomaly,
            }
        )

    def compute_numeric_anomalies(self):
        """
        Get the indices of the data where numeric differences exist between the
        cached and fresh data.
        """

        # handle missing data case
        if self.missing_fresh_data or self.missing_cache_data:
            if self.missing_fresh_data and self.missing_cache_data:
                missing_timebase_length = 0
            elif self.missing_fresh_data:
                missing_timebase_length = len(self.cache_column_data)
            else:
                missing_timebase_length = len(self.fresh_column_data)
            return np.ones(missing_timebase_length, dtype=bool), np.zeros(
                missing_timebase_length
            )

        cache_is_nan = pd.isnull(self.cache_column_data)
        fresh_is_nan = pd.isnull(self.fresh_column_data)

        # handle case where both arrays are all null
        if cache_is_nan.all() and fresh_is_nan.all():
            return np.zeros(len(self.fresh_column_data), dtype=bool), np.zeros(
                len(self.fresh_column_data)
            )

        relative_difference = safe_cast(
            np.where(
                self.cache_column_data != 0,
                np.abs(
                    (self.fresh_column_data - self.cache_column_data)
                    / self.cache_column_data
                ),
                np.where(self.fresh_column_data != 0, np.inf, np.nan),
            ),
            "float64",
        )  # necessary in case all produced values are nan

        numeric_anomalies_mask = np.where(
            np.isnan(relative_difference),
            False,
            relative_difference > config().testing.VAL_TOLERANCE,
        )
        nan_anomalies_mask = cache_is_nan != fresh_is_nan
        anomalies: pd.Series = numeric_anomalies_mask | nan_anomalies_mask

        return anomalies.to_numpy(), relative_difference


def assert_frame_equal_unordered(df1: pd.DataFrame, df2: pd.DataFrame):
    """Compare whether two dataframes have the same values."""
    df1_sorted = df1.sort_values(by=config().database.protected_columns).reset_index(
        drop=True
    )
    df2_sorted = df2.sort_values(by=config().database.protected_columns).reset_index(
        drop=True
    )
    pd.testing.assert_frame_equal(df1_sorted, df2_sorted, check_like=True)
