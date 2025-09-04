"""This file tests the TimestampMatcher class and related functions using the pytest framework."""

import pytest
import sys
from pathlib import Path
import pandas as pd

# Add the local source directory to Python path for development
libfolder = Path(str(Path(__file__).resolve())).parent.parent
sys.path.insert(0, str(libfolder / "src"))

from metobs_toolkit.timestampmatcher import get_likely_frequency
from metobs_toolkit.backend_collection.errorclasses import MetObsTimeSimplifyError


class TestTimestampMatcher:
    """Test class for TimestampMatcher functionality."""

    def test_get_likely_frequency_normal_case(self):
        """Test that get_likely_frequency works correctly with normal timestamps."""
        # Create valid timestamps with a regular frequency
        timestamps = pd.date_range('2023-01-01', periods=10, freq='1h')
        
        # Test median method
        freq = get_likely_frequency(timestamps, method="median")
        assert freq == pd.Timedelta('1h')
        
        # Test highest method
        freq = get_likely_frequency(timestamps, method="highest")
        assert freq == pd.Timedelta('1h')

    def test_get_likely_frequency_empty_timestamps(self):
        """Test that get_likely_frequency raises appropriate error for empty timestamps."""
        empty_timestamps = pd.DatetimeIndex([])
        
        with pytest.raises(MetObsTimeSimplifyError) as exc_info:
            get_likely_frequency(empty_timestamps, method="median")
        
        error_msg = str(exc_info.value)
        assert "Cannot estimate frequency from the provided timestamps" in error_msg
        assert "Number of timestamps: 0" in error_msg
        assert "Number of valid differences: 0" in error_msg

    def test_get_likely_frequency_single_timestamp(self):
        """Test that get_likely_frequency raises appropriate error for single timestamp."""
        single_timestamp = pd.DatetimeIndex(['2023-01-01 00:00:00'])
        
        with pytest.raises(MetObsTimeSimplifyError) as exc_info:
            get_likely_frequency(single_timestamp, method="median")
        
        error_msg = str(exc_info.value)
        assert "Cannot estimate frequency from the provided timestamps" in error_msg
        assert "Number of timestamps: 1" in error_msg
        assert "Number of valid differences: 0" in error_msg

    def test_get_likely_frequency_identical_timestamps(self):
        """Test that get_likely_frequency raises appropriate error for identical timestamps."""
        identical_timestamps = pd.DatetimeIndex([
            '2023-01-01 00:00:00', 
            '2023-01-01 00:00:00'
        ])
        
        with pytest.raises(MetObsTimeSimplifyError) as exc_info:
            get_likely_frequency(identical_timestamps, method="median")
        
        error_msg = str(exc_info.value)
        assert "Cannot estimate frequency from the provided timestamps" in error_msg
        assert "Number of timestamps: 2" in error_msg
        assert "Number of valid differences: 0" in error_msg

    def test_get_likely_frequency_irregular_timestamps(self):
        """Test that get_likely_frequency works with irregular timestamps."""
        # Create irregular timestamps
        timestamps = pd.DatetimeIndex([
            '2023-01-01 00:00:00',
            '2023-01-01 01:00:00', 
            '2023-01-01 02:30:00',
            '2023-01-01 04:00:00'
        ])
        
        # Should still work even with irregular timestamps
        freq = get_likely_frequency(timestamps, method="median")
        assert isinstance(freq, pd.Timedelta)
        assert freq > pd.Timedelta(0)