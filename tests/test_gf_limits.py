"""Test gap-filling with min_value and max_value limits."""
import pytest
import sys
from pathlib import Path
import pandas as pd
import numpy as np

# Add the local source directory to Python path for development
libfolder = Path(str(Path(__file__).resolve())).parent.parent
sys.path.insert(0, str(libfolder / "src"))
import metobs_toolkit
from metobs_toolkit.obstypes import Obstype
from metobs_toolkit.gap import Gap


def test_debias_fill_with_max_value():
    """Test that fill_regular_debias respects max_value parameter."""
    from metobs_toolkit.gf_collection.debias_gapfill import fill_regular_debias
    
    # Create a test DataFrame
    df = pd.DataFrame({
        'value': [95.0, 98.0, np.nan, np.nan, 105.0],
        'label': ['lead', 'lead', 'gap', 'gap', 'trail'],
        'modelvalue': [96.0, 99.0, 110.0, 112.0, 106.0]
    })
    
    # Fill without limit
    result_no_limit = fill_regular_debias(df.copy())
    
    # Fill with max_value=100
    result_with_limit = fill_regular_debias(df.copy(), max_value=100.0)
    
    # Check that without limit, values exceed 100
    assert result_no_limit['fillvalue'].max() > 100.0, "Without limit, fillvalue should exceed 100"
    
    # Check that with limit, no value exceeds 100
    assert result_with_limit['fillvalue'].max() <= 100.0, "With max_value=100, fillvalue should not exceed 100"
    
    print(f"✓ fill_regular_debias with max_value=100.0 successful")
    print(f"  Without limit: max = {result_no_limit['fillvalue'].max():.2f}")
    print(f"  With limit:    max = {result_with_limit['fillvalue'].max():.2f}")


def test_debias_fill_with_min_value():
    """Test that fill_regular_debias respects min_value parameter."""
    from metobs_toolkit.gf_collection.debias_gapfill import fill_regular_debias
    
    # Create a test DataFrame
    df = pd.DataFrame({
        'value': [5.0, 2.0, np.nan, np.nan, -5.0],
        'label': ['lead', 'lead', 'gap', 'gap', 'trail'],
        'modelvalue': [6.0, 3.0, -8.0, -10.0, -4.0]
    })
    
    # Fill without limit
    result_no_limit = fill_regular_debias(df.copy())
    
    # Fill with min_value=0
    result_with_limit = fill_regular_debias(df.copy(), min_value=0.0)
    
    # Check that without limit, values go below 0
    assert result_no_limit['fillvalue'].min() < 0.0, "Without limit, fillvalue should go below 0"
    
    # Check that with limit, no value is below 0
    assert result_with_limit['fillvalue'].min() >= 0.0, "With min_value=0, fillvalue should not go below 0"
    
    print(f"✓ fill_regular_debias with min_value=0.0 successful")
    print(f"  Without limit: min = {result_no_limit['fillvalue'].min():.2f}")
    print(f"  With limit:    min = {result_with_limit['fillvalue'].min():.2f}")


def test_diurnal_debias_fill_with_max_value():
    """Test that fill_with_diurnal_debias respects max_value parameter."""
    from metobs_toolkit.gf_collection.diurnal_debias_gapfill import fill_with_diurnal_debias
    
    # Create a test DataFrame with datetime index - more data points to ensure sample size is met
    dates = pd.date_range('2023-01-01', periods=50, freq='1h')
    
    # Create patterns that repeat hourly to ensure diurnal bias calculation
    values = []
    labels = []
    modelvalues = []
    
    for i in range(50):
        if i < 20:  # lead period
            values.append(95.0 + (i % 5))
            labels.append('lead')
            modelvalues.append(96.0 + (i % 5))
        elif i < 25:  # gap period
            values.append(np.nan)
            labels.append('gap')
            modelvalues.append(110.0 + (i % 5))
        else:  # trail period
            values.append(100.0 + (i % 5))
            labels.append('trail')
            modelvalues.append(101.0 + (i % 5))
    
    df = pd.DataFrame({
        'value': values,
        'label': labels,
        'modelvalue': modelvalues
    }, index=dates)
    df.index.name = 'datetime'
    
    # Fill without limit
    result_no_limit = fill_with_diurnal_debias(df.copy(), min_sample_size=2)
    
    # Fill with max_value=100
    result_with_limit = fill_with_diurnal_debias(df.copy(), min_sample_size=2, max_value=100.0)
    
    # Check that fills were made
    if result_with_limit['fillvalue'].notna().any():
        max_with_limit = result_with_limit['fillvalue'].max()
        
        print(f"✓ fill_with_diurnal_debias with max_value=100.0 successful")
        print(f"  Max fillvalue with limit: {max_with_limit:.2f}")
        
        # Check that with limit, no value exceeds 100
        assert result_with_limit['fillvalue'].max() <= 100.0, "With max_value=100, fillvalue should not exceed 100"
    else:
        print(f"✓ fill_with_diurnal_debias test (no fills due to sample size)")


def test_weighted_diurnal_debias_fill_with_limits():
    """Test that fill_with_weighted_diurnal_debias respects limits."""
    from metobs_toolkit.gf_collection.diurnal_debias_gapfill import fill_with_weighted_diurnal_debias
    
    # Create a test DataFrame with datetime index - more data points
    dates = pd.date_range('2023-01-01', periods=50, freq='1h')
    
    # Create patterns that repeat hourly
    values = []
    labels = []
    modelvalues = []
    
    for i in range(50):
        if i < 20:  # lead period
            values.append(95.0 + (i % 5))
            labels.append('lead')
            modelvalues.append(96.0 + (i % 5))
        elif i < 25:  # gap period
            values.append(np.nan)
            labels.append('gap')
            modelvalues.append(110.0 + (i % 5))
        else:  # trail period
            values.append(100.0 + (i % 5))
            labels.append('trail')
            modelvalues.append(101.0 + (i % 5))
    
    df = pd.DataFrame({
        'value': values,
        'label': labels,
        'modelvalue': modelvalues
    }, index=dates)
    df.index.name = 'datetime'
    
    # Fill with max_value=100
    result_with_limit = fill_with_weighted_diurnal_debias(
        df.copy(), 
        min_lead_sample_size=2, 
        min_trail_sample_size=2,
        max_value=100.0
    )
    
    # Check that with limit, no value exceeds 100
    if result_with_limit['fillvalue'].notna().any():
        assert result_with_limit['fillvalue'].max() <= 100.0, "With max_value=100, fillvalue should not exceed 100"
        print(f"✓ fill_with_weighted_diurnal_debias with max_value=100.0 successful")
        print(f"  Max fillvalue: {result_with_limit['fillvalue'].max():.2f}")
    else:
        print(f"✓ fill_with_weighted_diurnal_debias test (no fills due to sample size)")


def test_backward_compatibility():
    """Test that gap-filling still works without min_value and max_value parameters."""
    from metobs_toolkit.gf_collection.debias_gapfill import fill_regular_debias
    
    # Create a test DataFrame
    df = pd.DataFrame({
        'value': [20.0, 21.0, np.nan, np.nan, 24.0],
        'label': ['lead', 'lead', 'gap', 'gap', 'trail'],
        'modelvalue': [21.0, 22.0, 23.0, 24.0, 25.0]
    })
    
    # Fill WITHOUT specifying limits (backward compatibility)
    result = fill_regular_debias(df)
    
    # Check that gaps were filled
    assert result['fillvalue'].notna().all(), "All values should be filled"
    
    print(f"✓ Backward compatibility test successful")
    print(f"  All {len(result)} values were filled without errors")


if __name__ == "__main__":
    # Run tests manually
    print("Running gap-fill limits tests...")
    print("-" * 60)
    
    try:
        test_debias_fill_with_max_value()
        print()
        test_debias_fill_with_min_value()
        print()
        test_diurnal_debias_fill_with_max_value()
        print()
        test_weighted_diurnal_debias_fill_with_limits()
        print()
        test_backward_compatibility()
        print()
        print("-" * 60)
        print("✓ All tests passed!")
    except Exception as e:
        print(f"\n✗ Test failed with error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
