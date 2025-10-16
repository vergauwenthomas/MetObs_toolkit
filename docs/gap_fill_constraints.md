# Gap-Filling with Value Constraints

This document describes how to use the `min_value` and `max_value` parameters when filling gaps in observational data.

## Problem

When filling gaps using interpolation or model-based methods (like debiased model data), the filled values can sometimes be physically unrealistic. For example:

- Relative humidity (RH) values exceeding 100%
- Wind speed values becoming negative
- Temperature values outside physically possible ranges

## Solution

All gap-filling methods now support optional `min_value` and `max_value` parameters that constrain filled values to physically realistic ranges.

## Usage Examples

### Example 1: Constraining Humidity to 0-100%

```python
import metobs_toolkit

# Load your dataset
dataset = metobs_toolkit.Dataset()
dataset.import_data_from_file(...)

# Fill gaps for humidity with constraints
dataset.fill_gaps_with_diurnal_debiased_modeldata(
    target_obstype='humidity',
    modeldataset=era5_data,
    leading_period_duration='7days',
    trailing_period_duration='7days',
    min_debias_sample_size=30,
    min_value=0.0,    # Humidity cannot be negative
    max_value=100.0   # Humidity cannot exceed 100%
)
```

### Example 2: Constraining Wind Speed to Non-negative Values

```python
# Fill gaps for wind speed - wind speed cannot be negative
dataset.fill_gaps_with_debiased_modeldata(
    target_obstype='wind_speed',
    modeldataset=era5_data,
    leading_period_duration='3days',
    trailing_period_duration='3days',
    min_leading_records_total=50,
    min_trailing_records_total=50,
    min_value=0.0   # Wind speed cannot be negative
)
```

### Example 3: Using Constraints with Interpolation

```python
# Use polynomial interpolation with constraints
station.interpolate_gaps(
    target_obstype='humidity',
    method='polynomial',
    max_gap_duration_to_fill='5h',
    n_leading_anchors=3,
    n_trailing_anchors=3,
    min_value=0.0,
    max_value=100.0,
    method_kwargs={'order': 3}
)
```

### Example 4: At the SensorData Level

```python
# Get humidity sensor data
humidity_sensor = station.get_sensor('humidity')

# Fill gaps with constraints
humidity_sensor.interpolate_gaps(
    method='linear',
    max_gap_duration_to_fill='3h',
    min_value=0.0,
    max_value=100.0
)
```

## Supported Methods

The `min_value` and `max_value` parameters are supported by:

1. **Model-based gap filling:**
   - `raw_model_gapfill()` - Raw model data without bias correction
   - `debiased_model_gapfill()` - Model data with bias correction
   - `diurnal_debiased_model_gapfill()` - Model data with diurnal bias correction
   - `weighted_diurnal_debiased_model_gapfill()` - Model data with weighted diurnal bias correction

2. **Interpolation-based gap filling:**
   - `interpolate()` - All interpolation methods (linear, polynomial, spline, etc.)

## Implementation Details

- Values are clipped **after** the gap-filling computation but **before** updating the gap records
- If only `max_value` is specified, only the upper bound is enforced
- If only `min_value` is specified, only the lower bound is enforced
- Both parameters are optional and default to `None` (no constraints)
- The feature is backward compatible - existing code without these parameters will work unchanged

## Performance

The clipping operation uses pandas `.clip()` method which is highly optimized and has minimal performance impact.

## Related

- Issue: [gap fill generates unphysical data](https://github.com/vergauwenthomas/MetObs_toolkit/issues/XXX)
- Test file: `tests/test_gf_limits.py`
