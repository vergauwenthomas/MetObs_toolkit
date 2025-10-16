#!/usr/bin/env python
"""
Example: Gap-Filling with Value Constraints for Relative Humidity

This example demonstrates how to use min_value and max_value parameters
to prevent unphysical values when filling gaps in humidity data.
"""

import numpy as np
import pandas as pd

# Demonstrate the core functionality with the gap-filling functions directly
from metobs_toolkit.gf_collection.debias_gapfill import fill_regular_debias
from metobs_toolkit.gf_collection.diurnal_debias_gapfill import (
    fill_with_diurnal_debias,
    fill_with_weighted_diurnal_debias
)

print("=" * 70)
print("Gap-Filling with Value Constraints Example")
print("=" * 70)
print()

# Example 1: Regular debias gap-filling
print("Example 1: Regular Debiased Gap-Fill")
print("-" * 70)

# Create synthetic humidity data where model overpredicts
df = pd.DataFrame({
    'value': [95.0, 98.0, np.nan, np.nan, 105.0],  # Note: 105% is impossible!
    'label': ['lead', 'lead', 'gap', 'gap', 'trail'],
    'modelvalue': [96.0, 99.0, 110.0, 112.0, 106.0]  # Model overpredicts
})

print("Input data:")
print(df)
print()

# Fill WITHOUT constraints
result_unconstrained = fill_regular_debias(df.copy())
print("Without constraints:")
print(f"  Filled values: {result_unconstrained.loc[df['label'] == 'gap', 'fillvalue'].values}")
print(f"  Max value: {result_unconstrained['fillvalue'].max():.2f}")
print()

# Fill WITH constraints
result_constrained = fill_regular_debias(df.copy(), min_value=0.0, max_value=100.0)
print("With constraints (min=0, max=100):")
print(f"  Filled values: {result_constrained.loc[df['label'] == 'gap', 'fillvalue'].values}")
print(f"  Max value: {result_constrained['fillvalue'].max():.2f}")
print()

# Example 2: Demonstrate with diurnal debias
print("\nExample 2: Diurnal Debiased Gap-Fill")
print("-" * 70)

# Create more realistic time-series data
dates = pd.date_range('2023-07-01', periods=50, freq='1h')

values = []
labels = []
modelvalues = []

for i in range(50):
    if i < 20:  # lead period
        values.append(92.0 + (i % 5) * 1.5)  # Values around 92-99%
        labels.append('lead')
        modelvalues.append(93.0 + (i % 5) * 1.5)
    elif i < 25:  # gap period
        values.append(np.nan)
        labels.append('gap')
        modelvalues.append(105.0 + (i % 5) * 2.0)  # Model significantly overpredicts
    else:  # trail period
        values.append(96.0 + (i % 5) * 1.0)
        labels.append('trail')
        modelvalues.append(97.0 + (i % 5) * 1.0)

df_diurnal = pd.DataFrame({
    'value': values,
    'label': labels,
    'modelvalue': modelvalues
}, index=dates)
df_diurnal.index.name = 'datetime'

# Fill with diurnal debias
result_diurnal = fill_with_diurnal_debias(
    df_diurnal.copy(), 
    min_sample_size=2,
    max_value=100.0  # Constraint to physical maximum
)

gap_filled = result_diurnal.loc[result_diurnal.index[20:25]]
print("Gap-filled values (indices 20-24):")
print(gap_filled[['modelvalue', 'fillvalue']].to_string())
print()
print(f"Max filled value: {result_diurnal['fillvalue'].max():.2f}")
print(f"All values <= 100: {result_diurnal['fillvalue'].max() <= 100.0}")
print()

# Example 3: Practical usage tips
print("\nExample 3: Practical Usage Recommendations")
print("-" * 70)
print("""
For different observation types, consider these constraints:

1. Relative Humidity:
   - min_value=0.0, max_value=100.0
   - Prevents physically impossible values

2. Wind Speed:
   - min_value=0.0
   - Wind speed cannot be negative

3. Solar Radiation:
   - min_value=0.0
   - Radiation cannot be negative

4. Temperature:
   - min_value=-50.0, max_value=50.0 (or appropriate for your region)
   - Prevents extreme outliers from gap-filling errors

5. Precipitation:
   - min_value=0.0
   - Precipitation cannot be negative
""")

print("=" * 70)
print("Example completed successfully!")
print("=" * 70)
