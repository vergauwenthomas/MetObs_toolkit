# Type Hint Only Imports in MetObs Toolkit

This document summarizes modules in the `metobs_toolkit` package where classes are imported at the top level but are **only used for type hints** (return type annotations, parameter annotations, or docstring references) and have no functional use at runtime.

These imports could potentially be moved inside `if TYPE_CHECKING:` blocks to reduce import overhead and avoid unnecessary dependencies at runtime.

---

## Summary Table

| Module | Import | Source | Usage |
|--------|--------|--------|-------|
| `dataset.py` | `Axes` | `matplotlib.pyplot` | Return type annotation |
| `dataset.py` | `xrDataset` | `xarray` | Return type annotation |
| `station.py` | `Axes` | `matplotlib.pyplot` | Return type annotation |
| `station.py` | `xrDataset` | `xarray` | Return type annotation |
| `modeltimeseries.py` | `Axes` | `matplotlib.pyplot` | Return type annotation |
| `modeltimeseries.py` | `xrDataset` | `xarray` | Return type annotation |
| `sensordata.py` | `Axes` | `matplotlib.pyplot` | Return type annotation |
| `sensordata.py` | `xrDataset` | `xarray` | Return type annotation |

---

## Detailed Analysis by Module

### `dataset.py`

**Imports for type hints only:**
```python
from matplotlib.pyplot import Axes
from xarray import Dataset as xrDataset
```

**Usage:**
- `Axes`: Used as return type in `make_timeseries_plot()` (line 1087) and `make_geo_obs_plot()` (line 1219)
- `xrDataset`: Used as return type in `to_xr()` (line 375)

**Note:** These classes are never instantiated directly in this module. `Axes` objects are created via `plotting.create_axes()`, and xarray Datasets are created via `dataset_to_xr()`.

---

### `station.py`

**Imports for type hints only:**
```python
from matplotlib.pyplot import Axes
from xarray import Dataset as xrDataset
```

**Usage:**
- `Axes`: Used as return type in `make_timeseries_plot()` (line 1676) and `make_geo_obs_plot()` (line 1804)
- `xrDataset`: Used as return type in `to_xr()` (line 219)

**Note:** Same pattern as `dataset.py`. Axes created via helper functions, xarray Datasets via conversion functions.

---

### `modeltimeseries.py`

**Imports for type hints only:**
```python
from matplotlib.pyplot import Axes
from xarray import Dataset as xrDataset
```

**Usage:**
- `Axes`: Used as return type in `pd_plot()` (line 275) and `make_plot()` (line 285), and as parameter type annotation (line 293)
- `xrDataset`: Used as return type in `to_xr()` (line 222)

---

### `sensordata.py` ✅ (Already Optimized)

This module already uses the `TYPE_CHECKING` pattern correctly:

```python
from __future__ import annotations
from typing import Literal, Union, TYPE_CHECKING

# ...runtime imports...

if TYPE_CHECKING:
    import pytz
    import dateutil.tz
    import datetime
```

**Remaining type-hint-only imports that could be moved:**
```python
from matplotlib.pyplot import Axes
from xarray import Dataset as xrDataset
```

**Usage:**
- `Axes`: Used as return type in `pd_plot()` (line 792)
- `xrDataset`: Used as return type in `to_xr()` (line 243)

---

## Recommended Refactoring Pattern

For each module, apply the following pattern:

```python
from __future__ import annotations
from typing import TYPE_CHECKING

# ... other runtime imports ...

if TYPE_CHECKING:
    from matplotlib.pyplot import Axes
    from xarray import Dataset as xrDataset
```

### Benefits:
1. **Faster import times** - Heavy libraries like matplotlib and xarray won't be loaded just for type annotations
2. **Reduced memory footprint** - Especially important if only using a subset of the toolkit
3. **Cleaner dependency graph** - Makes it clear which imports are for runtime vs. static analysis
4. **Forward reference support** - Using `from __future__ import annotations` allows string-based annotations

### Caveat:
If any of these types are used in `isinstance()` checks or as base classes, they must remain as runtime imports. Based on the analysis, this is not the case for `Axes` or `xrDataset` in the listed modules.

---

## Modules NOT Requiring Changes

The following modules import classes that ARE used functionally (not just for type hints):

| Module | Import | Reason |
|--------|--------|--------|
| `analysis.py` | `Dataset`, `Station` | Used functionally in class methods |
| `gap.py` | `ModelTimeSeries`, `Obstype` | Used functionally |
| `site.py` | `GEEStaticDatasetManager` | Used functionally |
| `geedatasetmanagers.py` | `PathLike` | Used in `isinstance()` checks or path handling |
| `template.py` | All imports | All used functionally |

---

## Implementation Priority

1. **High Priority**: `dataset.py`, `station.py` - Most commonly imported modules
2. **Medium Priority**: `modeltimeseries.py` - Less frequently used directly
3. **Low Priority**: `sensordata.py` - Already partially optimized

---

*Generated: December 2025*
