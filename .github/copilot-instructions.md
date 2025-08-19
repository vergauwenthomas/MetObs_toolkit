# MetObs-toolkit AI Agent Instructions

## Architecture Overview

This is a scientific Python package for processing meteorological observations with a hierarchical data model:

- **Dataset**: Main container class holding multiple weather stations and their observations
- **Station**: Individual weather station with location metadata and time series data  
- **SensorData**: Single sensor measurements (e.g., temperature, humidity) for one station
- **Site**: Spatial metadata and Google Earth Engine (GEE) integration for each station
- **ModelDataset**: Gridded model output wrapper around xarray.Dataset for verification

Key entry points: `Dataset`, `Analysis`, `Verification`, `ModelDataset` classes from `src/metobs_toolkit/`

## Critical Workflows

### Data Import Pattern
All raw data flows through template-based parsing:
```python
dataset = metobs_toolkit.Dataset()
dataset.import_data_from_file(
    template_file=template_path,  # JSON describing data structure
    input_data_file=data_path,
    input_metadata_file=metadata_path
)
```

### Development Pipeline
Use `deploiment/develop_pipeline.sh` for full development workflow:
- Poetry dependency management and building
- Black code formatting (config in pyproject.toml)
- pytest with matplotlib image comparisons (`--mpl`)
- Documentation building with Sphinx
- Notebook testing with `--nbval-lax`

### Testing Architecture
- **Solution-based testing**: Tests store expected outputs as pickled objects in `tests/pkled_solutions/`
- **Image comparison**: Plotting tests use `@pytest.mark.mpl_image_compare` decorator
- **Notebook testing**: Example notebooks in `docs/examples/` and `docs/topics/` run as tests
- CI runs individual test files in parallel matrix strategy

## Project-Specific Patterns

### Google Earth Engine Integration
All GEE operations require authentication setup first:
```python
metobs_toolkit.connect_to_gee()  # One-time setup per session
dataset.get_LCZ()  # Extract Local Climate Zones
dataset.get_altitude()  # Extract elevation data
```

GEE data managers in `geedatasetmanagers.py` handle static (DEM, landcover) vs dynamic (weather model) datasets.

### Quality Control Pipeline
QC methods modify data in-place and track outliers:
```python
dataset.gross_value_check(target_obstype="temp", lower_threshold=-15, upper_threshold=35)
dataset.persistence_check(target_obstype="temp", timewindow="1h")
# Access results via dataset.outliersdf
```

### Logging Decorator Pattern
All public methods must use `@log_entry` decorator from `backend_collection.loggingmodule`:
```python
from metobs_toolkit.backend_collection.loggingmodule import log_entry

@log_entry
def my_method(self, arg1, arg2):
    """Method docstring required."""
    pass
```

### Collection-Based Organization
Code organized in `*_collection/` subdirectories by functionality:
- `qc_collection/`: Quality control algorithms
- `plot_collection/`: Visualization functions  
- `backend_collection/`: Utilities, error classes, data constructors
- `verif_collection/`: Model verification and scoring metrics

## Key Integrations

### xarray/Dask for Model Data
`ModelDataset` class wraps xarray.Dataset with required coordinates: `['lat', 'lon', 'validtime', 'reference_time']`. Field definitions in `nwp_collection/field_defenitions.py` map model variables to standard names/units.

### Property-Based DataFrames
Main data access through computed properties that construct DataFrames on-demand:
- `dataset.df` - main observations
- `dataset.metadf` - station metadata  
- `dataset.outliersdf` - QC flagged data
- `dataset.modeldatadf` - model timeseries

### Multiprocessing Support
QC operations support `use_mp=True` parameter for parallel processing across stations.

## Documentation Requirements

When adding new public methods:
1. Add numpy-style docstring
2. Include in appropriate `docs/reference/api/` file  
3. Add test in `tests/` following solution-based pattern
4. Apply `@log_entry` decorator
5. Update property docstrings with `@copy_doc` if building derived data
