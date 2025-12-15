# MetObs Toolkit - Argument Naming Inconsistencies

This document lists argument names that are inconsistent across different methods and functions in the metobs_toolkit package.

---

## 1. Observation Type Naming

**Issue**: Multiple different names are used for the observation type argument.

| Argument Name | Used In |
|---------------|---------|
| `target_obstype` | `Dataset.sync_records`, `Dataset.resample`, `Dataset.gross_value_check`, `Dataset.persistence_check`, `Dataset.repetitions_check`, `Dataset.step_check`, `Dataset.window_variation_check`, `Dataset.buddy_check`, `Dataset.buddy_check_with_safetynets`, `Dataset.get_qc_stats`, `Dataset.interpolate_gaps`, `Dataset.fill_gaps_with_*`, `Station.resample`, `Station.gross_value_check`, `Station.persistence_check`, `Station.repetitions_check`, `Station.step_check`, `Station.window_variation_check`, `Station.get_qc_stats`, `Station.fill_gaps_with_*`, `Station.interpolate_gaps` |
| `trgobstype` | `Analysis.aggregate_df`, `Analysis.plot_diurnal_cycle`, `Analysis.plot_diurnal_cycle_with_reference_station` |
| `obstype` | `Dataset.add_new_observationtype`, `Dataset.make_plot_of_modeldata`, `Dataset.make_plot`, `Dataset.convert_outliers_to_gaps`, `Station.get_sensor`, `Station.get_modeltimeseries`, `Station.convert_outliers_to_gaps`, `Station.make_plot_of_modeldata`, `Station.make_plot`, `SensorData.__init__`, `Gap.__init__` |
| `trg_obstype` | `WhiteSet.create_sensorwhitelist` |
| `obsname` | `Obstype.__init__` |

**Recommendation**: Standardize to `target_obstype` for methods that operate on a specific obstype, or `obstype` for methods that refer to an obstype object/name.

---

## 2. GEE Dataset Manager Naming

**Issue**: The parameter name for GEE dataset managers varies inconsistently.

| Argument Name | Used In |
|---------------|---------|
| `geestaticdatasetmanager` | `Dataset.get_static_gee_point_data`, `Dataset.get_static_gee_buffer_fraction_data`, `Station.get_static_gee_point_data` |
| `geestaticdataset` | `Station.get_static_gee_buffer_fraction_data` |
| `geedynamicdatasetmanager` | `Dataset.import_gee_data_from_file`, `Dataset.get_gee_timeseries_data`, `Station.get_gee_timeseries_data` |
| `geedatasetmanager` | `Dataset.make_gee_plot` |

**Recommendation**: Standardize to `geestaticdatasetmanager` and `geedynamicdatasetmanager` (or shorter `gee_static_manager` and `gee_dynamic_manager`).

---

## 3. File Path Arguments

**Issue**: Different naming conventions for file path parameters.

| Argument Name | Used In |
|---------------|---------|
| `filepath` | `Dataset.to_netcdf`, `Dataset.import_gee_data_from_file`, `Station.to_netcdf` |
| `target_file` | `Dataset.to_parquet`, `Dataset.to_csv`, `Station.to_parquet`, `Station.to_csv` |
| `outputfolder` + `filename` | `Dataset.make_gee_plot`, `Dataset.save_dataset_to_pkl`, `GEEStaticDatasetManager.make_gee_plot`, `GEEDynamicDatasetManager.make_gee_plot` |
| `target_folder` | `Dataset.save_dataset_to_pkl` |

**Recommendation**: Standardize to either `filepath` for full paths or `target_folder`/`target_file` pattern.

---

## 4. Datetime Arguments

**Issue**: Different naming for start/end datetime parameters.

| Argument Name | Used In |
|---------------|---------|
| `startdt`, `enddt` | `Analysis.subset_period` |
| `startdt_utc`, `enddt_utc` | `Dataset.get_gee_timeseries_data`, `Station.get_gee_timeseries_data`, `GEEDynamicDatasetManager.extract_timeseries_data` |

**Recommendation**: Use `startdt_utc`/`enddt_utc` consistently when UTC is required, or `startdt`/`enddt` with clear documentation about timezone handling.

---

## 5. Overwrite/Update Parameters

**Issue**: Different naming for parameters controlling overwrite behavior.

| Argument Name | Used In |
|---------------|---------|
| `overwrite` | `Dataset.save_dataset_to_pkl`, `Dataset.make_gee_plot`, `Dataset.get_static_gee_point_data`, `Dataset.get_static_gee_buffer_fraction_data`, `Dataset.get_LCZ`, `Dataset.get_altitude`, `Station.get_*`, `GEEStaticDatasetManager.make_gee_plot`, `GEEDynamicDatasetManager.make_gee_plot` |
| `overwrite_fill` | `Dataset.interpolate_gaps`, `Dataset.fill_gaps_with_*`, `Station.fill_gaps_with_*`, `Station.interpolate_gaps`, `SensorData.fill_gap_with_modeldata`, `SensorData.interpolate_gaps` |
| `force_update` | `Dataset.import_gee_data_from_file`, `Station.add_to_sensordata`, `Station.add_to_modeldata` |
| `update_stations` | `Dataset.get_landcover_fractions` |

**Recommendation**: Standardize to `overwrite` for general overwrite behavior and `overwrite_fill` specifically for gap filling operations.

---

## 6. Station Reference Arguments

**Issue**: Different naming for station reference parameters.

| Argument Name | Used In |
|---------------|---------|
| `ref_station` | `Analysis.plot_diurnal_cycle_with_reference_station` |
| `trg_station` | `WhiteSet.create_sensorwhitelist` |
| `stationname` | `Dataset.get_station`, `Station.__init__`, `SensorData.__init__`, `Gap.__init__` |

**Recommendation**: Use `stationname` for station identifiers, `ref_station` for reference stations in comparisons.

---

## 7. Color Arguments in Plotting

**Issue**: Different naming for color-related parameters.

| Argument Name | Used In |
|---------------|---------|
| `colorby` | `Dataset.make_plot`, `Station.make_plot`, `Analysis.plot_diurnal_cycle`, `Analysis.plot_diurnal_cycle_with_reference_station` |
| `linecolor` | `Station.make_plot_of_modeldata`, `Station.make_plot`, `ModelTimeSeries.make_plot` |
| `colormap` | `Dataset.make_plot_of_modeldata` |
| `colordict` | `Analysis.plot_diurnal_cycle`, `Analysis.plot_diurnal_cycle_with_reference_station` |

**Recommendation**: These serve different purposes (categorical coloring vs line colors vs colormaps), but naming could be more consistent. Consider `color` for single line color, `colorby` for categorical, `cmap` for colormaps.

---

## 8. Missing Parameters Across Similar Methods

### `initialize_gee` parameter inconsistency

| Method | Has `initialize_gee` |
|--------|---------------------|
| `Dataset.get_LCZ` | ✅ Yes |
| `Dataset.get_altitude` | ✅ Yes |
| `Dataset.get_landcover_fractions` | ✅ Yes |
| `Dataset.get_static_gee_point_data` | ✅ Yes |
| `Dataset.get_static_gee_buffer_fraction_data` | ✅ Yes |
| `Dataset.get_gee_timeseries_data` | ❌ No |
| `Dataset.import_gee_data_from_file` | ❌ No |
| `Dataset.make_gee_plot` | ❌ No |
| `Station.get_LCZ` | ✅ Yes |
| `Station.get_altitude` | ✅ Yes |
| `Station.get_landcover_fractions` | ❌ No |
| `Station.get_static_gee_point_data` | ✅ Yes |
| `Station.get_static_gee_buffer_fraction_data` | ✅ Yes |
| `Station.get_gee_timeseries_data` | ❌ No |

**Recommendation**: Add `initialize_gee` parameter to all GEE-related methods for consistency.

---

## 9. Model-related Arguments

**Issue**: `modelname` and `modelvariable` are present in `Station` gap fill methods but not in `Dataset` equivalents.

| Method | Has `modelname`/`modelvariable` |
|--------|--------------------------------|
| `Dataset.fill_gaps_with_raw_modeldata` | ❌ No |
| `Dataset.fill_gaps_with_debiased_modeldata` | ❌ No |
| `Dataset.fill_gaps_with_diurnal_debiased_modeldata` | ❌ No |
| `Dataset.fill_gaps_with_weighted_diurnal_debiased_modeldata` | ❌ No |
| `Station.fill_gaps_with_raw_modeldata` | ✅ Yes |
| `Station.fill_gaps_with_debiased_modeldata` | ✅ Yes |
| `Station.fill_gaps_with_diurnal_debiased_modeldata` | ✅ Yes |
| `Station.fill_gaps_with_weighted_diurnal_debiased_modeldata` | ✅ Yes |

**Recommendation**: Consider adding `modelname`/`modelvariable` to Dataset methods for filtering which model data to use.

---

## 10. QC Method Arguments

**Issue**: `use_mp` (multiprocessing) is only available in `Dataset` QC methods, not `Station`.

| Method | Has `use_mp` |
|--------|-------------|
| `Dataset.gross_value_check` | ✅ Yes |
| `Dataset.persistence_check` | ✅ Yes |
| `Dataset.repetitions_check` | ✅ Yes |
| `Dataset.step_check` | ✅ Yes |
| `Dataset.window_variation_check` | ✅ Yes |
| `Dataset.buddy_check` | ✅ Yes |
| `Dataset.buddy_check_with_safetynets` | ✅ Yes |
| `Station.*_check` | ❌ No |

**Note**: This is intentional as `use_mp` is for parallelizing across stations and doesn't apply to single Station objects.

---

## Summary of Recommended Changes

1. **Obstype naming**: Standardize `trgobstype` → `target_obstype` in `Analysis` class
2. **GEE managers**: Fix `geestaticdataset` → `geestaticdatasetmanager` in `Station.get_static_gee_buffer_fraction_data`
3. **Datetime**: Consider standardizing `startdt`/`enddt` vs `startdt_utc`/`enddt_utc`
4. **Station refs**: Standardize `trg_station` → `stationname` or `target_station` in `WhiteSet`
5. **GEE init**: Add `initialize_gee` to missing GEE methods for consistency
6. **Color args**: Consider renaming `colormap` to `cmap` for consistency with matplotlib conventions
