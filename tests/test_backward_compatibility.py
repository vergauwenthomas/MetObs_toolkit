"""Backward compatibility tests for the MetObs toolkit.

These tests verify that the current API remains compatible with usage patterns
from previous versions. They cover:

* Data import from all supported file formats and layout structures.
* Method signatures: all documented keyword argument names must still exist.
* Deprecated methods: must raise ``DeprecationWarning``, not silently disappear.
* Callable invocations: public methods must accept their documented keyword
  arguments without raising ``TypeError``.
* Pickle round-trip: a dataset saved to ``.pkl`` must be loadable and usable.
"""

import inspect
import sys
import tempfile
import pickle
from pathlib import Path

import pytest
import pandas as pd

# ---------------------------------------------------------------------------
# Path setup – ensure the local source tree is used
# ---------------------------------------------------------------------------
_LIBFOLDER = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_LIBFOLDER / "src"))

import metobs_toolkit
from solutionclass import datadir

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_DEMO_TEMPLATE = metobs_toolkit.demo_template
_DEMO_DATA = metobs_toolkit.demo_datafile
_DEMO_META = metobs_toolkit.demo_metadatafile

_SINGLE_STATION_TEMPLATE = datadir / "single_station_template.json"
_SINGLE_STATION_DATA = datadir / "single_station.csv"
_SINGLE_STATION_META = datadir / "single_station_metadata.csv"

_WIDE_TEMPLATE = datadir / "wide_test_template.json"
_WIDE_DATA = datadir / "wide_test_data.csv"

_PARQUET_TEMPLATE = datadir / "demo_template_for_parquet.json"
_PARQUET_DATA_UTC = datadir / "demo_data_with_timezone_utc.parquet"
_PARQUET_DATA_PARIS = datadir / "demo_data_with_timezone_paris.parquet"

_BASELINE_PKL_DIR = _LIBFOLDER / "tests" / "baseline_pickle_objects"
_LEGACY_PKL = _BASELINE_PKL_DIR / "dataset_created_with_and_for_v1.x.x.pkl"


# %%


def _make_demo_dataset() -> metobs_toolkit.Dataset:
    """Return a small, fully-loaded demo dataset (multi-station, long CSV)."""
    ds = metobs_toolkit.Dataset()
    ds.import_data_from_file(
        template_file=_DEMO_TEMPLATE,
        input_data_file=_DEMO_DATA,
        input_metadata_file=_DEMO_META,
    )
    return ds


def _param_names(method) -> set:
    """Return the set of parameter names (excluding 'self') for a callable."""
    return set(inspect.signature(method).parameters.keys()) - {"self"}


# ===========================================================================
# 1.  DATA IMPORT BACKWARD COMPATIBILITY
# ===========================================================================


class TestImportDataBackwardCompat:
    """Verify that all supported data formats and layouts still import cleanly."""

    def test_import_long_csv_with_metadata(self):
        """Long-format CSV with separate metadata file – the canonical use-case."""
        ds = metobs_toolkit.Dataset()
        ds.import_data_from_file(
            template_file=_DEMO_TEMPLATE,
            input_data_file=_DEMO_DATA,
            input_metadata_file=_DEMO_META,
        )
        assert len(ds.stations) > 0
        assert not ds.df.empty

    def test_import_long_csv_data_only(self):
        """Long-format CSV without a metadata file must still produce stations."""
        ds = metobs_toolkit.Dataset()
        ds.import_data_from_file(
            template_file=_DEMO_TEMPLATE,
            input_data_file=_DEMO_DATA,
        )
        assert len(ds.stations) > 0
        assert not ds.df.empty

    def test_import_metadata_only(self):
        """Metadata-only import (no data file) – used for spatial-only workflows."""
        ds = metobs_toolkit.Dataset()
        ds.import_data_from_file(
            template_file=_DEMO_TEMPLATE,
            input_metadata_file=_DEMO_META,
        )
        assert len(ds.stations) > 0
        assert not ds.metadf.empty

    def test_import_wide_csv(self):
        """Wide-format CSV (one column per station) must be parseable."""
        ds = metobs_toolkit.Dataset()
        ds.import_data_from_file(
            template_file=_WIDE_TEMPLATE,
            input_data_file=_WIDE_DATA,
        )
        assert len(ds.stations) > 0
        assert not ds.df.empty

    def test_import_single_station_csv_with_metadata(self):
        """Single-station CSV with matching metadata file."""
        ds = metobs_toolkit.Dataset()
        ds.import_data_from_file(
            template_file=_SINGLE_STATION_TEMPLATE,
            input_data_file=_SINGLE_STATION_DATA,
            input_metadata_file=_SINGLE_STATION_META,
        )
        assert len(ds.stations) == 1
        assert not ds.df.empty

    def test_import_single_station_csv_without_metadata(self):
        """Single-station CSV without metadata – station name comes from template."""
        ds = metobs_toolkit.Dataset()
        ds.import_data_from_file(
            template_file=_SINGLE_STATION_TEMPLATE,
            input_data_file=_SINGLE_STATION_DATA,
        )
        assert len(ds.stations) == 1

    def test_import_parquet_utc(self):
        """Parquet file with UTC timezone must be importable."""
        ds = metobs_toolkit.Dataset()
        ds.import_data_from_file(
            template_file=_PARQUET_TEMPLATE,
            input_data_file=_PARQUET_DATA_UTC,
        )
        assert len(ds.stations) > 0
        assert not ds.df.empty

    def test_import_parquet_paris_timezone_mismatch_raises(self):
        """Parquet with Europe/Paris timezone against a UTC template must raise
        MetObsTemplateError – this error path must remain backward-compatible."""
        from metobs_toolkit.backend_collection.errorclasses import MetObsTemplateError

        ds = metobs_toolkit.Dataset()
        with pytest.raises(MetObsTemplateError):
            ds.import_data_from_file(
                template_file=_PARQUET_TEMPLATE,
                input_data_file=_PARQUET_DATA_PARIS,
            )

    def test_import_with_freq_estimation_method_highest(self):
        """The 'highest' freq_estimation_method must still be accepted."""
        ds = metobs_toolkit.Dataset()
        ds.import_data_from_file(
            template_file=_DEMO_TEMPLATE,
            input_data_file=_DEMO_DATA,
            freq_estimation_method="highest",
        )
        assert len(ds.stations) > 0

    def test_import_with_explicit_tolerances(self):
        """Explicit string-form tolerances must still be accepted."""
        ds = metobs_toolkit.Dataset()
        ds.import_data_from_file(
            template_file=_DEMO_TEMPLATE,
            input_data_file=_DEMO_DATA,
            freq_estimation_simplify_tolerance="3min",
            origin_simplify_tolerance="10min",
            timestamp_tolerance="5min",
        )
        assert len(ds.stations) > 0

    def test_import_with_kwargs_data_read(self):
        """Extra CSV reader kwargs forwarded via kwargs_data_read must be accepted."""
        ds = metobs_toolkit.Dataset()
        # The demo CSV uses ';' as separator – pass it explicitly
        ds.import_data_from_file(
            template_file=_DEMO_TEMPLATE,
            input_data_file=_DEMO_DATA,
            kwargs_data_read={"sep": ";"},
        )
        assert len(ds.stations) > 0

    def test_import_with_timedelta_tolerances(self):
        """pd.Timedelta objects must be accepted for tolerance arguments."""
        ds = metobs_toolkit.Dataset()
        ds.import_data_from_file(
            template_file=_DEMO_TEMPLATE,
            input_data_file=_DEMO_DATA,
            freq_estimation_simplify_tolerance=pd.Timedelta("2min"),
            origin_simplify_tolerance=pd.Timedelta("5min"),
            timestamp_tolerance=pd.Timedelta("4min"),
        )
        assert len(ds.stations) > 0


# ===========================================================================
# 2.  PICKLE ROUND-TRIP (import_dataset_from_pkl)
# ===========================================================================


class TestPickleRoundTrip:
    """A dataset saved to .pkl with the current version must reload correctly."""

    def test_save_and_load_pkl(self, tmp_path):
        ds = _make_demo_dataset()
        pkl_path = tmp_path / "saved_dataset.pkl"
        ds.save_dataset_to_pkl(filepath=str(pkl_path), overwrite=True)
        loaded = metobs_toolkit.import_dataset_from_pkl(target_path=str(pkl_path))
        assert isinstance(loaded, metobs_toolkit.Dataset)
        assert len(loaded.stations) == len(ds.stations)

    def test_load_pkl_stations_intact(self, tmp_path):
        ds = _make_demo_dataset()
        pkl_path = tmp_path / "saved.pkl"
        ds.save_dataset_to_pkl(filepath=str(pkl_path), overwrite=True)
        loaded = metobs_toolkit.import_dataset_from_pkl(target_path=str(pkl_path))
        assert not loaded.df.empty
        assert not loaded.metadf.empty

    def test_load_legacy_pkl_from_v1(self):
        """A Dataset pickled with an older v1.x.x release must still load correctly.

        The file ``tests/baseline_pickle_objects/dataset_created_with_and_for_v1.x.x.pkl``
        was generated with a previous MetObs version and must remain loadable by the
        current version without errors.
        """
        assert _LEGACY_PKL.is_file(), (
            f"Baseline pickle file not found: {_LEGACY_PKL}\n"
            "This file must be present to run the legacy backward-compatibility test."
        )
        loaded = metobs_toolkit.import_dataset_from_pkl(target_path=str(_LEGACY_PKL))
        assert isinstance(
            loaded, metobs_toolkit.Dataset
        ), "Loading the legacy v1.x.x pickle did not return a Dataset instance."
        assert len(loaded.stations) > 0, "Legacy Dataset loaded without any stations."
        assert (
            not loaded.df.empty
        ), "Legacy Dataset loaded but observations DataFrame is empty."
        assert (
            not loaded.metadf.empty
        ), "Legacy Dataset loaded but metadata DataFrame is empty."
        assert (
            not loaded.outliersdf.empty
        ), "Legacy Dataset loaded but outliers DataFrame is empty."
        assert (
            not loaded.modeldatadf.empty
        ), "Legacy Dataset loaded but model data DataFrame is empty."
        assert (
            not loaded.gapsdf.empty
        ), "Legacy Dataset loaded but gaps DataFrame is empty."


# ===========================================================================
# 3.  METHOD SIGNATURE BACKWARD COMPATIBILITY
# ===========================================================================
#
# Each test checks that a known set of parameter names (as they appeared in
# previous releases) is still present in the current method signature.
# Adding new *optional* parameters is allowed; removing or renaming existing
# ones is a breaking change.
# ---------------------------------------------------------------------------


class TestSignatureBackwardCompat:
    """Parameter names that existed in previous versions must still be present."""

    # ---- Dataset.import_data_from_file ----

    def test_import_data_from_file_params(self):
        expected = {
            "template_file",
            "input_data_file",
            "input_metadata_file",
            "freq_estimation_method",
            "freq_estimation_simplify_tolerance",
            "origin_simplify_tolerance",
            "timestamp_tolerance",
            "kwargs_data_read",
            "kwargs_metadata_read",
            "templatefile_is_url",
        }
        actual = _param_names(metobs_toolkit.Dataset.import_data_from_file)
        missing = expected - actual
        assert not missing, f"Missing parameters in import_data_from_file: {missing}"

    # ---- Dataset.subset_by_stations ----

    def test_subset_by_stations_params(self):
        expected = {"stationnames"}
        actual = _param_names(metobs_toolkit.Dataset.subset_by_stations)
        missing = expected - actual
        assert not missing, f"Missing parameters in subset_by_stations: {missing}"

    # ---- Dataset.get_station ----

    def test_get_station_params(self):
        expected = {"stationname"}
        actual = _param_names(metobs_toolkit.Dataset.get_station)
        missing = expected - actual
        assert not missing, f"Missing parameters in get_station: {missing}"

    # ---- Dataset.get_info ----

    def test_get_info_params(self):
        expected = {"printout"}
        actual = _param_names(metobs_toolkit.Dataset.get_info)
        missing = expected - actual
        assert not missing, f"Missing parameters in Dataset.get_info: {missing}"

    # ---- Dataset.copy ----

    def test_copy_params(self):
        expected = {"deep"}
        actual = _param_names(metobs_toolkit.Dataset.copy)
        missing = expected - actual
        assert not missing, f"Missing parameters in Dataset.copy: {missing}"

    # ---- Dataset.sync_records ----

    def test_sync_records_params(self):
        expected = {"obstype", "timestamp_shift_tolerance", "freq_shift_tolerance"}
        actual = _param_names(metobs_toolkit.Dataset.sync_records)
        missing = expected - actual
        assert not missing, f"Missing parameters in sync_records: {missing}"

    # ---- Dataset.resample ----

    def test_resample_params(self):
        expected = {"target_freq", "obstype", "shift_tolerance", "origin"}
        actual = _param_names(metobs_toolkit.Dataset.resample)
        missing = expected - actual
        assert not missing, f"Missing parameters in Dataset.resample: {missing}"

    # ---- Dataset.rename_stations ----

    def test_rename_stations_params(self):
        expected = {"renamedict"}
        actual = _param_names(metobs_toolkit.Dataset.rename_stations)
        missing = expected - actual
        assert not missing, f"Missing parameters in rename_stations: {missing}"

    # ---- Dataset.save_dataset_to_pkl ----

    def test_save_dataset_to_pkl_params(self):
        expected = {"filepath", "overwrite"}
        actual = _param_names(metobs_toolkit.Dataset.save_dataset_to_pkl)
        missing = expected - actual
        assert not missing, f"Missing parameters in save_dataset_to_pkl: {missing}"

    # ---- Dataset.add_new_observationtype ----

    def test_add_new_observationtype_params(self):
        expected = {"obstype"}
        actual = _param_names(metobs_toolkit.Dataset.add_new_observationtype)
        missing = expected - actual
        assert not missing, f"Missing parameters in add_new_observationtype: {missing}"

    # ---- QC method signatures ----

    def test_gross_value_check_params(self):
        expected = {"obstype", "lower_threshold", "upper_threshold", "whiteset"}
        actual = _param_names(metobs_toolkit.Dataset.gross_value_check)
        missing = expected - actual
        assert not missing, f"Missing parameters in gross_value_check: {missing}"

    def test_persistence_check_params(self):
        expected = {"obstype", "timewindow", "whiteset"}
        actual = _param_names(metobs_toolkit.Dataset.persistence_check)
        missing = expected - actual
        assert not missing, f"Missing parameters in persistence_check: {missing}"

    def test_repetitions_check_params(self):
        expected = {"obstype", "whiteset"}
        actual = _param_names(metobs_toolkit.Dataset.repetitions_check)
        missing = expected - actual
        assert not missing, f"Missing parameters in repetitions_check: {missing}"

    def test_step_check_params(self):
        expected = {
            "obstype",
            "max_increase_per_second",
            "max_decrease_per_second",
            "whiteset",
        }
        actual = _param_names(metobs_toolkit.Dataset.step_check)
        missing = expected - actual
        assert not missing, f"Missing parameters in step_check: {missing}"

    def test_window_variation_check_params(self):
        expected = {
            "obstype",
            "timewindow",
            "min_records_per_window",
            "whiteset",
        }
        actual = _param_names(metobs_toolkit.Dataset.window_variation_check)
        missing = expected - actual
        assert not missing, f"Missing parameters in window_variation_check: {missing}"

    def test_buddy_check_params(self):
        expected = {
            "obstype",
            "spatial_buddy_radius",
            "min_sample_size",
            "min_std",
            "spatial_z_threshold",
            "N_iter",
            "instantaneous_tolerance",
            "lapserate",
            "whiteset",
        }
        actual = _param_names(metobs_toolkit.Dataset.buddy_check)
        missing = expected - actual
        assert not missing, f"Missing parameters in buddy_check: {missing}"

    def test_buddy_check_with_safetynets_params(self):
        expected = {
            "obstype",
            "spatial_buddy_radius",
            "safety_net_configs",
            "min_sample_size",
            "min_std",
            "spatial_z_threshold",
            "N_iter",
            "instantaneous_tolerance",
            "whiteset",
        }
        actual = _param_names(metobs_toolkit.Dataset.buddy_check_with_safetynets)
        missing = expected - actual
        assert (
            not missing
        ), f"Missing parameters in buddy_check_with_safetynets: {missing}"

    # ---- Gap handling ----

    def test_convert_outliers_to_gaps_params(self):
        expected = {"obstype"}
        actual = _param_names(metobs_toolkit.Dataset.convert_outliers_to_gaps)
        missing = expected - actual
        assert not missing, f"Missing parameters in convert_outliers_to_gaps: {missing}"

    def test_interpolate_gaps_params(self):
        expected = {"obstype", "method", "max_gap_duration_to_fill", "overwrite_fill"}
        actual = _param_names(metobs_toolkit.Dataset.interpolate_gaps)
        missing = expected - actual
        assert not missing, f"Missing parameters in interpolate_gaps: {missing}"

    # ---- QC stats ----

    def test_get_qc_stats_params(self):
        expected = {"obstype", "make_plot"}
        actual = _param_names(metobs_toolkit.Dataset.get_qc_stats)
        missing = expected - actual
        assert not missing, f"Missing parameters in get_qc_stats: {missing}"

    # ---- Obstype class ----

    def test_obstype_constructor_params(self):
        expected = {"name", "std_unit", "description"}
        actual = _param_names(metobs_toolkit.Obstype.__init__)
        missing = expected - actual
        assert not missing, f"Missing parameters in Obstype.__init__: {missing}"

    # ---- module-level functions ----

    def test_import_dataset_from_pkl_params(self):
        expected = {"target_path"}
        actual = _param_names(metobs_toolkit.import_dataset_from_pkl)
        missing = expected - actual
        assert not missing, f"Missing parameters in import_dataset_from_pkl: {missing}"


# ===========================================================================
# 4.  DEPRECATED METHODS
# ===========================================================================


class TestDeprecatedMethods:
    """Deprecated public methods must raise DeprecationWarning (not disappear)."""

    def test_buddy_check_with_LCZ_safety_net_raises_deprecation(self):
        """buddy_check_with_LCZ_safety_net was replaced by buddy_check_with_safetynets."""
        ds = _make_demo_dataset()
        with pytest.raises(DeprecationWarning):
            ds.buddy_check_with_LCZ_safety_net()

    def test_buddy_check_with_LCZ_safety_net_still_exists(self):
        """The deprecated method must still be an attribute of Dataset."""
        assert hasattr(metobs_toolkit.Dataset, "buddy_check_with_LCZ_safety_net")


# ===========================================================================
# 5.  KEYWORD ARGUMENT INVOCATION
# ===========================================================================
#
# Call public methods with their documented keyword arguments (as a previous-
# version user would) and assert no TypeError is raised. These tests exercise
# the entire call path, not just signature inspection.
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def demo_dataset():
    """Module-scoped dataset to avoid repeating the slow import."""
    return _make_demo_dataset()


@pytest.fixture(scope="module")
def qced_dataset():
    """Module-scoped dataset with QC applied, for gap/fill tests."""
    ds = _make_demo_dataset()
    ds.gross_value_check(obstype="temp", lower_threshold=-15.0, upper_threshold=39.0)
    ds.convert_outliers_to_gaps(obstype="temp")
    return ds


class TestKeywordArgumentInvocations:
    """Call methods with all documented kwargs to detect renamed/removed params."""

    def test_dataset_get_info_printout_false(self, demo_dataset):
        result = demo_dataset.get_info(printout=False)
        assert isinstance(result, str)

    def test_dataset_copy_deep_true(self, demo_dataset):
        copy = demo_dataset.copy(deep=True)
        assert isinstance(copy, metobs_toolkit.Dataset)
        assert len(copy.stations) == len(demo_dataset.stations)

    def test_dataset_copy_deep_false(self, demo_dataset):
        copy = demo_dataset.copy(deep=False)
        assert isinstance(copy, metobs_toolkit.Dataset)

    def test_subset_by_stations_with_stationnames(self, demo_dataset):
        names = [sta.name for sta in demo_dataset.stations[:3]]
        if len(names) < 2:
            pytest.skip("Need at least 2 stations for subset test.")
        subset = demo_dataset.subset_by_stations(stationnames=names)
        assert isinstance(subset, metobs_toolkit.Dataset)

    def test_get_station_by_stationname(self, demo_dataset):
        name = demo_dataset.stations[0].name
        station = demo_dataset.get_station(stationname=name)
        assert station.name == name

    def test_rename_stations_renamedict(self, demo_dataset):
        ds_copy = demo_dataset.copy(deep=True)
        original_name = ds_copy.stations[0].name
        ds_copy.rename_stations(renamedict={original_name: "renamed_station_compat"})
        assert any(sta.name == "renamed_station_compat" for sta in ds_copy.stations)

    def test_add_new_observationtype(self, demo_dataset):
        ds_copy = demo_dataset.copy(deep=True)
        new_obs = metobs_toolkit.Obstype(
            name="custom_compat_obs",
            std_unit="degC",
            description="test obstype for backward compat",
        )
        ds_copy.add_new_observationtype(obstype=new_obs)
        assert "custom_compat_obs" in ds_copy.obstypes

    def test_gross_value_check_kwargs(self, demo_dataset):
        ds_copy = demo_dataset.copy(deep=True)
        ds_copy.gross_value_check(
            obstype="temp",
            lower_threshold=-15.0,
            upper_threshold=39.0,
            use_mp=False,
        )

    def test_persistence_check_kwargs(self, demo_dataset):
        ds_copy = demo_dataset.copy(deep=True)
        ds_copy.persistence_check(
            obstype="temp",
            timewindow="1h",
            use_mp=False,
        )

    def test_repetitions_check_kwargs(self, demo_dataset):
        ds_copy = demo_dataset.copy(deep=True)
        ds_copy.repetitions_check(
            obstype="temp",
            use_mp=False,
        )

    def test_step_check_kwargs(self, demo_dataset):
        ds_copy = demo_dataset.copy(deep=True)
        ds_copy.step_check(
            obstype="temp",
            max_increase_per_second=8.0 / 3600.0,
            max_decrease_per_second=-10.0 / 3600.0,
            use_mp=False,
        )

    def test_window_variation_check_kwargs(self, demo_dataset):
        ds_copy = demo_dataset.copy(deep=True)
        ds_copy.window_variation_check(
            obstype="temp",
            timewindow="1h",
            use_mp=False,
        )

    def test_convert_outliers_to_gaps_obstype(self, demo_dataset):
        ds_copy = demo_dataset.copy(deep=True)
        ds_copy.gross_value_check(obstype="temp", use_mp=False)
        ds_copy.convert_outliers_to_gaps(obstype="temp")

    def test_get_qc_stats_kwargs(self, demo_dataset):
        ds_copy = demo_dataset.copy(deep=True)
        ds_copy.gross_value_check(obstype="temp", use_mp=False)
        with pytest.warns(DeprecationWarning, match="1.1.0"):
            result = ds_copy.get_qc_stats(obstype="temp", make_plot=False)
        assert isinstance(result, dict)

    def test_interpolate_gaps_kwargs(self, qced_dataset):
        ds_copy = qced_dataset.copy(deep=True)
        ds_copy.interpolate_gaps(
            obstype="temp",
            method="time",
            max_gap_duration_to_fill="3h",
        )

    def test_to_xr(self, demo_dataset):
        xrds = demo_dataset.to_xr()
        assert xrds is not None

    def test_to_parquet(self, demo_dataset, tmp_path):
        path = tmp_path / "compat_test.parquet"
        demo_dataset.to_parquet(filepath=str(path), overwrite=True)
        assert path.exists()

    def test_to_csv(self, demo_dataset, tmp_path):
        path = tmp_path / "compat_test.csv"
        demo_dataset.to_csv(filepath=str(path), overwrite=True)
        assert path.exists()

    def test_to_netcdf(self, demo_dataset, tmp_path):
        path = tmp_path / "compat_test.nc"
        demo_dataset.to_netcdf(filepath=str(path), overwrite=True)
        assert path.exists()

    def test_sync_records_obstype(self, demo_dataset):
        ds_copy = demo_dataset.copy(deep=True)
        ds_copy.sync_records(obstype="temp")

    def test_resample_target_freq(self, demo_dataset):
        ds_copy = demo_dataset.copy(deep=True)
        ds_copy.resample(target_freq="1h", obstype="temp")

    # ---- Station-level methods ----

    def test_station_get_info(self, demo_dataset):
        station = demo_dataset.stations[0]
        info = station.get_info()
        # get_info returns None (prints) or str
        assert info is None or isinstance(info, str)

    def test_station_get_sensor(self, demo_dataset):
        station = demo_dataset.stations[0]
        sensor = station.get_sensor("temp")
        assert sensor is not None

    def test_station_copy(self, demo_dataset):
        station = demo_dataset.stations[0]
        copy = station.copy(deep=True)
        assert copy.name == station.name

    # ---- Obstype construction ----

    def test_obstype_positional_construction(self):
        obs = metobs_toolkit.Obstype("my_obs", "degC", "A test observation")
        assert obs.name == "my_obs"

    def test_obstype_keyword_construction(self):
        obs = metobs_toolkit.Obstype(
            name="my_obs_kw",
            std_unit="degC",
            description="A test observation (kw)",
        )
        assert obs.name == "my_obs_kw"


# ===========================================================================
# 6.  PUBLIC API SURFACE – exported names must still exist
# ===========================================================================


class TestPublicAPISurface:
    """Names exported from the top-level package must remain importable."""

    @pytest.mark.parametrize(
        "name",
        [
            "Dataset",
            # Station is not directly exported; access via dataset.stations[0]
            "SensorData",
            "Obstype",
            "ModelObstype",
            "ModelObstype_Vectorfield",
            "Analysis",
            "ModelTimeSeries",
            "GEEStaticDatasetManager",
            "GEEDynamicDatasetManager",
            "WhiteSet",
            "import_dataset_from_pkl",
            "build_template_prompt",
            "connect_to_gee",
            "add_FileHandler",
            "add_StreamHandler",
            "default_GEE_datasets",
            "demo_datafile",
            "demo_metadatafile",
            "demo_template",
            "Settings",
            "__version__",
        ],
    )
    def test_name_exists_in_package(self, name):
        assert hasattr(
            metobs_toolkit, name
        ), f"'{name}' is no longer accessible from metobs_toolkit"

    def test_version_is_string(self):
        assert isinstance(metobs_toolkit.__version__, str)
        # sanity: non-empty version string
        assert metobs_toolkit.__version__.strip() != ""

    def test_demo_paths_exist(self):
        assert Path(metobs_toolkit.demo_datafile).exists()
        assert Path(metobs_toolkit.demo_metadatafile).exists()
        assert Path(metobs_toolkit.demo_template).exists()
