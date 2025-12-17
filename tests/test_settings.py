"""Tests for the Settings class."""

import pytest

import metobs_toolkit
from metobs_toolkit.settings_collection.settings import Settings, MetObsSettingsError


@pytest.fixture(autouse=True)
def reset_settings():
    """Reset settings to defaults before and after each test."""
    Settings.reset()
    yield
    Settings.reset()


class TestSettingsGet:
    """Tests for Settings.get() method."""

    def test_get_top_level_setting(self):
        """Test getting a top-level setting."""
        store_tz = Settings.get("store_tz")
        assert store_tz == "UTC"

    def test_get_nested_setting_with_dot_notation(self):
        """Test getting a nested setting using dot notation."""
        label = Settings.get("label_def.goodrecord.label")
        assert label == "ok"

    def test_get_nonexistent_key_with_default(self):
        """Test getting a nonexistent key returns the default value."""
        result = Settings.get("nonexistent_key", "fallback")
        assert result == "fallback"

    def test_get_nonexistent_key_without_default_raises(self):
        """Test getting a nonexistent key without default raises error."""
        with pytest.raises(MetObsSettingsError):
            Settings.get("nonexistent_key")

    def test_get_version(self):
        """Test getting the version setting."""
        version = Settings.get("version")
        assert version is not None
        assert isinstance(version, str)

    def test_get_log_level(self):
        """Test getting the log level setting."""
        log_level = Settings.get("log_level")
        assert log_level == "WARNING"


class TestSettingsSet:
    """Tests for Settings.set() method."""

    def test_set_top_level_setting(self):
        """Test setting a top-level setting."""
        Settings.set("store_tz", "Europe/Brussels")
        assert Settings.get("store_tz") == "Europe/Brussels"

    def test_set_nested_setting_with_dot_notation(self):
        """Test setting a nested setting using dot notation."""
        Settings.set("label_def.goodrecord.color", "red")
        assert Settings.get("label_def.goodrecord.color") == "red"

    def test_set_new_key(self):
        """Test setting a new key that doesn't exist."""
        Settings.set("new_setting", "new_value")
        assert Settings.get("new_setting") == "new_value"

    def test_set_new_nested_key(self):
        """Test setting a new nested key that doesn't exist."""
        Settings.set("new.nested.setting", "nested_value")
        assert Settings.get("new.nested.setting") == "nested_value"


class TestSettingsReset:
    """Tests for Settings.reset() method."""

    def test_reset_all_settings(self):
        """Test resetting all settings to defaults."""
        # Modify some settings
        Settings.set("store_tz", "Europe/Brussels")
        Settings.set("log_level", "DEBUG")

        # Reset all
        Settings.reset()

        # Check they're back to defaults
        assert Settings.get("store_tz") == "UTC"
        assert Settings.get("log_level") == "WARNING"

    def test_reset_specific_setting(self):
        """Test resetting a specific setting to its default."""
        # Modify a setting
        Settings.set("store_tz", "Europe/Brussels")
        Settings.set("log_level", "DEBUG")

        # Reset only store_tz
        Settings.reset("store_tz")

        # Check store_tz is reset but log_level is not
        assert Settings.get("store_tz") == "UTC"
        assert Settings.get("log_level") == "DEBUG"

    def test_reset_nested_setting(self):
        """Test resetting a nested setting to its default."""
        original_label = Settings.get("label_def.goodrecord.label")

        # Modify nested setting
        Settings.set("label_def.goodrecord.label", "modified")
        assert Settings.get("label_def.goodrecord.label") == "modified"

        # Reset it
        Settings.reset("label_def.goodrecord.label")
        assert Settings.get("label_def.goodrecord.label") == original_label


class TestSettingsToDict:
    """Tests for Settings.to_dict() method."""

    def test_to_dict_returns_dict(self):
        """Test that to_dict returns a dictionary."""
        config = Settings.to_dict()
        assert isinstance(config, dict)

    def test_to_dict_returns_copy(self):
        """Test that to_dict returns a copy, not the original."""
        config = Settings.to_dict()
        config["store_tz"] = "Modified"

        # Original should be unchanged
        assert Settings.get("store_tz") == "UTC"


class TestSettingsGetInfo:
    """Tests for Settings.get_info() method."""

    def test_get_info_printout_returns_none(self):
        """Test that get_info with printout=True returns None."""
        result = Settings.get_info(printout=True)
        assert result is None

    def test_get_info_no_printout_returns_string(self):
        """Test that get_info with printout=False returns a string."""
        result = Settings.get_info(printout=False)
        assert isinstance(result, str)
        assert len(result) > 0

    def test_get_info_contains_settings(self):
        """Test that get_info output contains setting information."""
        result = Settings.get_info(printout=False)
        assert "MetObs Toolkit Settings" in result
        assert "store_tz" in result
        assert "log_level" in result


class TestSettingsSingleton:
    """Tests for Settings singleton pattern."""

    def test_singleton_same_instance(self):
        """Test that Settings always returns the same instance."""
        instance1 = Settings()
        instance2 = Settings()
        assert instance1 is instance2

    def test_settings_accessible_via_metobs_toolkit(self):
        """Test that Settings is accessible via metobs_toolkit module."""
        assert metobs_toolkit.Settings is Settings

class TestSettingsStrRepr:
    """Tests for Settings __str__ and __repr__ methods."""

    def test_str_representation(self):
        """Test string representation of Settings."""
        settings = Settings()
        assert str(settings) == "MetObs Settings object"

    def test_repr_representation(self):
        """Test repr representation of Settings."""
        settings = Settings()
        assert repr(settings) == "MetObs Settings object"
