#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Global configuration class for metobs_toolkit.

This module provides a singleton Settings class that stores default configuration
values that can be accessed and modified throughout the package.

@author: thoverga
"""

from __future__ import annotations

import logging
from typing import Dict, Any, Optional, Union
from copy import deepcopy

# import settings modules

from metobs_toolkit.settings_collection.version import __version__
from metobs_toolkit.settings_collection.label_defenitions import (
    label_defs,
    gapfill_label_group,
    failed_gapfill_label_group,
    qc_label_group,
    scatter,
    line,
    vline,
)
from metobs_toolkit.settings_collection.plotting_defaults import default_plot_settings

logger = logging.getLogger("<metobs_toolkit>")


class MetObsSettingsError(Exception):
    """Raised when an error occurs in the MetObsSettings."""

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)


class Settings:
    """
    Singleton configuration class for metobs_toolkit.

    This class stores global settings that are used across multiple functions
    and methods in the package. Settings stored here affect behavior throughout
    the toolkit (e.g., timezone handling, label definitions, logging).

    Note: Settings specific to a single function or method are NOT stored here.
    Those are defined as default argument values in the function/method signature.

    The class uses the singleton pattern to ensure only one instance exists.
    Settings are accessed and modified using class methods, so no instantiation
    is required by the user.

    Examples
    --------
    >>> import metobs_toolkit
    >>> # Get current setting
    >>> metobs_toolkit.Settings.get("store_tz")
    'UTC'

    """

    _instance: Optional["Settings"] = None
    _initialized: bool = False

    # Default configuration values
    _defaults: Dict[str, Any] = {
        "version": __version__,
        # Label defenitions
        "label_def": label_defs,
        "gapfill_label_group": gapfill_label_group,
        "failed_gapfill_label_group": failed_gapfill_label_group,
        "qc_label_group": qc_label_group,
        # Logging defaults
        "log_level": "WARNING",
        "log_format": "LOG:: %(levelname)s - %(message)s",
        # Data storage settings,
        "store_tz": "UTC",
        # Printing
        "print_config": {
            "max_width": 80,
            "item_indent": " " * 2,
            "title_char": "=",
        },
        # Plotting defaults
        "plotting_settings": default_plot_settings,
    }

    _config: Dict[str, Any] = {}

    def __new__(cls) -> "Settings":
        """Ensure only one instance exists (singleton pattern)."""
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self) -> None:
        """Initialize settings with defaults (only once)."""
        if not Settings._initialized:
            Settings._config = deepcopy(Settings._defaults)
            Settings._initialized = True

    def __str__(self) -> str:
        """Return a string representation of the Settings."""
        return f"MetObs Settings object"

    def __repr__(self) -> str:
        """Return a human-readable representation of the Settings."""
        return f"MetObs Settings object"

    @classmethod
    def get(cls, key: str, default: Any = None) -> Any:
        """
        Get a configuration value.

        Parameters
        ----------
        key : str
            The configuration key to retrieve. Use dot notation for nested
            keys (e.g., "label_def.goodrecord.label").
        default : Any, optional
            Default value if key is not found. Default is None.

        Returns
        -------
        Any
            The configuration value.

        Examples
        --------
        >>> Settings.get("store_tz")
        'UTC'
        >>> Settings.get("label_def.goodrecord.label")
        'ok'
        >>> Settings.get("nonexistent_key", "fallback")
        'fallback'
        """
        # Ensure initialized
        cls()

        # Handle dot notation for nested keys
        keys = key.split(".")
        value = cls._config

        try:
            for k in keys:
                value = value[k]
            logger.debug(f"Settings: Retrieved: {key} --> {value}")
            return value
        except KeyError:
            if default is not None:
                logger.warning(
                    f"Settings: {key} not found. Returning default: {default}"
                )
                return default
            raise MetObsSettingsError(
                f"Settings: Key '{key}' not found and no default provided."
            )
        except Exception as e:
            raise MetObsSettingsError(f"Settings: Error retrieving key '{key}': {e}")

    @classmethod
    def set(cls, key: str, value: Any) -> None:
        """
        Set a configuration value.

        Parameters
        ----------
        key : str
            The configuration key to set. Use dot notation for nested
            keys (e.g., "label_def.goodrecord.color").
        value : Any
            The value to set.

        Examples
        --------
        >>> Settings.set("store_tz", "Europe/Brussels")
        >>> Settings.set("log_level", "DEBUG")
        """
        # Ensure initialized
        cls()

        # Handle dot notation for nested keys
        keys = key.split(".")
        config = cls._config

        for k in keys[:-1]:
            if k not in config:
                config[k] = {}
            config = config[k]

        config[keys[-1]] = value
        logger.debug(f"Settings: {key} set to {value}")

    @classmethod
    def reset(cls, key: Optional[str] = None) -> None:
        """
        Reset configuration to defaults.

        Parameters
        ----------
        key : str, optional
            Specific key to reset. If None, reset all settings to defaults.

        Examples
        --------
        >>> Settings.reset()  # Reset all settings
        >>> Settings.reset("plot.dpi")  # Reset only plot.dpi
        """
        if key is None:
            cls._config = deepcopy(cls._defaults)
            logger.debug("Settings: All settings reset to defaults")
        else:
            default_value = cls._get_default(key)
            if default_value is not None:
                cls.set(key, deepcopy(default_value))
                logger.debug(f"Settings: {key} reset to default value: {default_value}")

    @classmethod
    def _get_default(cls, key: str) -> Any:
        """Get default value for a key."""
        keys = key.split(".")
        value = cls._defaults

        try:
            for k in keys:
                value = value[k]
            return value
        except (KeyError, TypeError):
            return None

    @classmethod
    def to_dict(cls) -> Dict[str, Any]:
        """
        Return all settings as a dictionary.

        Returns
        -------
        dict
            Copy of all current settings.

        Examples
        --------
        >>> config = Settings.to_dict()
        >>> config["log_level"]
        'WARNING'
        """
        cls()
        return deepcopy(cls._config)

    @classmethod
    def get_info(cls, printout: bool = True) -> Union[str, None]:
        """
        Get a formatted string with all current settings.

        Parameters
        ----------
        printout : bool, optional
            If True, print the info string. If False, return it. Default is True.

        Returns
        -------
        str or None
            If printout is False, returns the info string. Otherwise, returns None.

        """
        # Import here to avoid circular imports
        from metobs_toolkit.backend_collection import printing_collection as printing

        cls()

        def _format_nested_dict(d: dict, indent_level: int = 1) -> str:
            """Recursively format nested dictionary."""
            result = ""
            for key, val in d.items():
                if isinstance(val, dict):
                    # Print key as a header, then recurse
                    result += printing.print_fmt_line(f"{key}:", identlvl=indent_level)
                    result += _format_nested_dict(val, indent_level + 1)
                else:
                    result += printing.print_fmt_line(
                        f"{key}: {val}", identlvl=indent_level
                    )
            return result

        infostr = ""
        infostr += printing.print_fmt_title("MetObs Toolkit Settings")
        infostr += _format_nested_dict(cls._config, indent_level=1)

        if printout:
            print(infostr)
            return None
        else:
            return infostr

    # ------------------------------------------
    #   Methods Specific to plotting
    # ------------------------------------------
    @classmethod
    def _get_color_from_label(cls, label) -> str:
        cls()
        return {
            group["label"]: group["plotkwargs"]["color"]
            for group in cls.get("label_def").values()
        }.get(label, "")

    @classmethod
    def _label_to_qccheckmap(cls) -> Dict[str, str]:
        cls()
        return {val["label"]: key for key, val in cls.get("label_def").items()}

    @classmethod
    def _flag_plot_as_scatter(cls, label: str) -> bool:
        cls()
        labelmap = cls._label_to_qccheckmap()
        plot_as = cls.get(f"label_def.{labelmap[label]}.plot_as", None)
        if plot_as == scatter:
            return True
        else:
            return False

    @classmethod
    def _flag_plot_as_line(cls, label: str) -> bool:
        cls()
        labelmap = cls._label_to_qccheckmap()
        plot_as = cls.get(f"label_def.{labelmap[label]}.plot_as", None)
        if plot_as == line:
            return True
        else:
            return False

    @classmethod
    def _flag_plot_as_vline(cls, label: str) -> bool:
        cls()
        labelmap = cls._label_to_qccheckmap()
        plot_as = cls.get(f"label_def.{labelmap[label]}.plot_as", None)
        if plot_as == vline:
            return True
        else:
            return False

    # ------------------------------------------
    #    Methods specific to xarray conv
    # ------------------------------------------
    @classmethod
    def _label_to_numericmap(cls) -> Dict[str, int]:
        cls()
        return {
            val["label"]: int(val["numeric_val"])
            for val in cls.get("label_def").values()
        }
