#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test script to verify GEE authentication works in the current environment.

This script helps debug GEE authentication issues by:
1. Checking for credentials files
2. Testing ee.Initialize()
3. Making a simple GEE API call

Usage:
    python test_gee_auth.py
"""

import os
import sys
from pathlib import Path
import ee

# Add source to path
repodir = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(repodir / "src"))


def check_credentials_files():
    """Check if GEE credentials files exist."""
    print("=" * 70)
    print("Checking for GEE credentials files...")
    print("=" * 70)

    home = Path.home()
    locations = [
        home / ".config" / "earthengine" / "credentials",
        home / ".config" / "earthengine" / "credentials.json",
        home / ".config" / "EarthEngine" / "credentials",
    ]

    found = False
    for loc in locations:
        if loc.exists():
            print(f"✓ Found credentials: {loc}")
            print(f"  Size: {loc.stat().st_size} bytes")
            print(f"  Permissions: {oct(loc.stat().st_mode)[-3:]}")
            return True
        else:
            print(f"✗ Not found: {loc}")

    if not found:
        print("\n⚠ No GEE credentials found!")
        print("\nTo set up credentials, run:")
        print("  python -c 'import metobs_toolkit; metobs_toolkit.connect_to_gee()'")
        return False

    print()
    return True


def check_environment_variables():
    """Check relevant environment variables."""
    print("=" * 70)
    print("Checking environment variables...")
    print("=" * 70)

    env_vars = [
        "EARTHENGINE_CREDENTIALS_PATH",
        "GOOGLE_APPLICATION_CREDENTIALS",
        "EARTHENGINE_SERVICE_ACCOUNT_CREDENTIALS",
    ]

    for var in env_vars:
        value = os.getenv(var)
        if value:
            print(f"✓ {var} = {value}")
        else:
            print(f"✗ {var} not set")

    print()


def test_ee_direct():
    """Test ee module directly."""
    print("=" * 70)
    print("Testing ee module directly...")
    print("=" * 70)

    try:

        print("✓ ee module imported successfully")

        # Try to initialize
        try:
            ee.Initialize()
            print("✓ ee.Initialize() succeeded")

            # Try a simple API call
            try:
                image = ee.Image("USGS/SRTMGL1_003")
                info = image.getInfo()
                print("✓ Simple GEE API call succeeded")
                print(f"  Test image bands: {info.get('bands', 'N/A')}")
                return True
            except Exception as e:
                print(f"✗ GEE API call failed: {e}")
                return False

        except Exception as e:
            print(f"✗ ee.Initialize() failed: {e}")
            return False

    except ImportError as e:
        print(f"✗ Failed to import ee: {e}")
        return False


def test_metobs_toolkit():
    """Test metobs_toolkit GEE connection."""
    print("=" * 70)
    print("Testing metobs_toolkit.connect_to_gee()...")
    print("=" * 70)

    try:
        import metobs_toolkit

        print("✓ metobs_toolkit imported successfully")

        try:
            metobs_toolkit.connect_to_gee()
            print("✓ metobs_toolkit.connect_to_gee() succeeded")
            return True
        except Exception as e:
            print(f"✗ metobs_toolkit.connect_to_gee() failed: {e}")
            return False

    except ImportError as e:
        print(f"✗ Failed to import metobs_toolkit: {e}")
        return False


def main():
    """Run all tests."""
    print("\n")
    print("*" * 70)
    print("*" + " " * 68 + "*")
    print("*" + "  GEE Authentication Test".center(68) + "*")
    print("*" + " " * 68 + "*")
    print("*" * 70)
    print("\n")

    # Run checks
    creds_ok = check_credentials_files()
    check_environment_variables()
    ee_ok = test_ee_direct()
    print()
    toolkit_ok = test_metobs_toolkit()

    # Summary
    print()
    print("=" * 70)
    print("Summary")
    print("=" * 70)
    print(f"Credentials found: {'✓ Yes' if creds_ok else '✗ No'}")
    print(f"ee module works: {'✓ Yes' if ee_ok else '✗ No'}")
    print(f"metobs_toolkit works: {'✓ Yes' if toolkit_ok else '✗ No'}")
    print("=" * 70)

    if creds_ok and ee_ok and toolkit_ok:
        print("\n✓ All tests passed! GEE authentication is working correctly.")
        return 0
    else:
        print("\n✗ Some tests failed. See output above for details.")
        return 1


if __name__ == "__main__":
    sys.exit(main())