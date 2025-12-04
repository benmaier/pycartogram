#!/usr/bin/env python
"""Test if cCartogram imports correctly after LC_RPATH fix."""

try:
    import cCartogram
    print("cCartogram imported successfully!")
    print(f"Module location: {cCartogram.__file__}")
except ImportError as e:
    print(f"Import failed: {e}")
