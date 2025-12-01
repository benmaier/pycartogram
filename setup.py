# Legacy setup.py for backwards compatibility.
# All configuration is now in pyproject.toml.
# This file enables `pip install -e .` for older pip versions.

from setuptools import setup

if __name__ == "__main__":
    setup()
