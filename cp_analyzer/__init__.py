"""
CP-Analyzer: A Constant-Potential Platform for Electrocatalytic Free Energy Calculations
"""

import shutil
from pathlib import Path
import warnings
import sys
import os

__version__ = "2.0.0"


def initialize_config() -> None:
    """Copy default config to current working directory if not exists"""
    cwd = Path.cwd()
    package_dir = Path(__file__).parent

    default_config = package_dir / "config_default.py"
    target_config = cwd / "config.py"

    if not default_config.exists():
        warnings.warn("config_default.py not found in package.")
        return

    if target_config.exists():
        print(f"✅ Using existing config.py in current directory.")
        return

    try:
        shutil.copy(default_config, target_config)
        print(f"✅ Default config.py has been created in current directory: {target_config}")
        print("   Please edit config.py before running analysis.")
    except Exception as e:
        print(f"⚠️ Failed to create config.py: {e}")


def get_config_path():
    """Return config.py path in current working directory"""
    return Path.cwd() / "config.py"