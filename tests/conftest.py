"""Test configuration for the SBTN_LEAF package."""

from pathlib import Path
import atexit
import os
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_PATH = PROJECT_ROOT / "src"

if str(SRC_PATH) not in sys.path:
    sys.path.insert(0, str(SRC_PATH))


_PREVIOUS_CWD = Path.cwd()
os.chdir(SRC_PATH)


@atexit.register
def _restore_working_directory():
    os.chdir(_PREVIOUS_CWD)
