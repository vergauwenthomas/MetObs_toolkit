import os
import sys
from pathlib import Path


lib_folder = Path(__file__).resolve().parents[2]
sys.path.insert(0,str(lib_folder))

import metobs_toolkit