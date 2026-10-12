"""Core WEC configuration works without plotting or workbook packages."""

import subprocess
import sys


def test_core_import_and_array_mcr_without_optional_packages():
    # A separate interpreter catches imports triggered by the public package,
    # even when the developer's environment happens to contain the extras.
    script = """
import sys

class BlockOptional:
    def find_spec(self, fullname, path=None, target=None):
        if fullname.split('.')[0] in {'matplotlib', 'openpyxl', 'control'}:
            raise ModuleNotFoundError(fullname)

sys.meta_path.insert(0, BlockOptional())
from wecsim import WEC, mcr_wave_statistics
assert WEC('base install').name == 'base install'
cases = mcr_wave_statistics([[0, 8], [2.5, 1]], [1_200_000])
assert len(cases) == 1 and cases[0].wave_height == 2.5
"""
    subprocess.run([sys.executable, "-c", script], check=True,
                   capture_output=True, text=True)
