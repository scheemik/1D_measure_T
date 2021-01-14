# Switchboard test script
"""
Description:
This is a script to test the switchboard for the Dedalus experiment
It will relations between different values to insure the assumptions needed
    for the experiment are satisfied
"""


"""
Each test should do the following:
1. Create inputs
2. Execute code being tested and capture the output
3. Compare the output with an expected result
"""

import sys
sys.path.append("../") # Adds higher directory to python modules path
import importlib
switchboard_module = "_code_files.switchboard"
sbp = importlib.import_module(switchboard_module)
# sbp = __import__("_code_files/switchboard.py")

def test_always_passes():
    assert True

def test_always_fails():
    assert False

def test_uppercase():
    assert "loud noises".upper() == "LOUD NOISES"

def test_params():
    assert sbp.nz == 1024
