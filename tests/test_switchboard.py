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

import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from _code_files import switchboard as sbp
# sys.path.append("../") # Adds higher directory to python modules path
# import importlib
# switchboard_module = "_code_files.switchboard"
# sbp = importlib.import_module(switchboard_module)
# sbp = __import__("_code_files/switchboard.py")

def test_always_passes():
    assert True

# def test_always_fails():
#     assert False

def test_params():
    assert sbp.mL == 2.0

def test_switchboard():
    assert sbp.nz == 1024
