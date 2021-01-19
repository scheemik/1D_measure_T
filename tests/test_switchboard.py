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

import numpy as np

# Funny workaround to get a relative module import for the switchboard
import os, sys
currentdir = os.path.dirname(os.path.realpath(__file__))
parentdir = os.path.dirname(currentdir)
sys.path.append(parentdir)
from _code_files import switchboard as sbp

###############################################################################
import pytest

###############################################################################
# Set up fixtures
@pytest.fixture
def example_array():
    return [1, 2, 3]

###############################################################################
# example tests
def test_array_sum(example_array):
    assert sum(example_array) == 6, "Should be 6"

@pytest.mark.myfirstmark
def test_numbers():
    assert 2.0 == 2.0

@pytest.mark.parametrize("multiples_of_2", [
    2, 4, 8, 16, 32,
])
def test_is_divisible_by_2(multiples_of_2):
    assert multiples_of_2%2 == 0, "This should be divisible by 2"

###############################################################################
# Physical parameter tests

def test_dispersion_relation():
    """
    Checking to make sure the dispersion relation is satifsied within a certain tolerance
    Using eq (3.54) from Sutherland 2010
    """
    LHS = sbp.omega**2
    RHS = sbp.N_0**2 * sbp.k**2 / (sbp.k**2 + sbp.m**2)
    tolerance = 1E-14       # usually fine down to 10^-15 or 10^-16
    assert abs(LHS - RHS) < tolerance, "Dispersion relation not satisfied"

def test_propagation_angle():
    """
    Checking to make sure theta matches with omega and N_0 within a certain tolerance
    Using eq (3.57) from Sutherland 2010
    """
    LHS = sbp.omega
    RHS = sbp.N_0 * np.cos(sbp.theta)
    tolerance = 1E-14       # usually fine down to 10^-15 or 10^-16
    assert abs(LHS - RHS) < tolerance, "Propagation angle doesn't match"

###############################################################################
# Domain parameter tests

def test_Lz():
    """
    Check to make sure Lz is an integer number of wavelengths
        to insure periodic boundary conditions are applicable
    """
    n_lambda = sbp.Lz / sbp.lam_z
    assert n_lambda - int(n_lambda) == 0, "Z domain is not an integer number of lambda"

###############################################################################
# Computational parameter tests

def test_nz():
    """
    Checking to make sure the number of points in the z direction is both
        an integer and a power of 2 to allow fast Fourier transforms
    """
    nz = sbp.nz
    assert isinstance(nz, int), "nz should be an integer"
    # Use bit operation to check for a power of 2
    assert nz != 0 and (nz & (nz-1) == 0), "nz should be a power of 2"
