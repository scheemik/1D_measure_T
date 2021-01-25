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
# example tests
@pytest.mark.myfirstmark
def test_numbers():
    assert 2.0 == 2.0

@pytest.mark.parametrize("multiples_of_2", [
    2, 4, 8, 16, 32,
])
def test_is_divisible_by_2(multiples_of_2):
    assert multiples_of_2%2 == 0, "This should be divisible by 2"

###############################################################################
# Set up fixtures
@pytest.fixture(params=[sbp.N_0])
def arr_N_0(request):
    return request.param

kL_s = np.linspace(0.0, 5.0, 16)
@pytest.fixture(params=kL_s)
def arr_kL(request):
    return request.param

@pytest.fixture(params=[sbp.theta, 0.3])
def arr_theta(request):
    return request.param

@pytest.fixture(params=[sbp.lam_z, 0.5])
def arr_lam_z(request):
    return request.param

@pytest.fixture(params=[0, 1, 2, 3, 4, 5])
def arr_n_layers(request):
    return request.param

@pytest.fixture(params=[sbp.R_i, 0, 0.2, 0.5, 0.7])
def arr_R_i(request):
    return request.param

@pytest.fixture(params=[sbp.nz])
def arr_nz(request):
    return request.param

@pytest.fixture(params=[sbp.p_n_steps])
def arr_p_n_steps(request):
    return request.param

@pytest.fixture(params=[sbp.p_o_steps])
def arr_p_o_steps(request):
    return request.param

@pytest.fixture(params=[sbp.T_cutoff])
def arr_T_cutoff(request):
    return request.param

@pytest.fixture(params=[sbp.T_keep])
def arr_T_keep(request):
    return request.param

###############################################################################
# Physical parameter tests

def test_dispersion_relation(arr_N_0, arr_theta, arr_lam_z):
    """
    Checking to make sure the dispersion relation is satifsied within a certain tolerance
    Using eq (3.54) from Sutherland 2010
    """
    N_0 = arr_N_0
    omega, m, k, k_total, lam_x = sbp.calc_wave_params(arr_N_0, arr_theta, arr_lam_z)
    # Set up both sides of the equation
    LHS = omega**2
    RHS = N_0**2 * k**2 / (k**2 + m**2)
    tolerance = 1E-14       # usually fine down to 10^-15 or 10^-16
    assert abs(LHS - RHS) < tolerance, "Dispersion relation not satisfied"

def test_propagation_angle(arr_N_0, arr_theta, arr_lam_z):
    """
    Checking to make sure theta matches with omega and N_0 within a certain tolerance
    Using eq (3.57) from Sutherland 2010
    """
    N_0 = arr_N_0
    theta = arr_theta
    omega, m, k, k_total, lam_x = sbp.calc_wave_params(arr_N_0, arr_theta, arr_lam_z)
    # Set up both sides of the equation
    LHS = omega
    RHS = N_0 * np.cos(theta)
    tolerance = 1E-14       # usually fine down to 10^-15 or 10^-16
    assert abs(LHS - RHS) < tolerance, "Propagation angle doesn't match"

###############################################################################
# Background profile tests

def test_0_layers(arr_R_i):
    """
    Check to make sure the interface ratio is set appropriately when the
        number of layers is zero
    """
    R_i = sbp.calc_R_i(arr_R_i, 0)
    assert R_i == 0, "Interface ratio for zero layers set incorrectly"

def test_top_window():
    """
    Check to make sure the parameters for the gaussian forcing windows
        are calculated correctly for the top window
    """
    z0_dis = 0.0
    lam_z  = 1.0
    a, b   = 3, 1
    a_bf, b_bf, c_bf = sbp.calc_bf_win_params(z0_dis, 1, lam_z, a, b)
    assert a_bf ==  3.0
    assert b_bf ==  1.0
    assert c_bf >   z0_dis
    assert c_bf ==  1.5

def test_bottom_window():
    """
    Check to make sure the parameters for the gaussian forcing windows
        are calculated correctly for the bottom window
    """
    zf_dis = -5.0
    lam_z  = 1.0
    a, b   = 3, 1
    a_bf, b_bf, c_bf = sbp.calc_bf_win_params(zf_dis, -1, lam_z, a, b)
    assert a_bf ==  3.0
    assert b_bf ==  1.0
    assert c_bf <   zf_dis
    assert c_bf == -6.5

###############################################################################
# Domain parameter tests

def test_Lz():
    """
    Check to make sure Lz is an integer number of vertical wavelengths
        to insure periodic boundary conditions are applicable
    """
    n_lambda = sbp.Lz / sbp.lam_z
    assert n_lambda - int(n_lambda) == 0, "Simulation domain is not an integer number of lambda_z"

def test_dis_domain(arr_N_0, arr_theta, arr_lam_z, arr_kL, arr_n_layers, arr_R_i):
    """
    Check to make sure the display domain is an integer number of vertical wavelengths
    """
    omega, m, k, k_total, lam_x = sbp.calc_wave_params(arr_N_0, arr_theta, arr_lam_z)
    L = sbp.calc_layer_thickness(arr_kL, k)
    z_I, z0_str, zf_str, z_T, zf_dis, Lz_dis = sbp.calc_structure_depths(sbp.z0_dis, arr_lam_z, L, arr_n_layers, arr_R_i)
    n_lambda = Lz_dis / arr_lam_z
    assert n_lambda - int(n_lambda) == 0, "Display domain is not an integer number of lambda_z"

###############################################################################
# Computational parameter tests

def test_nz(arr_nz):
    """
    Check to make sure the number of points in the z display domain is both
        an integer and a power of 2 to allow fast Fourier transforms
    """
    nz = arr_nz
    assert isinstance(nz, int), "nz should be an integer"
    # Use bit operation to check for a power of 2
    assert nz != 0 and (nz & (nz-1) == 0), "nz should be a power of 2"

def test_nz_dealias(arr_nz, arr_kL, arr_n_layers, arr_R_i):
    """
    Check to make sure the total number of points in the z direction
        multiplied by the dealias factor is an integer
    """
    L = sbp.calc_layer_thickness(arr_kL, sbp.k)
    z_I, z0_str, zf_str, z_T, zf_dis, Lz_dis = sbp.calc_structure_depths(sbp.z0_dis, sbp.lam_z, L, arr_n_layers, arr_R_i)
    z0, zf, Lz = sbp.calc_sim_domain(sbp.z0_dis, zf_dis, sbp.a_bf, sbp.a_sp)
    nz_sim, dz = sbp.calc_nz_sim(arr_nz, Lz, sbp.Lz_dis)
    nz_da  = nz_sim * sbp.dealias
    assert nz_da - int(nz_da) == 0, "nz_sim*dealias should be an integer"

def test_n_steps_per_oscillation():
    """
    Check to make sure the number of total timesteps in the simulation is
        larger than the number of timesteps per oscillation period
    If this isn't true, n_T = 2^-n, so it will be less than 1
    """
    n_T = sbp.calc_oscillation_periods(sbp.p_n_steps, sbp.p_o_steps)
    assert n_T > 1, "Total timesteps must be greater than timesteps per oscillation"
