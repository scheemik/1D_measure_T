"""
Runs an idealized steady state test of various functions.
Run with $ python3 steady_state_test.py

Overview of post-processing operations:
* Import simulation data
    * Plot full wavefield
* Trim data in space (restrict depths)
* Trim data in time (remove transients at beginning)
    * Plot trimmed wavefield
* Perform complex demodulation
    * Plot data in spectral form
* Isolate both directions of the wave
    * Plot down and upward waves (real parts)
* Find amplitude by multiplying by complex conjugate
    * Plot amplitude vs depth
* Calculate incident and transmitted wave amplitudes
* Calculate transmission coefficient
    * Compare to analytical value

"""

import h5py
import numpy as np
import matplotlib
#matplotlib.use('Agg')
from matplotlib import ticker
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from dedalus.extras.plot_tools import quad_mesh

run_name = "steady_state_test"

###############################################################################
# Helper functions
#   This import assumes the helper functions are in the same directory as the core code
import helper_functions as hf
import helper_functions_CD as hfCD

###############################################################################
# Import SwitchBoard Parameters (sbp)
#   This import assumes the switchboard is in the same directory as the core code
import switchboard as sbp
# Physical parameters
nu          = sbp.nu            # [m^2/s]       Viscosity (momentum diffusivity)
f_0         = sbp.f_0           # [s^-1]        Reference Coriolis parameter
g           = sbp.g             # [m/s^2]       Acceleration due to gravity
# Problem parameters
mL          = sbp.mL            # []            Ratio of vertical wavelength to stratification length
N_0         = sbp.N_0           # [rad/s]       Reference stratification
lam_x       = sbp.lam_x         # [m]           Horizontal wavelength
lam_z       = sbp.lam_z         # [m]           Vertical wavelength
k           = sbp.k             # [m^-1]        Horizontal wavenumber
m           = sbp.m             # [m^-1]        Vertical wavenumber
k_total     = sbp.k_total       # [m^-1]        Total wavenumber
theta       = sbp.theta         # [rad]         Propagation angle from vertical
omega       = sbp.omega         # [rad s^-1]    Wave frequency
T           = sbp.T             # [s]           Wave period

# Parameters
nz          = sbp.nz
nz_sim      = sbp.nz_sim
z0_dis      = sbp.z0_dis
zf_dis      = sbp.zf_dis
lz_lambda = abs(z0_dis-zf_dis)/lam_z
print("Lz_dis/lambda=",lz_lambda)
z0          = sbp.z0
zf          = sbp.zf

# Steps in time and space
dt          = sbp.dt
dz          = sbp.dz#abs(z0_dis - zf_dis)/nz

# Mesuring depths for incident and transient waves
z_I         = sbp.z_I
z_T         = sbp.z_T

plt_fd      = sbp.plot_full_domain
T_cutoff    = sbp.T_cutoff

###############################################################################
# Create ideal, steady-state, periodic boundary simulation data

# Set amplitude coefficients for up and down wave components
A           = 1.0
B           = 1.0

# Find space and time axes (z, t)
#z = np.arange(zf, z0, dz) # careful about the order of endpoints
z = sbp.z
# z and psi arrays come out sorted from most positive to most negative on z axis
#   This flips things around (ud = up / down)
z = np.flip(z)
# psi = np.flipud(psi)

stop_sim_time = sbp.sim_time_stop
nt = stop_sim_time / dt
t = np.linspace(0.0, stop_sim_time, int(nt))
# Find wavenumbers
kz_dis = np.fft.fftfreq(nz, dz)
kz_sim = np.fft.fftfreq(nz_sim, dz)
# Make a space time meshgrid
tm, zm = np.meshgrid(t, z)

# Form of wave field following Mercier et al 2008 eq 5
#   Up and down terms separated for ease of access later
up  = A*np.cos(omega*tm - m*zm)
dn  = B*np.cos(omega*tm + m*zm)
# Program psi field directly
psi = up + dn

# Form complex-valued version of wave field following Mercier et al 2008 eq 7
#   Equivalent to taking FT in time of psi, masking the negative frequencies,
#    and multiplying the rest by 2 to maintain the same amplitude, then iFT
#   The real part of these should equal the corresponding fields above
#   Again, up and down separated for ease of access later
up_c = A*np.exp(1j*(omega*tm - m*zm))
dn_c = B*np.exp(1j*(omega*tm + m*zm))
psi_hat = up_c + dn_c

###############################################################################
# Creating stratification background profile
BP_array = hf.BP_n_layers(0, z, sbp.z0_str, sbp.zf_str)

###############################################################################
# Plot vertical profiles and full wavefield

if sbp.plot_windows:
    hf.plot_v_profiles(z, BP_array, sbp.win_bf_array, sbp.win_sp_array, mL=mL, theta=theta, omega=omega, z_I=z_I, z_T=z_T, z0_dis=z0_dis, zf_dis=zf_dis, title_str=run_name, filename='ss_1D_windows.png')

if sbp.plot_spacetime:
    hf.plot_z_vs_t(z, t, T, psi, BP_array, mL, theta, omega, z_I=z_I, z_T=z_T, z0_dis=z0_dis, zf_dis=zf_dis, plot_full_domain=True, T_cutoff=T_cutoff, title_str=run_name, filename='ss_1D_wave.png')


###############################################################################
# Trim data

# Trim in space
z_tr, psi_tr = hf.trim_data_z(z, psi, z0_dis, zf_dis)
BP_tr        = hf.BP_n_layers(0, z_tr, sbp.z0_str, sbp.zf_str)

# Trim in time
t_tr, psi_tr = hf.trim_data_t(t, psi_tr, T_cutoff, T)

# Plot trimmed wavefield
if sbp.plot_spacetime:
    hf.plot_z_vs_t(z_tr, t_tr, T, psi_tr, BP_tr, mL, theta, omega, z_I=z_I, z_T=z_T, z0_dis=z0_dis, zf_dis=zf_dis, plot_full_domain=plt_fd, T_cutoff=T_cutoff, title_str=run_name, filename='ss_1D_wave_tr.png')

###############################################################################
# Plot spectral form of data

freqs, ks, spec_data = hfCD.z_t_to_k_omega(t_tr, z_tr, psi_tr, dt, dz)

if sbp.plot_spectra:
    hf.plot_spectral(ks, freqs, spec_data.real, spec_data.imag, mL, theta, omega, c_map='viridis', title_str=run_name, filename='ss_1D_spectra.png')

# raise SystemExit(0)
###############################################################################
# Perform complex demodulation

t_then_z = True
up_field, dn_field = hfCD.Complex_Demodulate(t_then_z, t, z, kz_sim, psi, dt, omega)
tr_up_field, tr_dn_field = hfCD.Complex_Demodulate(t_then_z, t_tr, z_tr, kz_dis, psi_tr, dt, omega)

###############################################################################
# Plotting up and downward propagating waves

# Trimmed case
if sbp.plot_up_dn:
    hf.plot_z_vs_t(z_tr, t_tr, T, tr_up_field.real, BP_tr, mL, theta, omega, z_I=z_I, z_T=z_T, z0_dis=z0_dis, zf_dis=zf_dis,  plot_full_domain=plt_fd, T_cutoff=None, title_str=run_name+' up', filename='ss_1D_up_field_tr.png')
    hf.plot_z_vs_t(z_tr, t_tr, T, tr_dn_field.real, BP_tr, mL, theta, omega, z_I=z_I, z_T=z_T, z0_dis=z0_dis, zf_dis=zf_dis, plot_full_domain=plt_fd, T_cutoff=None, title_str=run_name+' dn', filename='ss_1D_dn_field_tr.png')

# Full domain
if sbp.plot_up_dn:
    hf.plot_z_vs_t(z, t, T, up_field.real, BP_array, mL, theta, omega, z_I=z_I, z_T=z_T, z0_dis=z0_dis, zf_dis=zf_dis,  plot_full_domain=True, T_cutoff=None, title_str=run_name+' up', filename='ss_1D_up_field.png')
    hf.plot_z_vs_t(z, t, T, dn_field.real, BP_array, mL, theta, omega, z_I=z_I, z_T=z_T, z0_dis=z0_dis, zf_dis=zf_dis, plot_full_domain=True, T_cutoff=None, title_str=run_name+' dn', filename='ss_1D_dn_field.png')

###############################################################################
# Multiply the downward wavefield by it's complex-conjugate to get AA^*
#   Plot this amplitude across time at two specific depths:
if sbp.plot_amplitude:
    hf.plot_A_of_I_T(z_tr, t_tr, T, tr_dn_field, mL, theta, omega, z_I, z_T, plot_full_domain=plt_fd, T_cutoff=T_cutoff, title_str=run_name, filename='ss_1D_A_of_I_T.png')
# Plot this amplitude (averaged across time) as a function of depth
if sbp.plot_amplitude:
    hf.plot_AA_for_z(z_tr, BP_tr, hf.AAcc(tr_dn_field), mL, theta, omega, z_I=z_I, z_T=z_T, z0_dis=z0_dis, zf_dis=zf_dis, plot_full_domain=plt_fd, title_str=run_name, filename='ss_1D_AA_for_z.png')

###############################################################################
# Measuring the transmission coefficient

I_, T_, AAcc_I, AAcc_T = hf.measure_T(tr_dn_field, z_tr, z_I, z_T, T_skip=None, T=T, t=t_tr)
big_T = T_/I_
print("Transmission coefficient is:", big_T)

# raise SystemExit(0)
###############################################################################
# Profiling the code
# profile_it = False
# if profile_it == True:
#     import cProfile
#     cProfile.run('Complex_Demodulate(t_then_z, t, z, kz, psi, dt)', 'restats')
#     import pstats
#     from pstats import SortKey
#     p = pstats.Stats('restats')
#     p.sort_stats(SortKey.CUMULATIVE).print_stats(10)

###############################################################################
# Extra plots

plot_CD_checks = False
if plot_CD_checks:
    # Add up and down fields to see if they reproduce the original psi field
    up_plus_dn = plot_up.real + plot_dn.real
    hf.plot_z_vs_t(z, t, T, up_plus_dn, plot_BP_, mL, theta, omega, z0_dis=z0_dis, zf_dis=zf_dis, plot_full_domain=plt_fd, nT=None, title_str=run_name, filename='ss_1D_up_plus_dn.png')
    # Plot the difference, which ideally should be zero everywhere
    CD_diff = psi.real - up_plus_dn
    hf.plot_z_vs_t(z, t, T, CD_diff, plot_BP_, mL, theta, omega, z0_dis=z0_dis, zf_dis=zf_dis, plot_full_domain=plt_fd, nT=None, title_str=run_name, filename='ss_1D_CD_diff.png')
