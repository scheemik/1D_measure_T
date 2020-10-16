# Switchboard for Dedalus experiment
"""
Description:
This is the switchboard for the Dedalus experiment. This file contains parameter values called by multiple scripts. This insures they all get the same values
"""

import numpy as np
from dedalus import public as de
# import sys
# sys.path.append("../") # Adds higher directory to python modules path
# import helper_functions as hf

###############################################################################
# Main parameters, the ones I'll change a lot. Many more below

# Relevant parameters
nz              = 1024          # [] number of grid points in the z direction
mL      = 1                     # [] vertical wave number times step length
theta   = None                  # [] angle between wave's propagation and horizontal (or vertical?)

# Time parameters
p_n_steps   = 10                # [] power of the number of timesteps for the simulation
p_o_steps   = 6                 # [] power of the number of timesteps per oscillation period
T_cutoff    = 3                 # [] number of oscillations to cut from beginning to leave just steady state
#
n_steps     = int(2**p_n_steps) # [] number of timesteps for the simulation
o_steps     = int(2**p_o_steps) # [] number of timesteps per oscillation period
p_n_T       =p_n_steps-p_o_steps# [] power of the number of oscillation periods
if p_n_T <= 0:
    print("Total timesteps must be greater than timesteps per oscillation")
    raise SystemExit(0)
n_T         = int(2**p_n_T)     # [] number of oscillation periods

# Run parameters (to be deleted)
stop_n_periods  = 28            # [] oscillation periods (28, 57)
extend_to_pwr_2 = True          # [] Extend simulation time so nt is a power of 2

# Displayed domain parameters
z0_dis = 0.0                    # [m] the top of the displayed z domain

# Problem parameters
A       = 2.0e-4                # []            Amplitude of boundary forcing
N_0     = 1.0                   # [rad/s]       Reference stratification
f_0     = 0.000                 # [s^-1]        Reference Coriolis parameter
set_case= 1                     # Picks combination of variables to set in switch below
if set_case == 1:
    lam_z   = 1.0 / 4.0             # [m]           Vertical wavelength
    lam_x   = lam_z                 # [m]           Horizontal wavelength
    #
    m       = 2*np.pi / lam_z       # [m^-1]        Vertical wavenumber
    k       = 2*np.pi / lam_x       # [m^-1]        Horizontal wavenumber
    k_total = np.sqrt(k**2 + m**2)  # [m^-1]        Total wavenumber
    theta   = np.arctan(m/k)        # [rad]         Propagation angle from vertical
    omega   = N_0 * np.cos(theta)   # [rad s^-1]    Wave frequency
elif set_case == 2:
    lam_z   = Lz_dis / 8.0          # [m]           Vertical wavelength
    omega   = 0.7071                # [rad s^-1]    Wave frequency
    #
    m       = 2*np.pi / lam_z       # [m^-1]        Vertical wavenumber
    theta   = np.arccos(omega/N_0)  # [rad]         Propagation angle from vertical
    k       = m/np.tan(theta)       # [m^-1]        Horizontal wavenumber
    k_total = np.sqrt(k**2 + m**2)  # [m^-1]        Total wavenumber
    lam_x   = 2*np.pi / k           # [m]           Horizontal wavelength
elif set_case == 3:
    lam_z   = Lz_dis / 8.0          # [m]           Vertical wavelength
    theta   = 0.7854 # 45deg        # [rad]         Propagation angle from vertical
    #
    m       = 2*np.pi / lam_z       # [m^-1]        Vertical wavenumber
    k       = m/np.tan(theta)       # [m^-1]        Horizontal wavenumber
    k_total = np.sqrt(k**2 + m**2)  # [m^-1]        Total wavenumber
    lam_x   = 2*np.pi / k           # [m]           Horizontal wavelength
    omega   = N_0 * np.cos(theta)   # [rad s^-1]    Wave frequency

T       = 2*np.pi / omega       # [s]           Wave period

###############################################################################
# ON / OFF Switches

# Determine whether adaptive time stepping is on or off
adapt_dt                = False
if extend_to_pwr_2 == True:
    adapt_dt = False

# Terms in equations of motion
viscous_term            = True
pressure_term           = True
advection_term          = False
buoyancy_term           = True
diffusivity_term        = True
rotation_term           = False
# Diffusion / dissipation of reflections
use_sponge              = True
use_rayleigh_friction   = False
boundary_forcing_region = True  # If False, waves will be forced over entire domain

# Plotting parameters
plot_spacetime = True
plot_wavespace = True
plot_freqspace = False
plot_amplitude = True
plot_windows   = True
plot_up_dn     = True
# If true, plot will include full simulated domain, if false, just the display domain
plot_full_domain = True

# # Measurements
# take_ef_comp  = False # Energy flux terms recorded separately
# # Records snapshots of total vertical energy flux
# take_ef_snaps = False # Total energy flux recorded

###############################################################################
# Background profile in N_0
n_layers = 1
layer_th = mL/m
L       = layer_th

###############################################################################
# Buffers and important depths

# Buffers
#   distance between ends of vertical structure extent and other points
dis_buff    = 2*lam_z           # [m] buffer from vertical structure to display domain
IT_buff     = dis_buff/2.0      # [m] buffer from vertical structure to measure I and T

# Depths
z0_str  = z0_dis - dis_buff     # [m] top of vertical structure
z_I     = z0_str + IT_buff      # [m] depth at which to measure Incident wave
zf_str  = z0_str - L            # [m] bottom of vertical structure
z_T     = zf_str - IT_buff      # [m] depth at which to measure Transmitted wave
zf_buff = dis_buff+(lam_z-L)    # [m] extra buffer to make sure dis domain is integer number of lam_z
zf_dis  = zf_str - zf_buff      # [m] bottom of display domain

###############################################################################
# Boundary forcing window parameters
b_bf    = 1*lam_z               # [m] full width at half max of forcing window
a_bf    = 3*b_bf                # [m] forcing area, height above display domain
c_bf    = z0_dis + 0.5*a_bf     # [m] center of forcing area
tau_bf  = 1.0e-0                # [s] time constant for boundary forcing

# Sponge layer window parameters
b_sp    = 1*lam_z               # [m] full width at half max of sponge window
a_sp    = 3*b_sp                # [m] sponge area, height below display domain
c_sp    = zf_dis - 0.5*a_sp     # [m] center of sponge area
tau_sp  = 1.0e-0                # [s] time constant for sponge layer

###############################################################################
# Simulated domain parameters
z0     = z0_dis + a_bf          # [m] top of simulated domain
zf     = zf_dis - a_sp          # [m] bottom of simulated domain
Lz     = abs(zf - z0)           # [m] length of simulated domain
dz     = Lz / nz                # [m] spacing between each grid point
dealias= 3/2                    # [] dealiasing factor

# Bases and domain
z_basis = de.Fourier('z', nz, interval=(z0, zf), dealias=dealias)
domain = de.Domain([z_basis], grid_dtype=np.complex128)#float64)
# Z grid
z_da = domain.grid(0, scales=domain.dealias)
z = domain.grid(0)
# Getting wavenumbers
ks = z_basis.wavenumbers

###############################################################################
# Windows for simulation

# Boundary forcing window 2
if boundary_forcing_region == True:
    # the -4*ln(2) is to insure the FWHM is easily identifiable
    win_bf_array = np.exp(-4*np.log(2)*((z - c_bf)/b_bf)**2)     # Gaussian
else:
    win_bf_array = z*0.0 + 1.0  # Forcing over entire domain
    use_sponge = False          # Having sp layer and full domain forcing causes problems

# Sponge layer window 2
if use_sponge == True:
    # the -4*ln(2) is to insure the FWHM is easily identifiable
    win_sp_array = np.exp(-4*np.log(2)*((z - c_sp)/b_sp)**2)     # Gaussian
else:
    win_sp_array = z * 0.0      # No sponge

###############################################################################
# Run parameters
dt              = T/o_steps     # [s] initial time step size (should be around 0.125)
snap_dt         = 32*dt         # [s] time step size for snapshots
snap_max_writes = 100           # [] max number of writes per snapshot file
fh_mode         = 'overwrite'   # file handling mode, either 'overwrite' or 'append'
# Stopping conditions for the simulation
sim_time_stop = T*(T_cutoff+n_T)# [s] number of simulated seconds until the sim stops
stop_wall_time = 180 * 60.0     # [s] length in minutes * 60 = length in seconds, sim stops if exceeded
stop_iteration = np.inf         # [] number of iterations before the simulation stops

###############################################################################
# Physical parameters
nu          = 1.0E-6        # [m^2/s] Viscosity (momentum diffusivity)
kappa       = 1.4E-7        # [m^2/s] Thermal diffusivity
Prandtl     = nu / kappa    # [] Prandtl number, about 7.0 for water at 20 C
Rayleigh    = 1e6
g           = 9.81          # [m/s^2] Acceleration due to gravity

###############################################################################
# Plotting parameters

# Dark mode on or off (ideally would make plots that have white text and alpha background)
dark_mode = False
cmap = 'RdBu_r'
# import colorcet as cc
# cmap = cc.CET_D4

# Presentation mode on or off (increases size of fonts and contrast of colors)
presenting = False

# Vertical profile and Wave field animation
# If True, plots b, p, u, and w. If false, plots profile and w
plot_all_variables = False
# If True, the sponge layer plot will be plotted to the right of the animation
plot_sponge        = False
# If True, the Rayleigh friction plot will replace background profile
plot_rf            = False
plot_twin          = False

# Auxiliary snapshot plots
plot_ef_total = False
plot_ef_comps = False

# Miscellaneous
# Fudge factor to make plots look nicer
buffer = 0.04
# Extra buffer for a constant vertical profile
extra_buffer = 0.5
# Display ratio of vertical profile plot
vp_dis_ratio = 2.0 # Profile plot gets skinnier as this goes up
# The number of ticks on the top color bar
n_clrbar_ticks = 3
# Overall font size of plots
font_size   = 12
scale       = 2.5
dpi         = 100

# Animation parameters
fps = 20

###############################################################################
# Snapshot parameters
snapshots_dir   = 'snapshots'
snap_dt         = 0.25
snap_max_writes = 25

# Background profile snapshot parameters
take_bp_snaps   = True
# Sponge layer snapshot parameters
take_sl_snaps   = True
# Rayleigh friction snapshot parameters
take_rf_snaps   = True

# Define all vertical profile snapshots in an array of dictionaries
#   Meant for profiles that are constant in time
take_vp_snaps = True
vp_snap_dir = 'vp_snapshots'
vp_snap_dicts = [
           {'take_vp_snaps':   take_bp_snaps,
            'vp_task':         "N0*BP",
            'vp_task_name':    'bp'},

           {'take_vp_snaps':   take_sl_snaps,
            'vp_task':         "SL",
            'vp_task_name':    'sl'},

           {'take_vp_snaps':   take_rf_snaps,
            'vp_task':         "RF",
            'vp_task_name':    'rf'}
            ]

# Auxiliary snapshot directory
aux_snap_dir = 'aux_snapshots'

###############################################################################
# Initial conditions
psi_initial = 0.0

###############################################################################
# CFL parameters
CFL_cadence     = 10
CFL_safety      = 1
CFL_max_change  = 1.5
CFL_min_change  = 0.5
CFL_max_dt      = 0.125
CFL_threshold   = 0.05
###############################################################################
# Flow properties
flow_cadence    = 10
flow_property   = "(k*u + m*w)/omega"
flow_name       = 'Lin_Criterion'
###############################################################################
# Logger parameters
endtime_str     = 'Sim end period: %f'
logger_cadence  = 100
iteration_str   = 'Iteration: %i, t/T: %e, dt/T: %e'
flow_log_message= 'Max linear criterion = {0:f}'
###############################################################################
