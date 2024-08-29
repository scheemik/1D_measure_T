# Switchboard for Dedalus experiment
"""
Author: Mikhail Schee

Description:
This is the switchboard for the Dedalus experiment. This file contains parameter values called by multiple scripts. This insures they all get the same values
"""

import numpy as np
from dedalus import public as de
# import sys
# sys.path.append("../") # Adds higher directory to python modules path
import params

###############################################################################
# Main parameters, the ones I'll change a lot. Many more below

# Problem parameters
N_0     = 1.0                   # [rad/s]   Reference stratification
kL      = params.kL             # []        Vertical wave number times step length

theta   = params.theta          # [rad]     Propagation angle from horizontal (if arctan(k/m))
lam_z   = 1.0                   # [m]       Vertical wavelength (used as vertical scale)

A       = 2.0e-4                # []        Amplitude of boundary forcing
f_0     = 0.000                 # [s^-1]    Reference Coriolis parameter

# Background profile in N_0
n_layers = 1                    # []        Number of mixed layers in the background stratification
interface_ratio = 1.0           # []        Ratio between the thickness of an interface to a layer

# Vertical space parameters (z)
nz      = 1024                  # [] number of grid points in z display domain (expecting power of 2)

# Time parameters
p_T_keep    = 3                 # [] number of steady state oscillations to keep = 2 ** p_T_keep
p_o_steps   = 6                 # [] timesteps per oscillation period = 2 ** p_o_steps

# Domain parameters
z0_dis      = 0.0               # [m] Top of the displayed z domain
dealias     = 3/2               # []  Dealiasing factor
snap_rate   = 4                 # []  Snapshots record every `snap_rate` timesteps

###############################################################################
# Physical parameters
nu          = 1.0E-6        # [m^2/s]   Viscosity (momentum diffusivity)
kappa       = 1.4E-7        # [m^2/s]   Thermal diffusivity
Prandtl     = nu / kappa    # []        Prandtl number, about 7.0 for water at 20 C
Rayleigh    = 1e6           # []        Rayleigh number
g           = 9.81          # [m/s^2]   Acceleration due to gravity

###############################################################################
# ON / OFF Switches

# Determine whether adaptive time stepping is on or off
adapt_dt                = False

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
if boundary_forcing_region == False:
    # Having sp layer and full domain forcing causes problems
    use_sponge = False

# Plotting parameters
plot_spacetime  = True
plot_spectra    = False
plot_amplitude  = False
plot_windows    = True
plot_up_dn      = False
plot_untrimmed  = True
# If true, plot will include full simulated domain, if false, just the display domain
plot_full_x     = True
plot_full_y     = True
dark_mode       = True

###############################################################################
###############################################################################
# Calculate wave parameters
def calc_wave_params(N_0, theta, lam_z):
    omega   = N_0 * np.cos(theta)   # [rad s^-1]    Wave frequency
    m       = 2*np.pi / lam_z       # [m^-1]        Vertical wavenumber
    k       = m/np.tan(theta)       # [m^-1]        Horizontal wavenumber
    k_total = np.sqrt(k**2 + m**2)  # [m^-1]        Total wavenumber
    lam_x   = 2*np.pi / k           # [m]           Horizontal wavelength
    return omega, m, k, k_total, lam_x
omega, m, k, k_total, lam_x = calc_wave_params(N_0, theta, lam_z)

def calc_period(omega):
    T       = 2*np.pi / omega       # [s]           Wave period
    return T
T = calc_period(omega)

def calc_wave_speeds(omega, k, m):
    # Cushman-Roisin and Beckers eq 13.8 pg 399
    c_ph = omega / np.sqrt(k**2 + m**2)     # [m/s] Phase speed
    # Cushman-Roisin and Beckers eq 13.9 pg 400
    c_gx = omega*m**2 / (k*(k**2 + m**2))   # [m/s] Horizontal group speed
    # Cushman-Roisin and Beckers eq 13.9 pg 400
    c_gz = - omega*m / (k**2 + m**2)        # [m/s] Vertical group speed
    return c_ph, c_gx, c_gz
c_ph, c_gx, c_gz = calc_wave_speeds(omega, k, m)

###############################################################################
# Calculate background profile in N_0
def calc_layer_thickness(kL, k):
    L = kL/k                        # [m]       Layer thickness
    return L
L = calc_layer_thickness(kL, k)

def calc_R_i(interface_ratio, n_layers):
    # If there are no mixed layers,
    #   set interface ratio to zero to avoid negative lengths in structure
    if n_layers == 0:
        return 0.0
    else:
        return interface_ratio      # []        Ratio between layer and interface thicknesses
R_i = calc_R_i(interface_ratio, n_layers)
###############################################################################
# Calculate buffers and important depths
def calc_structure_depths(z0_dis, lam_z, L, n_layers, R_i):
    # Find distance between ends of vertical structure extent and other points
    dis_buff    = 2*lam_z           # [m] buffer from vertical structure to display domain
    IT_buff     = dis_buff/2.0      # [m] buffer from vertical structure to measure I and T

    # The total length of the vertical structure
    Lz_str  = n_layers*L + (n_layers-1)*R_i*L
    # Find a correction (structure buffer)
    #   to ensure the display domain is an integer number of lam_z
    foo      = Lz_str / lam_z
    str_buff = (1- (foo - int(foo))) * lam_z

    # Depths
    z0_str  = z0_dis - dis_buff     # [m] top of vertical structure
    z_I     = z0_str + IT_buff      # [m] depth at which to measure Incident wave
    zf_str  = z0_str - Lz_str       # [m] bottom of vertical structure
    z_T     = zf_str - IT_buff      # [m] depth at which to measure Transmitted wave
    zf_buff = dis_buff + str_buff   # [m] add extra buffer so dis domain is integer number of lam_z
    zf_dis  = zf_str - zf_buff      # [m] bottom of display domain
    Lz_dis = abs(zf_dis - z0_dis)   # [m] length of display domain
    # Return relevant depths
    return z_I, z0_str, zf_str, z_T, zf_dis, Lz_dis
z_I, z0_str, zf_str, z_T, zf_dis, Lz_dis = calc_structure_depths(z0_dis, lam_z, L, n_layers, R_i)

###############################################################################
# Calculate boundary forcing window parameters
def calc_bf_win_params(zE_dis, t_or_b, lam_z, a, b):
    """
    zE_dis          One edge of the display domain (top or bottom)
    t_or_b          Either 1 for top of display domain or -1 for bottom
    """
    b_bf    = b*lam_z                   # [m] full width at half max of gaussian forcing window
    a_bf    = a*b_bf                    # [m] forcing extent, outside display domain
    c_bf    = zE_dis + t_or_b*a_bf/2    # [m] center of forcing area, t_or_b is +1 or -1
    return a_bf, b_bf, c_bf
# Boundary forcing from above
a_bf, b_bf, c_bf = calc_bf_win_params(z0_dis,  1, lam_z, 3, 1)
# Sponge layer below
a_sp, b_sp, c_sp = calc_bf_win_params(zf_dis, -1, lam_z, 3, 1)

# Boundary forcing window parameters
tau_bf  = 1.0e-0                # [s] time constant for boundary forcing
# Sponge layer window parameters
tau_sp  = 1.0e-0                # [s] time constant for sponge layer

###############################################################################
# Calculate simulation timing

# Calculate the transient time of the simulation based on vertical group speed
def calc_t_tr(c_bf, z0_dis, zf_dis, c_gz):
    # Calculate distance a wave would propagate during transient time
    #   Assuming wave comes from generation point, travels down to
    #   the bottom of the display domain, and back to the top of
    #   the display domain
    d_prop = abs(z0_dis - c_bf) + 2*abs(zf_dis - z0_dis)
    # Use the vertical group speed to find transient time (require a positive number)
    t_tr = abs(d_prop / c_gz)
    return t_tr
t_tr = calc_t_tr(c_bf, z0_dis, zf_dis, c_gz)

def calc_tr_timesteps(t_tr, T, p_o_steps):
    # The number of time steps required for the transient period at the
    #   start of the simulation
    # Find number of oscillation periods that fit in the transient period (doesn't need to be an integer)
    T_tr    = t_tr / T
    # Find number of time steps per oscillation period
    o_steps = int(2**p_o_steps)
    # Find number of time steps in the transient period (needs to be an integer)
    nt_tr   = int(T_tr * o_steps)
    return nt_tr, o_steps, T_tr
nt_tr, o_steps, T_cutoff = calc_tr_timesteps(t_tr, T, p_o_steps)

def calc_keep_timesteps(p_T_keep, p_o_steps):
    # The number of time steps to be kept at the end of the simulation
    #   as 2^(p_T_keep) * 2^(p_o_steps)
    nt_keep     = int(2**(p_T_keep + p_o_steps))
    return nt_keep
nt_keep = calc_keep_timesteps(p_T_keep, p_o_steps)

def calc_t_keep(nt_keep, p_o_steps, T):
    # Calculate the simulation time in the period I want to keep
    # Find the number of time steps per oscillation period
    o_steps = int(2**p_o_steps)
    # Time is period time divided by steps/period times steps
    t_keep = nt_keep * T / o_steps
    return t_keep
t_keep = calc_t_keep(nt_keep, p_o_steps, T)

def calc_timesteps(nt_tr, nt_keep):
    # Total time steps in the simulation equals the timesteps in the transient period
    #   plus the timestepst that I want to keep at the end
    n_steps = nt_tr + nt_keep
    return n_steps
n_steps = calc_timesteps(nt_tr, nt_keep)

###############################################################################
# Simulated domain parameters
def calc_sim_domain(z0_dis, zf_dis, a_bf, a_sp):
    z0     = z0_dis + a_bf          # [m] top of simulated domain
    zf     = zf_dis - a_sp          # [m] bottom of simulated domain
    Lz     = abs(zf - z0)           # [m] length of simulated domain
    return z0, zf, Lz
z0, zf, Lz = calc_sim_domain(z0_dis, zf_dis, a_bf, a_sp)

def calc_nz_sim(nz, Lz, Lz_dis):
    dz     = Lz_dis / nz            # [m] spacing between each grid point
    nz_sim = int(Lz/dz)             # []  number of points in the simulated domain
    # make sure nz_sim is an even number
    if nz_sim%2 != 0:
        nz_sim += 1
    return nz_sim, dz
nz_sim, dz = calc_nz_sim(nz, Lz, Lz_dis)

# Bases and domain
grid_data_type = np.complex128
z_basis = de.Fourier('z', nz_sim, interval=(z0, zf), dealias=dealias)
domain = de.Domain([z_basis], grid_dtype=grid_data_type)#float64)
# Z grid
z_da = domain.grid(0, scales=domain.dealias)
z = domain.grid(0)
# Getting wavenumbers
ks = z_basis.wavenumbers

###############################################################################
# Windows for simulation

# Boundary forcing window 2
def calc_bf_array(z, c_bf, b_bf, boundary_forcing_region):
    if boundary_forcing_region == True:
        # the -4*ln(2) is to insure the FWHM is easily identifiable
        win_bf_array = np.exp(-4*np.log(2)*((z - c_bf)/b_bf)**2)     # Gaussian
    else:
        win_bf_array = z*0.0 + 1.0  # Forcing over entire domain
    return win_bf_array
win_bf_array = calc_bf_array(z, c_bf, b_bf, boundary_forcing_region)

# Sponge layer window 2
def calc_sp_array(z, c_sp, b_sp, use_sponge):
    if use_sponge == True:
        # the -4*ln(2) is to insure the FWHM is easily identifiable
        win_sp_array = np.exp(-4*np.log(2)*((z - c_sp)/b_sp)**2)     # Gaussian
    else:
        win_sp_array = z * 0.0      # No sponge
    return win_sp_array
win_sp_array = calc_sp_array(z, c_sp, b_sp, use_sponge)

###############################################################################
# Run parameters
dt              = T/o_steps     # [s] initial time step size (should be around 0.125)
snap_dt         = dt*snap_rate  # [s] time step size for snapshots
nt_snap    = nt_keep//snap_rate # []  number of time steps in the snapshots
snap_max_writes = 100           # []  max number of writes per snapshot file
snapshots_dir   = 'snapshots'   # name of directory in which to put snapshot files
fh_mode         = 'overwrite'   # file handling mode, either 'overwrite' or 'append'
# Stopping conditions for the simulation
sim_time_stop  = t_tr + t_keep  # [s] number of simulated seconds until the sim stops
stop_wall_time = 180 * 60.0     # [s] length in minutes * 60 = length in seconds, sim stops if exceeded
stop_iteration = np.inf         # []  number of iterations before the simulation stops

###############################################################################
# Set plotting parameters
import colorcet as cc
if dark_mode:
    plt_style = 'dark_background'
    cmap      = cc.cm.bkr
else:
    plt_style = 'default'
    cmap      = 'RdBu_r' # cc.cm.CET_D9

# Presentation mode on or off (increases size of fonts and contrast of colors)
presenting = False

# # Miscellaneous
# # Fudge factor to make plots look nicer
# buffer = 0.04
# # Extra buffer for a constant vertical profile
# extra_buffer = 0.5
# # Display ratio of vertical profile plot
# vp_dis_ratio = 2.0 # Profile plot gets skinnier as this goes up
# # The number of ticks on the top color bar
# n_clrbar_ticks = 3
# # Overall font size of plots
# font_size   = 12
# scale       = 2.5
# dpi         = 100

# Animation parameters
fps = 20

###############################################################################
# Auxiliary snapshot parameters

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
