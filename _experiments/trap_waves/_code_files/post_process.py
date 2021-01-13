"""
Performs post-processing actions. Run with $ python3 post_process.py NAME PLOT_CHECKS snapshots/*.h5

Usage:
    post_process.py NAME PLOT_CHECKS <files>... [--output=<dir>]

Options:
    NAME            # name of the experiment run from -n
    PLOT_CHECKS     # True or False whether to make the plots or not
    <files>         # h5 snapshot files

"""
###############################################################################
###############################################################################
# Import standard(ish) libraries and functions
import h5py
import numpy as np
import importlib
import matplotlib
#matplotlib.use('Agg')
from matplotlib import ticker
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from dedalus.extras.plot_tools import quad_mesh
#plt.ioff()

# Parse input parameters
from docopt import docopt
args = docopt(__doc__)
run_name    = args['NAME']      # Simulation name, used to route filepaths of plots
plot_checks = args['PLOT_CHECKS'].lower() == 'true'
h5_files    = args['<files>']

###############################################################################
# Import auxiliary files

# Import SwitchBoard Parameters (sbp)
import switchboard as sbp
# switchboard_module = run_name + "." + run_name + "_" + switchboard
# sbp = importlib.import_module(switchboard_module)

import helper_functions as hf
import helper_functions_CD as hfCD

###############################################################################
# Importing parameters from auxiliary files

# Problem parameters
mL          = sbp.mL            # []            Ratio of vertical wavelength to stratification length
k           = sbp.k             # [m^-1]        Horizontal wavenumber
m           = sbp.m             # [m^-1]        Vertical wavenumber
theta       = sbp.theta         # [rad]         Propagation angle from vertical
omega       = sbp.omega         # [rad s^-1]    Wave frequency
T           = sbp.T             # [s]           Wave period

# Parameters
tasks = ['psi']

# Relevant depths
z0_dis      = sbp.z0_dis        # [m]           Top of the displayed z domain
zf_dis      = sbp.zf_dis        # [m]           Bottom of the displayed z domain
z_I         = sbp.z_I           # [m]           depth at which to measure Incident wave
z_T         = sbp.z_T           # [m]           depth at which to measure Transmitted wave

# Grid spacing
dt          = sbp.dt
dz          = sbp.dz

plt_f_x     = sbp.plot_full_x
plt_f_y     = sbp.plot_full_y
T_cutoff    = sbp.T_cutoff
n_layers    = sbp.n_layers
layer_th    = sbp.layer_th

###############################################################################
# Additional post-processing helper functions

def get_h5_data(tasks, h5_files):
    for task in tasks:
        # At the moment, this only works for one h5 file
        for filename in h5_files:
            with h5py.File(filename, mode='r') as f:
                # The [()] syntax returns all data from an h5 object
                psi   = f['tasks'][task]
                # psi_c = f['tasks'][task]
                # Need to transpose into the correct orientation
                #   Also need to convert to np.array for plotting function
                psi_array   = np.transpose(np.array(psi[()]))
                # psi_c_array = np.transpose(np.array(psi_c[()]))
                # Grab the scales t, z, and kz
                t  = np.array(f['scales']['sim_time'])
                z  = np.array(f['scales']['z']['1.0'])
                kz = np.array(f['scales']['kz'])
    return t, z, kz, psi_array#, psi_c_array


###############################################################################
###############################################################################
# Get the data from the snapshot files
t, z, kz, psi = get_h5_data(tasks, h5_files)
# z and psi arrays come out sorted from most positive to most negative on z axis
#   This flips things around (ud = up / down)
z = np.flip(z)
psi = np.flipud(psi)

###############################################################################
# Trim data

# Trim in space
z_tr, psi_tr = hf.trim_data_z(z, psi, z0_dis, zf_dis)

# Trim in time
# print("T_cutoff = ", T_cutoff)
# print("t = ", t)
t_tr, psi_tr = hf.trim_data_t(t, psi_tr, T_cutoff, T)

###############################################################################
# Perform complex demodulation

t_then_z = True
# find relevant wavenumbers (k) for the trimmed z range
kz_tr = np.fft.fftfreq(len(z_tr), dz)
# I don't know why, but the up and down fields switch when I trim the data
tr_dn_field, tr_up_field = hfCD.Complex_Demodulate(t_then_z, t_tr, z_tr, kz_tr, psi_tr, dt, omega)

###############################################################################
# Measuring the transmission coefficient

I_, T_, AAcc_I, AAcc_T = hf.measure_T(tr_dn_field, z_tr, z_I, z_T, T_skip=None, T=T, t=t_tr)
big_T = T_/I_
print("(n_layers =",n_layers,"(mL =",mL,", theta =",theta,")")
print("Simulated transmission coefficient is:", big_T)
print("AnaEq 2.4 transmission coefficient is:", hf.SY_eq2_4(theta, mL))


# Write out results to file
import csv
csv_file = "sim_data.csv"
with open(csv_file, 'a') as datafile:
    csvwriter = csv.writer(datafile)
    csvwriter.writerow([run_name, mL, theta, big_T])

###############################################################################
###############################################################################
###############################################################################
# Plotting checks

if plot_checks == False:
    raise SystemExit(0)

BP_array = hf.BP_n_layers(z, sbp.z0_str, n_layers, sbp.L, sbp.R_i)
#hf.BP_n_layers(n_layers, z, sbp.z0_str, sbp.zf_str)
foo, BP_tr   = hf.trim_data_z(z, BP_array, z0_dis, zf_dis)

filename_prefix = run_name #+ '/' + run_name

# Plot windows
if sbp.plot_windows:
    hf.plot_v_profiles(z, BP_array, sbp.win_bf_array, sbp.win_sp_array, mL, theta, omega, z_I=z_I, z_T=z_T, z0_dis=z0_dis, zf_dis=zf_dis, plot_full_x=True, plot_full_y=True, title_str=run_name, filename=filename_prefix+'_windows.png')

# Plot wavefield
if sbp.plot_spacetime:
    hf.plot_z_vs_t(z, t, T, psi.real, BP_array, mL, theta, omega, z_I=z_I, z_T=z_T, z0_dis=z0_dis, zf_dis=zf_dis, plot_full_x=plt_f_x, plot_full_y=plt_f_y, T_cutoff=T_cutoff, title_str=run_name, filename=filename_prefix+'_wave.png')

# Plot trimmed wavefield
if sbp.plot_spacetime:
    hf.plot_z_vs_t(z_tr, t_tr, T, psi_tr.real, BP_tr, mL, theta, omega, z_I=z_I, z_T=z_T, z0_dis=z0_dis, zf_dis=zf_dis, plot_full_x=True, plot_full_y=False, T_cutoff=T_cutoff, title_str=run_name, filename=filename_prefix+'_wave_tr.png')

###############################################################################
# Plot spectral form of data

freqs, ks, spec_data = hfCD.z_t_to_k_omega(t_tr, z_tr, psi_tr, dt, dz)

if sbp.plot_spectra:
    hf.plot_spectral(ks, freqs, spec_data.real, spec_data.imag, mL, theta, omega, c_map='viridis', title_str=run_name, filename=filename_prefix+'_wave_spectra.png')
    # plot spectra lines
    k_data = np.fft.fft(psi_tr, axis=0)
    f_data = np.fft.fft(psi_tr, axis=1)
    hf.plot_k_f_spectra(z_tr, dz, t_tr, dt, T, ks, freqs, k_data, f_data, mL, theta, omega, z_I, z_T, plot_full_x=True, plot_full_y=True, T_cutoff=T_cutoff+1, title_str=run_name, filename=filename_prefix+'_k_and_f.png')

###############################################################################
# Plotting up and downward propagating waves

# Trimmed case
if sbp.plot_up_dn:
    hf.plot_z_vs_t(z_tr, t_tr, T, tr_up_field.real, BP_tr, mL, theta, omega, z_I=z_I, z_T=z_T, z0_dis=z0_dis, zf_dis=zf_dis,  plot_full_x=plt_f_x, plot_full_y=plt_f_y, T_cutoff=None, title_str=run_name+' up', filename=filename_prefix+'_up_tr.png')
    hf.plot_z_vs_t(z_tr, t_tr, T, tr_dn_field.real, BP_tr, mL, theta, omega, z_I=z_I, z_T=z_T, z0_dis=z0_dis, zf_dis=zf_dis, plot_full_x=plt_f_x, plot_full_y=plt_f_y, T_cutoff=None, title_str=run_name+' dn', filename=filename_prefix+'_dn_tr.png')

up_field, dn_field = hfCD.Complex_Demodulate(t_then_z, t, z, kz, psi, dt, omega)
# Full domain
if sbp.plot_up_dn:
    hf.plot_z_vs_t(z, t, T, up_field.real, BP_array, mL, theta, omega, z_I=z_I, z_T=z_T, z0_dis=z0_dis, zf_dis=zf_dis,  plot_full_x=plt_f_x, plot_full_y=plt_f_y, T_cutoff=T_cutoff, title_str=run_name+' up', filename=filename_prefix+'_up.png')
    hf.plot_z_vs_t(z, t, T, dn_field.real, BP_array, mL, theta, omega, z_I=z_I, z_T=z_T, z0_dis=z0_dis, zf_dis=zf_dis, plot_full_x=plt_f_x, plot_full_y=plt_f_y, T_cutoff=T_cutoff, title_str=run_name+' dn', filename=filename_prefix+'_dn.png')

###############################################################################
# Multiply the downward wavefield by it's complex-conjugate to get AA^*
#   Plot this amplitude across time at two specific depths:
if sbp.plot_amplitude:
    hf.plot_A_of_I_T(z_tr, t_tr, T, tr_dn_field, mL, theta, omega, z_I, z_T, plot_full_x=plt_f_x, plot_full_y=plt_f_y, T_cutoff=T_cutoff, title_str=run_name, filename=filename_prefix+'_A_of_I_T.png')
# Plot this amplitude (averaged across time) as a function of depth
if sbp.plot_amplitude:
    hf.plot_AA_for_z(z_tr, BP_tr, hf.AAcc(tr_dn_field).real, mL, theta, omega, z_I=z_I, z_T=z_T, z0_dis=z0_dis, zf_dis=zf_dis, plot_full_x=plt_f_x, plot_full_y=plt_f_y, title_str=run_name, filename=filename_prefix+'_AA_for_z.png')

if sbp.plot_amplitude:
    hf.plot_z_vs_t(z_tr, t_tr, T, hf.AAcc(tr_dn_field).real, BP_tr, mL, theta, omega, z_I=z_I, z_T=z_T, z0_dis=z0_dis, zf_dis=zf_dis, plot_full_x=plt_f_x, plot_full_y=plt_f_y, T_cutoff=None, title_str=run_name+' AA', filename=filename_prefix+'_1D_AA.png')

raise SystemExit(0)

###############################################################################
# Profiling the code
profile_it = False
if profile_it == True:
    import cProfile
    cProfile.run('Complex_Demodulate(t_then_z, t, z, kz, psi, dt)', 'restats')
    import pstats
    from pstats import SortKey
    p = pstats.Stats('restats')
    p.sort_stats(SortKey.CUMULATIVE).print_stats(10)
###############################################################################

plot_CD_checks = False
if plot_CD_checks:
    # Add up and down fields to see if they reproduce the original psi field
    up_plus_dn = up_field.real + dn_field.real
    hf.plot_z_vs_t(z, t, T, up_plus_dn, BP_array, mL, theta, omega, z0_dis=z0_dis, zf_dis=zf_dis, z_I=z_I, z_T=z_T, plot_full_x=plt_f_x, plot_full_y=plt_f_y, nT=T_cutoff, title_str=run_name, filename=filename_prefix+'_up_plus_dn.png')
    # Plot the difference, which ideally should be zero everywhere
    CD_diff = psi.real - up_plus_dn
    hf.plot_z_vs_t(z, t, T, CD_diff, BP_array, mL, theta, omega, z0_dis=z0_dis, zf_dis=zf_dis, z_I=z_I, z_T=z_T, plot_full_x=plt_f_x, plot_full_y=plt_f_y, nT=T_cutoff, title_str=run_name, filename=filename_prefix+'_CD_diff.png')
