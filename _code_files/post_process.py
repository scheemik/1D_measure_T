"""
Author: Mikhail Schee

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
from dedalus import public as de
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

###############################################################################
# Importing parameters from auxiliary files

# Problem parameters
kL          = sbp.kL            # []            Ratio of horizontal wavelength to stratification length
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
dt          = sbp.snap_dt
dz          = sbp.dz

# Setup parameters
T_cutoff    = sbp.T_cutoff
nt_keep     = sbp.nt_snap
n_layers    = sbp.n_layers

# Plotting parameters
plot_spacetime  = sbp.plot_spacetime
plot_spectra    = sbp.plot_spectra
plot_amplitude  = sbp.plot_amplitude
plot_windows    = sbp.plot_windows
plot_up_dn      = sbp.plot_up_dn
plot_untrimmed  = sbp.plot_untrimmed
plt_f_x         = sbp.plot_full_x
plt_f_y         = sbp.plot_full_y
plt_style       = sbp.plt_style
cmap            = sbp.cmap

###############################################################################
# Additional post-processing helper functions

# This function is built to only work with one h5 file, need to merge before using
#   Also, only works for 1 task that has 1 spatial and 1 temporal dimension
#   It's a pretty specific function
def import_h5_data(hf_files, zf_dis, z0_dis, nt_keep, dtype=sbp.grid_data_type, dealias_f=sbp.dealias):
    with h5py.File(hf_files[0], mode='r') as f:
        # Get task (one and only file i/o here)
        psi_data = f['tasks']['psi']
        # Get dimensions
        t = psi_data.dims[0]['sim_time']
        z = psi_data.dims[1][0]
        # Trim domain of data
        z_tr, t_tr, psi, psi_tr_data = hf.trim_data(z, t, psi_data, z0_dis, zf_dis, nt_keep)
        # Bases and domain
        t_tr_basis = de.Fourier('t', len(t_tr), interval=(t_tr[0], t_tr[-1]), dealias=dealias_f)
        z_tr_basis = de.Fourier('z', len(z_tr), interval=(z_tr[0], z_tr[-1]), dealias=dealias_f)
        domain_tr = de.Domain([t_tr_basis,z_tr_basis], grid_dtype=dtype)
        # set field
        psi_tr = domain_tr.new_field(name = 'psi_tr')
        psi_tr['g'] = psi_tr_data

        if plot_untrimmed:
            # Bases and domain
            t_basis = de.Fourier('t', len(t), interval=(t[0], t[-1]), dealias=dealias_f)
            z_basis = de.Fourier('z', len(z), interval=(z[0], z[-1]), dealias=dealias_f)
            domain = de.Domain([t_basis,z_basis], grid_dtype=dtype)
            # set field
            psi = domain.new_field(name = 'psi')
            psi['g'] = psi_data
        else:
            psi    = None
            domain = None
        z = np.flip(np.array(z))
    return psi_tr, domain_tr, psi, domain, t, z

def Complex_Demodulation(ft, kz, psi):
    # Filter in time
    n_ft = ft.shape[0] # ft is a 2D array with 1 column
    # Add 1 because Nyquist is dropped
    if hf.is_power_of_2(n_ft+1) == False:
        print('Warning: n_ft is not a power of 2:',n_ft)
        print('  Are you plotting untrimmed data?')
    # Need to copy so I leave the original field untouched
    psi_up = psi.copy()
    # Not sure why, but copying sets scales to (1.5,1.50)
    psi_up.set_scales((1.0,1.0))
    psi.set_scales((1.0,1.0))
    # Remove all negative frequencies
    selection1 = (ft<0) * (kz==kz)
    psi_up['c'][selection1] = 0
    # Multiply remaining by 2 to compensate for lost magnitude
    double_op = 2*psi_up
    double_op.operate(psi_up)

    # Filter in space
    n_kz = kz.shape[1] # kz is a 2D array with 1 row
    # Add 1 because Nyquist is dropped
    if hf.is_power_of_2(n_kz+1) == False:
        print('Warning: n_kz is not a power of 2:',n_kz)
        print('  Are you plotting untrimmed data?')
    psi_dn = psi_up.copy()
    # Not sure why, but copying sets scales to (1.5,1.50)
    psi_up.set_scales((1.0,1.0))
    psi_dn.set_scales((1.0,1.0))
    # Filter for the up propagating waves
    selection2 = (ft==ft) * (kz<0)
    psi_up['c'][selection2] = 0
    # Filter for the down propagating waves
    selection3 = (ft==ft) * (kz>0)
    psi_dn['c'][selection3] = 0
    return psi_up, psi_dn

def get_domain_scales(domain):
    # For grid space
    t = domain.grid(0)[:,0]
    z = domain.grid(1)[0,:]
    # For coefficient space
    f = domain.elements(0)
    k = domain.elements(1)
    return t, z, f, k

def get_plt_field(field):
    return np.flipud(np.transpose(field['g'].real))

# raise SystemExit(0)
###############################################################################
###############################################################################
# Get the data from the snapshot files
psi_tr, domain_tr, psi, domain, t, z = import_h5_data(h5_files, zf_dis, z0_dis, nt_keep)

# Perform complex demodulation
t_tr, z_tr, f_tr, k_tr = get_domain_scales(domain_tr)
print(z_tr)
tr_psi_up, tr_psi_dn = Complex_Demodulation(f_tr, k_tr, psi_tr)

###############################################################################
# Get some arrays ready for measuring T
plt_tr_dn = get_plt_field(tr_psi_dn)
plt_tr_z = np.flip(z_tr[:])
plt_tr_t = t_tr[:]

###############################################################################
# Measuring the transmission coefficient
I_, T_, AAcc_I, AAcc_T = hf.measure_T(plt_tr_dn, plt_tr_z, z_I, z_T, T_skip=None, T=T, t=plt_tr_t)
big_T = T_/I_
print("(n_layers =",n_layers,", kL =",kL,", theta =",theta,")")
print("Simulated transmission coefficient is:", big_T)
print("AnaEq 2.4 transmission coefficient is:", hf.SY_eq2_4(theta, kL))

# Write out results to file
import csv
csv_file = "sim_data.csv"
with open(csv_file, 'a') as datafile:
    csvwriter = csv.writer(datafile)
    csvwriter.writerow([run_name, kL, theta, big_T])

###############################################################################
###############################################################################
###############################################################################
# Plotting checks

if plot_checks == False:
    raise SystemExit(0)

plt.style.use(plt_style)

plt_z = z[:]

BP_array = hf.BP_n_layers(plt_z, sbp.z0_str, n_layers, sbp.L, sbp.R_i)
win_bf_array = sbp.calc_bf_array(plt_z, sbp.c_bf, sbp.b_bf, sbp.boundary_forcing_region)
win_sp_array = sbp.calc_sp_array(plt_z, sbp.c_sp, sbp.b_sp, sbp.use_sponge)
foo, BP_tr   = hf.trim_data_z(z, BP_array, z0_dis, zf_dis)

filename_prefix = run_name #+ '/' + run_name

# Plot windows
if sbp.plot_windows:
    hf.plot_v_profiles(plt_z, BP_array, win_bf_array, win_sp_array, kL, theta, omega, z_I=z_I, z_T=z_T, z0_dis=z0_dis, zf_dis=zf_dis, plot_full_x=True, plot_full_y=True, title_str=run_name, filename=filename_prefix+'_windows.png')

###############################################################################
# Plotting with trimmed data

# Plot trimmed wavefield
if sbp.plot_spacetime:
    plt_psi_tr = get_plt_field(psi_tr)
    hf.plot_z_vs_t(plt_tr_z, plt_tr_t, T, plt_psi_tr, BP_tr, kL, theta, omega, z_I=z_I, z_T=z_T, z0_dis=z0_dis, zf_dis=zf_dis, plot_full_x=True, plot_full_y=False, T_cutoff=T_cutoff, c_map=cmap, title_str=run_name, filename=filename_prefix+'_wave_tr.png')

# Plot trimmed up and downward propagating waves
if sbp.plot_up_dn:
    plt_tr_up = get_plt_field(tr_psi_up)
    hf.plot_z_vs_t(plt_tr_z, plt_tr_t, T, plt_tr_up, BP_tr, kL, theta, omega, z_I=z_I, z_T=z_T, z0_dis=z0_dis, zf_dis=zf_dis,  plot_full_x=plt_f_x, plot_full_y=plt_f_y, T_cutoff=None, title_str=run_name+' up', c_map=cmap, filename=filename_prefix+'_up_tr.png')
    hf.plot_z_vs_t(plt_tr_z, plt_tr_t, T, plt_tr_dn, BP_tr, kL, theta, omega, z_I=z_I, z_T=z_T, z0_dis=z0_dis, zf_dis=zf_dis, plot_full_x=plt_f_x, plot_full_y=plt_f_y, T_cutoff=None, title_str=run_name+' dn', c_map=cmap, filename=filename_prefix+'_dn_tr.png')

###############################################################################
# Plotting with untrimmed data

if sbp.plot_untrimmed:
    t, z, f, k = get_domain_scales(domain)
    plot_psi = get_plt_field(psi)
    # Plot untrimmed wavefield
    if sbp.plot_spacetime:
        hf.plot_z_vs_t(plt_z, t[:], T, plot_psi, BP_array, kL, theta, omega, c_gz=sbp.c_gz, c_bf=sbp.c_bf, z_I=z_I, z_T=z_T, z0_dis=z0_dis, zf_dis=zf_dis, plot_full_x=plt_f_x, plot_full_y=plt_f_y, T_cutoff=T_cutoff, c_map=cmap, title_str=run_name, filename=filename_prefix+'_wave.png')

    # Plot untrimmed up and downward propagating waves
        # This will lead to warning about nt not being a power of 2
    # Full domain
    if sbp.plot_up_dn:
        psi_up, psi_dn = Complex_Demodulation(f, k, psi)
        plt_up = get_plt_field(psi_up)
        plt_dn = get_plt_field(psi_dn)
        # up_field, dn_field = hfCD.Complex_Demodulate(t_then_z, t, z, kz, plot_psi, dt, omega)
        hf.plot_z_vs_t(plt_z, t, T, plt_up, BP_array, kL, theta, omega, z_I=z_I, z_T=z_T, z0_dis=z0_dis, zf_dis=zf_dis,  plot_full_x=plt_f_x, plot_full_y=plt_f_y, T_cutoff=T_cutoff, c_map=cmap, title_str=run_name+' up', filename=filename_prefix+'_up.png')
        hf.plot_z_vs_t(plt_z, t, T, plt_dn, BP_array, kL, theta, omega, z_I=z_I, z_T=z_T, z0_dis=z0_dis, zf_dis=zf_dis, plot_full_x=plt_f_x, plot_full_y=plt_f_y, T_cutoff=T_cutoff, c_map=cmap, title_str=run_name+' dn', filename=filename_prefix+'_dn.png')

###############################################################################
# Plot spectral form of data

if sbp.plot_spectra:
    # plt.clf()
    from dedalus.extras.plot_tools import plot_bot_2d
    figkw = {'figsize':(6,4), 'dpi':100}
    psi_tr['c']
    # Plot log magnitude of spectral coefficients
    log_mag = lambda xmesh, ymesh, data: (xmesh, ymesh, np.log10(np.abs(data)))
    plot_bot_2d(psi_tr, func=log_mag, clim=(-20, 0), cmap='viridis', title="log10(abs(f['c'])", figkw=figkw);
    # plt.show()
    plt.savefig(filename_prefix+'_k_and_f.png')
    # freqs, ks, spec_data = hfCD.z_t_to_k_omega(t_tr, z_tr, psi_tr, dt, dz)
    # hf.plot_spectral(ks, freqs, spec_data.real, spec_data.imag, kL, theta, omega, c_map='viridis', title_str=run_name, filename=filename_prefix+'_wave_spectra.png')
    # # plot spectra lines
    # k_data = np.fft.fft(psi_tr, axis=0)
    # f_data = np.fft.fft(psi_tr, axis=1)
    # hf.plot_k_f_spectra(z_tr, dz, t_tr, dt, T, ks, freqs, k_data, f_data, kL, theta, omega, z_I, z_T, plot_full_x=True, plot_full_y=True, T_cutoff=T_cutoff+1, title_str=run_name, filename=filename_prefix+'_k_and_f.png')

###############################################################################
# Multiply the downward wavefield by it's complex-conjugate to get AA^*

# Plot this amplitude across time at two specific depths:
if sbp.plot_amplitude:
    hf.plot_A_of_I_T(z_tr, t_tr, T, tr_dn_field, kL, theta, omega, z_I, z_T, plot_full_x=plt_f_x, plot_full_y=plt_f_y, T_cutoff=T_cutoff, title_str=run_name, filename=filename_prefix+'_A_of_I_T.png')

# Plot this amplitude (averaged across time) as a function of depth
if sbp.plot_amplitude:
    hf.plot_AA_for_z(z_tr, BP_tr, hf.AAcc(tr_dn_field).real, kL, theta, omega, z_I=z_I, z_T=z_T, z0_dis=z0_dis, zf_dis=zf_dis, plot_full_x=plt_f_x, plot_full_y=plt_f_y, title_str=run_name, filename=filename_prefix+'_AA_for_z.png')

# Plot this amplitude across time and space
if sbp.plot_amplitude:
    hf.plot_z_vs_t(z_tr, t_tr, T, hf.AAcc(tr_dn_field).real, BP_tr, kL, theta, omega, z_I=z_I, z_T=z_T, z0_dis=z0_dis, zf_dis=zf_dis, plot_full_x=plt_f_x, plot_full_y=plt_f_y, T_cutoff=None, c_map=cmap, title_str=run_name+' AA', filename=filename_prefix+'_1D_AA.png')

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
    hf.plot_z_vs_t(z, t, T, up_plus_dn, BP_array, kL, theta, omega, z0_dis=z0_dis, zf_dis=zf_dis, z_I=z_I, z_T=z_T, plot_full_x=plt_f_x, plot_full_y=plt_f_y, nT=T_cutoff, title_str=run_name, filename=filename_prefix+'_up_plus_dn.png')
    # Plot the difference, which ideally should be zero everywhere
    CD_diff = psi.real - up_plus_dn
    hf.plot_z_vs_t(z, t, T, CD_diff, BP_array, kL, theta, omega, z0_dis=z0_dis, zf_dis=zf_dis, z_I=z_I, z_T=z_T, plot_full_x=plt_f_x, plot_full_y=plt_f_y, nT=T_cutoff, title_str=run_name, filename=filename_prefix+'_CD_diff.png')
