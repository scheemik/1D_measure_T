"""
Performs post-processing actions. Run with $ python3 post_process.py snapshots/*.h5

Usage:
    post_process.py <files>... [--output=<dir>]

Options:
    <files>         # h5 snapshot files

"""

import h5py
import numpy as np
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
h5_files = args['<files>']

###############################################################################
# Import SwitchBoard Parameters (sbp)
#   This import assumes the switchboard is in the same directory as the core code
import switchboard as sbp
# Physical parameters
nu          = sbp.nu            # [m^2/s] Viscosity (momentum diffusivity)
f_0         = sbp.f_0           # [s^-1]        Reference Coriolis parameter
g           = sbp.g             # [m/s^2] Acceleration due to gravity
# Problem parameters
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
tasks = ['psi']
skip_nT = 0
nz = sbp.nz
z0_dis = sbp.z0_dis
zf_dis = sbp.zf_dis
dt = sbp.dt
dz = sbp.dz
plot_f_d=sbp.plot_full_domain
nT=sbp.nT
n_steps = sbp.n_steps
step_th = sbp.step_th

###############################################################################
# Helper functions
#   This import assumes the helper functions are in the same directory as the core code
import helper_functions as hf

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
                # Just plotting the real part for now
                #psi_real = psi_array.real
                # Grab the scales t, z, and kz
                t  = np.array(f['scales']['sim_time'])
                z  = np.array(f['scales']['z']['1.0'])
                kz = np.array(f['scales']['kz'])
    return t, z, kz, psi_array#, psi_c_array

# fourier transform in time, filter negative freq's, inverse fourier transform
def FT_in_time(t, z, data, dt):
    # FT in time of the data (axis 1 is time)
    ftd = np.fft.fft(data, axis=1)
    # find relevant frequencies
    freq = np.fft.fftfreq(len(t), dt)
    # Filter out negative frequencies
    for j in range(len(freq)):
        if freq[j] < 0.0:
            # Gets rid of negative freq's
            ftd[:,j] = 0.0
        else:
            # Corrects for lost amplitude
            ftd[:,j] = ftd[:,j] * 2.0
    # inverse fourier transform in time of the data
    iftd = np.fft.ifft(ftd, axis=1)
    #   a complex valued signal where iftd.real == data, or close enough
    return iftd

# fourier transform in spatial dimension (z)
#   similar to FT in time, but switch dimensions around
def FT_in_space(t, k_zs, data):
    # FT in space (z) of the data (axis 0 is z) for positive wave numbers
    fzdp = np.fft.fft(data, axis=0)
    # make a copy for the negative wave numbers
    fzdn = fzdp.copy()
    # Filter out one half of wavenumbers to separate up and down
    #   Looping only over wavenumbers because their positions don't change with t
    for i in range(len(k_zs)):#k_grid.shape[0]):
        if k_zs[i] > 0.0:
            # for down, remove values for positive wave numbers
            fzdn[i,:] = 0.0
        else:
            # for up, remove values for negative wave numbers
            fzdp[i,:] = 0.0
    # inverse fourier transform in space (z)
    ifzdp = np.fft.ifft(fzdp, axis=0)
    ifzdn = np.fft.ifft(fzdn, axis=0)
    return ifzdp, ifzdn

# fourier transform in spatial dimension (z)
#   similar to FT in time, but switch dimensions around
def IFT_in_space(t, k_zs, data_up):
    # make a copy for the negative wave numbers
    data_dn = data_up.copy()
    # Filter out one half of wavenumbers to separate up and down
    #   Looping only over wavenumbers because their positions don't change with t
    for i in range(len(k_zs)):#k_grid.shape[0]):
        if k_zs[i] > 0.0:
            # for down, remove values for positive wave numbers
            data_dn[i][:] = 0.0
        else:
            # for up, remove values for negative wave numbers
            data_up[i][:] = 0.0
    # inverse fourier transform in space (z)
    ifzdp = np.fft.ifft(data_up, axis=0)
    ifzdn = np.fft.ifft(data_dn, axis=0)
    return ifzdp, ifzdn

###############################################################################
# Get the data from the snapshot files
t, z, kz, psi = get_h5_data(tasks, h5_files)
#c_data = hf.sort_k_coeffs(raw_data, nz)

BP_array = hf.BP_n_steps(sbp.n_steps, sbp.z, z0_dis, zf_dis, sbp.step_th)

###############################################################################
# Sanity check plots

if sbp.plot_spacetime:
    hf.plot_z_vs_t(z, t, T, psi.real, BP_array, k, m, omega, z0_dis, zf_dis, plot_full_domain=sbp.plot_full_domain, nT=sbp.nT)#, title_str=run_name)

# if sbp.plot_wavespace:
#     hf.plot_k_vs_t(z, t, T, psi.real, psi.imag, k, m, omega, title_str='psi', filename='f_1D_psi_r_i.png')
#
# if sbp.plot_wavespace:
#     hf.plot_k_vs_t(hf.sort_k_coeffs(kz,1024), t, T, psi_c.real, psi_c.imag, k, m, omega, title_str='psi_c', filename='f_1D_psic_r_i.png')

if sbp.plot_amplitude:
    hf.plot_A_vs_t(t, T, psi.real, sbp.A, k, m, omega, nT=sbp.nT)#, title_str=run_name)

# raise SystemExit(0)

###############################################################################
# Complex demodulation

def Complex_Demodulate(t_then_z, t, z, kz, data, dt):
    if t_then_z == True:
        ## Step 1
        ift_t_y = FT_in_time(t, z, data, dt)
        ## Step 2
        ift_z_y_p, ift_z_y_n = FT_in_space(t, kz, ift_t_y)
        # Get up and down fields as F = |mag_f| * exp(i*phi_f)
        up_field = ift_z_y_p.real * np.exp(np.real(1j * ift_z_y_p.imag))
        dn_field = ift_z_y_n.real * np.exp(np.real(1j * ift_z_y_n.imag))
    else:
        ## Step 1
        ift_z_y_p, ift_z_y_n = FT_in_space(t, kz, data)
        ## Step 2
        up_f = FT_in_time(t, z, ift_z_y_p, dt)
        dn_f = FT_in_time(t, z, ift_z_y_n, dt)
        # Get up and down fields as F = |mag_f| * exp(i*phi_f)
        up_field = up_f.real * np.exp(np.real(1j * up_f.imag))
        dn_field = dn_f.real * np.exp(np.real(1j * dn_f.imag))
    return up_field, dn_field

t_then_z = False
up_field, dn_field = Complex_Demodulate(t_then_z, t, z, kz, psi, dt)

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
# Measuring the transmission coefficient

big_T = hf.measure_T(dn_field, z, -0.25, -0.75, dz)
print("Transmission coefficient is:", big_T)

###############################################################################
# More plotting for up and down waves

if sbp.plot_up_dn:
    hf.plot_z_vs_t(z, t, T, up_field, BP_array, k, m, omega, z0_dis, zf_dis, plot_full_domain=plot_f_d, nT=nT, filename='f_1D_up_field.png')
    hf.plot_z_vs_t(z, t, T, dn_field, BP_array, k, m, omega, z0_dis, zf_dis, plot_full_domain=plot_f_d, nT=nT, filename='f_1D_dn_field.png')
