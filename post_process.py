"""
Performs post-processing actions. Run with $ python3 post_process.py snapshots/*.h5

Usage:
    post_process.py NAME <files>... [--output=<dir>]

Options:
    NAME        # name of the experiment run from -n
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
run_name = args['NAME']
h5_files = args['<files>']

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
tasks = ['psi']
nz          = sbp.nz
z0_dis      = sbp.z0_dis
zf_dis      = sbp.zf_dis

dt          = sbp.dt
dz          = sbp.dz

plot_f_d    = sbp.plot_full_domain
T_skip      = sbp.T_skip
n_steps     = sbp.n_steps
step_th     = sbp.step_th

z_I         = sbp.z_I
z_T         = sbp.z_T

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
                # Grab the scales t, z, and kz
                t  = np.array(f['scales']['sim_time'])
                z  = np.array(f['scales']['z']['1.0'])
                kz = np.array(f['scales']['kz'])
    return t, z, kz, psi_array#, psi_c_array

def filter_neg_freqs(ftd, freq):
    # Find number of frequencies
    nf = len(freq)
    # Filter out negative frequencies
    for j in range(nf):
        if freq[j] < 0.0:
            # Gets rid of negative freq's
            ftd[:,j] = 0.0
        else:
            # Corrects for lost amplitude
            ftd[:,j] = ftd[:,j] * 2.0
    return ftd

def apply_band_pass(ftd, freq, omega, bp_wid=1):
    """
    Applies rectangular band pass to data already FT'd in time. Fragile function

    ftd         data as output by np.fft.fft
    freq        array of frequencies as output by np.fft.fftfreq
    omega       frequency of forced waves
    bp_wid      number extra indices to include on either side of idx_om
    """
    # Find number of frequencies
    nf = len(freq)
    # Find index for band pass, only consider indices of positive freqs
    idx_om = hf.find_nearest_index(freq[0:nf//2], omega)
    # Filter out frequencies left of band pass
    ftd[:,0:idx_om-1-bp_wid]  = 0.0
    # Filter out frequencies left of omega
    ftd[:,idx_om+1+bp_wid:-1] = 0.0
    # Correct for lost amplitude
    ftd[:,idx_om-bp_wid:idx_om+bp_wid] = ftd[:,idx_om-bp_wid:idx_om+bp_wid] * 2.0
    return ftd

# fourier transform in time, band pass around omega, inverse fourier transform
def FT_in_time(t, z, data, dt, omega):
    # FT in time of the data (axis 1 is time)
    ftd = np.fft.fft(data, axis=1)
    # find relevant frequencies
    freq = np.fft.fftfreq(len(t), dt)
    # Apply filtering on frequencies
    use_bp = False
    if use_bp:
        # Apply band pass
        ftd = apply_band_pass(ftd, freq, omega, bp_wid=63)
    else:
        # Filter out just the negative frequencies
        ftd = filter_neg_freqs(ftd, freq)
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
            # for up, remove values for positive wave numbers
            fzdp[i,:] = 0.0
        else:
            # for down, remove values for negative wave numbers
            fzdn[i,:] = 0.0
    # inverse fourier transform in space (z)
    ifzdp = np.fft.ifft(fzdp, axis=0)
    ifzdn = np.fft.ifft(fzdn, axis=0)
    return ifzdp, ifzdn

###############################################################################
# Get the data from the snapshot files
t, z, kz, psi = get_h5_data(tasks, h5_files)
# z and psi arrays come out sorted from most positive to most negative on z axis
#   This flips things around (ud = up / down)
z = np.flip(z)
psi = np.flipud(psi)

BP_array = hf.BP_n_steps(sbp.n_steps, sbp.z, sbp.z0_str, sbp.zf_str)

###############################################################################
# Sanity check plots

if sbp.plot_spacetime:
    hf.plot_z_vs_t(z, t, T, psi.real, BP_array, mL, theta, omega, z0_dis=z0_dis, zf_dis=zf_dis, z_I=z_I, z_T=z_T, plot_full_domain=sbp.plot_full_domain, nT=T_skip, title_str=run_name)

# if sbp.plot_wavespace:
#     hf.plot_k_vs_t(z, t, T, psi.real, psi.imag, k, m, omega, title_str='psi', filename='f_1D_psi_r_i.png')
#
# if sbp.plot_wavespace:
#     hf.plot_k_vs_t(hf.sort_k_coeffs(kz,1024), t, T, psi_c.real, psi_c.imag, k, m, omega, title_str='psi_c', filename='f_1D_psic_r_i.png')

# raise SystemExit(0)

###############################################################################
# Complex demodulation

def Complex_Demodulate(t_then_z, t, z, kz, data, dt, omega):
    if t_then_z == True:
        ## Step 1
        ift_t_y = FT_in_time(t, z, data, dt, omega)
        ## Step 2
        ift_z_y_p, ift_z_y_n = FT_in_space(t, kz, ift_t_y)
        # Get up and down fields as F = |mag_f| * exp(i*phi_f)
        up_field = ift_z_y_p.real * np.exp(np.real(1j * ift_z_y_p.imag))
        dn_field = ift_z_y_n.real * np.exp(np.real(1j * ift_z_y_n.imag))
    else:
        ## Step 1
        ift_z_y_p, ift_z_y_n = FT_in_space(t, kz, data)
        ## Step 2
        up_f = FT_in_time(t, z, ift_z_y_p, dt, omega)
        dn_f = FT_in_time(t, z, ift_z_y_n, dt, omega)
        # Get up and down fields as F = |mag_f| * exp(i*phi_f)
        up_field = up_f.real * np.exp(np.real(1j * up_f.imag))
        dn_field = dn_f.real * np.exp(np.real(1j * dn_f.imag))
    return up_field, dn_field

t_then_z = True
up_field, dn_field = Complex_Demodulate(t_then_z, t, z, kz, psi, dt, omega)

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

big_T = hf.measure_T(dn_field, z, z_I, z_T, T_skip=T_skip, T=T, t=t)
print("Transmission coefficient is:", big_T)

###############################################################################
# More plotting for up and down waves

if sbp.plot_amplitude:
    hf.plot_A_of_I_T(z, t, T, dn_field, z_I, z_T, dz, mL, theta, omega, nT=T_skip, title_str=run_name)

if sbp.plot_amplitude:
    hf.plot_AA_for_z(BP_array, dn_field, z, mL, theta, omega, T_skip=T_skip, T=T, t=t, z0_dis=z0_dis, zf_dis=zf_dis, z_I=z_I, z_T=z_T, title_str=run_name, filename='f_1D_AA_for_z.png')

if sbp.plot_up_dn:
    hf.plot_z_vs_t(z, t, T, up_field.real, BP_array, mL, theta, omega, z0_dis=z0_dis, zf_dis=zf_dis, z_I=z_I, z_T=z_T, plot_full_domain=plot_f_d, nT=T_skip, title_str=run_name, filename='f_1D_up_field.png')
    hf.plot_z_vs_t(z, t, T, dn_field.real, BP_array, mL, theta, omega, z0_dis=z0_dis, zf_dis=zf_dis, z_I=z_I, z_T=z_T, plot_full_domain=plot_f_d, nT=T_skip, title_str=run_name, filename='f_1D_dn_field.png')

plot_CD_checks = False
if plot_CD_checks:
    # Add up and down fields to see if they reproduce the original psi field
    up_plus_dn = up_field.real + dn_field.real
    hf.plot_z_vs_t(z, t, T, up_plus_dn, BP_array, mL, theta, omega, z0_dis=z0_dis, zf_dis=zf_dis, z_I=z_I, z_T=z_T, plot_full_domain=plot_f_d, nT=T_skip, title_str=run_name, filename='f_1D_up_plus_dn.png')
    # Plot the difference, which ideally should be zero everywhere
    CD_diff = psi.real - up_plus_dn
    hf.plot_z_vs_t(z, t, T, CD_diff, BP_array, mL, theta, omega, z0_dis=z0_dis, zf_dis=zf_dis, z_I=z_I, z_T=z_T, plot_full_domain=plot_f_d, nT=T_skip, title_str=run_name, filename='f_1D_CD_diff.png')
