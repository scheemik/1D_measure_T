"""
Runs an idealized steady state test of various functions.
Run with $ python3 steady_state_test.py

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
z0_dis      = sbp.z0_dis
zf_dis      = sbp.zf_dis
lz_lambda = abs(z0_dis-zf_dis)/lam_z
print("Lz/lambda=",lz_lambda)

dt          = sbp.dt
dz          = abs(z0_dis - zf_dis)/nz

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
z = np.linspace(zf_dis, z0_dis, nz) # careful about the order of endpoints
#stop_sim_time, nt  = hf.extended_stop_time(sbp.sim_time_stop, dt)
stop_sim_time = sbp.sim_time_stop
nt = stop_sim_time / dt
t = np.linspace(0.0, stop_sim_time, int(nt))
# Find wavenumbers
kz = np.fft.fftfreq(len(z), dz)
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
# Additional post-processing helper functions

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

def is_power_of_2(x):
    n = np.log2(x)
    return (np.floor(n) == np.ceil(n))

# fourier transform in time, band pass around omega, inverse fourier transform
def FT_in_time(t, z, data, dt, omega):
    # Check to make sure the time dimension matches
    nt = len(t)
    print("nt=",nt)
    nt_data = data.shape[1]
    if nt != nt_data:
        print('MISMATCH IN nt!')
    if is_power_of_2(nt) == False:
        print('nt is not a power of 2')
    # FT in time of the data (axis 1 is time)
    ftd = np.fft.fft(data, axis=1)
    # find relevant frequencies
    freq = np.fft.fftfreq(nt, dt)
    # Apply filtering on frequencies
    use_bp = False
    if use_bp:
        # Apply band pass
        ftd = apply_band_pass(ftd, freq, omega, bp_wid=1)
    else:
        # Filter out just the negative frequencies
        ftd = filter_neg_freqs(ftd, freq)
    # inverse fourier transform in time of the data
    iftd = np.fft.ifft(ftd, axis=1)
    #   a complex valued signal where iftd.real == data, or close enough
    return iftd, ftd, freq

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
# z and psi arrays come out sorted from most positive to most negative on z axis
#   This flips things around (ud = up / down)
# z = np.flip(z)
# psi = np.flipud(psi)

BP_array = hf.BP_n_layers(0, sbp.z, sbp.z0_str, sbp.zf_str)

###############################################################################
# Complex demodulation

def Complex_Demodulate(t_then_z, t, z, kz, data, dt, omega):
    if t_then_z == True:
        ## Step 1
        ift_t_y = FT_in_time(t, z, data, dt, omega)[0]
        ## Step 2
        up_field, dn_field = FT_in_space(t, kz, ift_t_y)
    else:
        ## Step 1
        ift_z_y_p, ift_z_y_n = FT_in_space(t, kz, data)
        ## Step 2
        up_field = FT_in_time(t, z, ift_z_y_p, dt, omega)[0]
        dn_field = FT_in_time(t, z, ift_z_y_n, dt, omega)[0]
    return up_field, dn_field

# Trimming the data
tr_z, tr_t, tr_psi = hf.trim_data(z, t, psi, z0_dis=None, zf_dis=None, T_cutoff=T_cutoff, T=T)

t_then_z = True
up_field, dn_field = Complex_Demodulate(t_then_z, t, z, kz, psi, dt, omega)
tr_up_field, tr_dn_field = Complex_Demodulate(t_then_z, tr_t, tr_z, kz, tr_psi, dt, omega)

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

I_, T_, AAcc_I, AAcc_T = hf.measure_T(dn_field, z, z_I, z_T, T_skip=None, T=T, t=t)
big_T = T_/I_
print("Transmission coefficient is:", big_T)

###############################################################################
# What to plot
show_trimmed = True
if show_trimmed:
    plot_up = tr_up_field
    plot_dn = tr_dn_field
    plot_psi = tr_psi
    plot_z = tr_z
    plot_t = tr_t
else:
    plot_up     = up_field
    plot_dn     = dn_field
    plot_psi    = psi
    plot_z = z
    plot_t = t

###############################################################################
# Plotting and stuff

if sbp.plot_windows:
    hf.plot_v_profiles(plot_z, BP_array, sbp.win_bf_array, sbp.win_sp_array, mL=mL, theta=theta, omega=omega, z_I=z_I, z_T=z_T, z0_dis=z0_dis, zf_dis=zf_dis, title_str=run_name, filename='ss_1D_windows.png')

if sbp.plot_spacetime:
    hf.plot_z_vs_t(plot_z, plot_t, T, plot_psi.real, BP_array, mL, theta, omega, z_I=z_I, z_T=z_T, z0_dis=z0_dis, zf_dis=zf_dis, plot_full_domain=plt_fd, T_cutoff=T_cutoff, title_str=run_name, filename='ss_1D_wave.png')

if sbp.plot_freqspace:
    foobar, psi_FT_t, freqs = FT_in_time(t, z, psi, dt, omega)
    # fig, axes = plt.subplots(nrows=1, ncols=1)
    # axes.plot(freqs, psi_FT_t[200])
    # plt.show()
    hf.plot_freq_space(z, freqs, psi_FT_t.real, psi_FT_t.imag, mL, theta, omega, plot_full_domain=plt_fd, title_str=run_name, filename='ss_1D_freq_spectra.png')

if sbp.plot_amplitude:
    hf.plot_A_of_I_T(plot_z, plot_t, T, plot_dn, mL, theta, omega, z_I, z_T, T_cutoff=T_cutoff, title_str=run_name, filename='ss_1D_A_of_I_T.png')

if sbp.plot_amplitude:
    hf.plot_AA_for_z(plot_z, BP_array, hf.AAcc(plot_dn), mL, theta, omega, z_I=z_I, z_T=z_T, z0_dis=z0_dis, zf_dis=zf_dis, plot_full_domain=plt_fd, title_str=run_name, filename='ss_1D_AA_for_z.png')

if sbp.plot_up_dn:
    hf.plot_z_vs_t(plot_z, plot_t, T, plot_up.real, BP_array, mL, theta, omega, z0_dis=z0_dis, zf_dis=zf_dis,  plot_full_domain=plt_fd, T_cutoff=None, title_str=run_name+' up', filename='ss_1D_up_field.png')
    hf.plot_z_vs_t(plot_z, plot_t, T, plot_dn.real, BP_array, mL, theta, omega, z0_dis=z0_dis, zf_dis=zf_dis, plot_full_domain=plt_fd, T_cutoff=None, title_str=run_name+' dn', filename='ss_1D_dn_field.png')

plot_CD_checks = False
if plot_CD_checks:
    # Add up and down fields to see if they reproduce the original psi field
    up_plus_dn = plot_up.real + plot_dn.real
    hf.plot_z_vs_t(z, t, T, up_plus_dn, BP_array, mL, theta, omega, z0_dis=z0_dis, zf_dis=zf_dis, plot_full_domain=plt_fd, nT=None, title_str=run_name, filename='ss_1D_up_plus_dn.png')
    # Plot the difference, which ideally should be zero everywhere
    CD_diff = psi.real - up_plus_dn
    hf.plot_z_vs_t(z, t, T, CD_diff, BP_array, mL, theta, omega, z0_dis=z0_dis, zf_dis=zf_dis, plot_full_domain=plt_fd, nT=None, title_str=run_name, filename='ss_1D_CD_diff.png')
