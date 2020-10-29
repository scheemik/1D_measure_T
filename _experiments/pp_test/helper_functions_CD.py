# Helper functions for Dedalus experiment
"""
Description:
This contains helper functions for the Dedalus code so the same version of
    functions can be called by multiple scripts
This script contains just the functions related to Complex Demodulation
"""

import numpy as np
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as clrs
from matplotlib import ticker
from dedalus.extras.plot_tools import quad_mesh, pad_limits

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

###############################################################################
# Get data in spectral form

def z_t_to_k_omega(t, z, data, dt, dz):
    """
    Takes in data (z,t) and returns spectral form (k, omega)
    """
    # FT in time of the data (axis 1 is time)
    ft_t      = np.fft.fft(data, axis=1)
    # find relevant frequencies
    freqs     = np.fft.fftfreq(len(t), dt)
    # FT in space (z) of the data (axis 0 is z)
    spec_data = np.fft.fft(ft_t, axis=0)
    # find relevant wavenumbers (k)
    ks        = np.fft.fftfreq(len(z), dz)
    return freqs, ks, spec_data
