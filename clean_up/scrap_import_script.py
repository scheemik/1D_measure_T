import h5py
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('dark_background')
import colorcet as cc
from dedalus import public as de

# Import SwitchBoard Parameters (sbp)
import switchboard as sbp
import helper_functions as hf
import helper_functions_CD as hfCD

plot_untrimmed = True

def Complex_Demodulation(ft, kz, psi_up):
    # Filter in time
    n_ft = ft.shape[0] # ft is a 2D array with 1 column
    # Add 1 because Nyquist is dropped
    if hf.is_power_of_2(n_ft+1) == False:
        print('Warning: n_ft is not a power of 2:',n_ft)
        print('  Are you plotting untrimmed data?')
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

# Inspired by this post from Jeff https://groups.google.com/g/dedalus-users/c/mTZliZszO_M/m/-_I76EJDAQAJ

# In order to plot vs. time, I needed to fudge a `t_basis` which only works because
#   I force my simulations to have constant dt, and the Fourier basis is evenly spaced

def import_h5_data(zf_dis, z0_dis, nt_keep, dtype=sbp.grid_data_type, dealias_f=sbp.dealias):
    with h5py.File("consolidated_analysis1.h5", mode='r') as f:
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
    return psi_tr, domain_tr, psi, domain
    # return psi_tr, t_tr, z_tr, kz_tr, psi, t, z, kz

psi_tr, domain_tr, psi, domain = import_h5_data(sbp.zf_dis, sbp.z0_dis, sbp.nt_snap)

def get_domain_scales(domain):
    # For grid space
    t = domain.grid(0)[:,0]
    z = domain.grid(1)[0,:]
    # For coefficient space
    f = domain.elements(0)
    k = domain.elements(1)
    return t, z, f, k

def plot_some_stuff(psi, domain):
    t, z, f, k = get_domain_scales(domain)
    psi_up, psi_dn = Complex_Demodulation(f, k, psi)
    # plotting
    from dedalus.extras.plot_tools import plot_bot_2d
    figkw = {'figsize':(6,4), 'dpi':100}
    # Plot log magnitude of spectral coefficients
    psi_dn['c']
    # psi_tr.towards_coeff_space()
    log_mag = lambda xmesh, ymesh, data: (xmesh, ymesh, np.log10(np.abs(data)))
    plot_bot_2d(psi_dn, func=log_mag, clim=(-20, 0), cmap='viridis', title="log10(abs(f['c'])", figkw=figkw);
    plt.show()

    plot_psi_dn = np.flipud(np.transpose(psi_dn['g'].real))
    plt.pcolormesh(t[:], np.flip(z[:]), plot_psi_dn, cmap=cc.cm.bkr);
    plt.show()

    plot_psi_up = np.flipud(np.transpose(psi_up['g'].real))
    plt.pcolormesh(t[:], np.flip(z[:]), plot_psi_up, cmap=cc.cm.bkr);
    plt.show()

if plot_untrimmed:
    plot_some_stuff(psi, domain)
else:
    plot_some_stuff(psi_tr, domain_tr)
