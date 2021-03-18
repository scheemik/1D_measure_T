import h5py
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('dark_background')
import colorcet as cc
from dedalus import public as de
from dedalus.core import operators as ops
from dedalus.extras.plot_tools import quad_mesh, pad_limits

# Import SwitchBoard Parameters (sbp)
import switchboard as sbp
import helper_functions as hf
import helper_functions_CD as hfCD

plot_untrimmed = False

def trim_data(z, t, psi_data, z0_dis, zf_dis, nt_keep):
    # Transpose psi and flip z to allow trimming (see find_nearest_index)
    z = np.flip(np.array(z))
    psi = np.transpose(np.array(psi_data[()]))
    # Trim data in space
    z_tr, psi_tr_data = hf.trim_data_z(z, psi, z0_dis, zf_dis)
    # Trim data in time
    t_tr, psi_tr_data = hf.trim_data_t(t, psi_tr_data, nt_keep)
    # Transpose and flip arrays back to make domain
    psi_tr_data = np.transpose(psi_tr_data)
    z_tr = np.flip(z_tr)
    return z_tr, t_tr, psi, psi_tr_data

def Complex_Demodulation(ft, kz, psi_up):
    # Filter in time
    # Remove all negative frequencies
    selection1 = (ft<0) * (kz==kz)
    psi_up['c'][selection1] = 0
    # Multiply remaining by 2 to compensate for lost magnitude
    double_op = 2*psi_up
    double_op.operate(psi_up)

    # Filter in space
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
        z_tr, t_tr, psi, psi_tr_data = trim_data(z, t, psi_data, z0_dis, zf_dis, nt_keep)

        # Bases and domain
        t_tr_basis = de.Fourier('t', len(t_tr), interval=(t_tr[0], t_tr[-1]), dealias=dealias_f)
        z_tr_basis = de.Fourier('z', len(z_tr), interval=(z_tr[0], z_tr[-1]), dealias=dealias_f)
        domain_tr = de.Domain([t_tr_basis,z_tr_basis], grid_dtype=dtype)
        # set field
        psi_tr = domain_tr.new_field(name = 'psi_tr')
        psi_tr['g'] = psi_tr_data
        t_tr = domain_tr.grid(0)[:,0]
        z_tr = domain_tr.grid(1)[0,:]
        kz_tr = z_tr_basis.wavenumbers
        ft_el = domain_tr.elements(0)
        kz_el = domain_tr.elements(1)

        psi_up, psi_dn = Complex_Demodulation(ft_el, kz_el, psi_tr)

        from dedalus.extras.plot_tools import plot_bot_2d
        figkw = {'figsize':(6,4), 'dpi':100}
        # Plot log magnitude of spectral coefficients
        psi_dn['c']
        # psi_tr.towards_coeff_space()
        log_mag = lambda xmesh, ymesh, data: (xmesh, ymesh, np.log10(np.abs(data)))
        plot_bot_2d(psi_dn, func=log_mag, clim=(-20, 0), cmap='viridis', title="log10(abs(f['c'])", figkw=figkw);
        plt.show()

        # psi_tr['g']
        # psi_tr.towards_grid_space()
        plot_psi_dn = np.flipud(np.transpose(psi_dn['g'].real))
        plt.pcolormesh(t_tr[:], np.flip(z_tr[:]), plot_psi_dn, cmap=cc.cm.bkr);
        plt.show()

        plot_psi_up = np.flipud(np.transpose(psi_up['g'].real))
        plt.pcolormesh(t_tr[:], np.flip(z_tr[:]), plot_psi_up, cmap=cc.cm.bkr);
        plt.show()

        if plot_untrimmed:
            # Bases and domain
            t_basis = de.Fourier('t', len(t), interval=(t[0], t[-1]), dealias=dealias_f)
            z_basis = de.Fourier('z', len(z), interval=(z[0], z[-1]), dealias=dealias_f)
            domain = de.Domain([t_basis,z_basis], grid_dtype=dtype)
            # set field
            psi = domain.new_field(name = 'psi')
            psi['g'] = psi_data
            t  = np.array(t)
            kz = z_basis.wavenumbers
        else:
            kz  = None
    return psi_tr, t_tr, z_tr, kz_tr, psi, t, z, kz

psi_tr, t_tr, z_tr, kz_tr, psi, t, z, kz = import_h5_data(sbp.zf_dis, sbp.z0_dis, sbp.nt_snap)

# # Plot data
# plt.figure(figsize=(7,4), dpi=100)
# plt.pcolormesh(t[:], z[:], np.transpose(psi['g'].real), shading='nearest', cmap=cc.cm.bkr)
# plt.colorbar()
# plt.xlabel('t')
# plt.ylabel('z')
# plt.title('Test plot of psi')
# plt.tight_layout()
# plt.show()
