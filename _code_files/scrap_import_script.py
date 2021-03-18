import h5py
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('dark_background')
import colorcet as cc
from dedalus import public as de
from dedalus.extras.plot_tools import quad_mesh, pad_limits

# Import SwitchBoard Parameters (sbp)
import switchboard as sbp
import helper_functions as hf
import helper_functions_CD as hfCD

plot_untrimmed = False

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
        # print(kz_tr)
        # print(kz_el)
        # print("Result : ",np.array_equal(kz_tr,kz_el[0]))
        # Write out results to file
        # import csv
        # csv_file = "wavenumbers.csv"
        # with open(csv_file, 'a') as datafile:
        #     csvwriter = csv.writer(datafile)
        #     csvwriter.writerow(kz_tr)
        #     csvwriter.writerow(kz_el[0])

        #CD
        # Filter in time
        # Remove all negative frequencies
        selection1 = (ft_el<0) * (kz_el==kz_el)
        psi_tr['c'][selection1] = 0
        # Multiply remaining by 2 to compensate for lost magnitude
        double_op = 2*psi_tr
        double_op.operate(psi_tr)

        # Filter in space
        # Take the up (down?) propagating waves
        selection2 = (ft_el==ft_el) * (kz_el>0)
        psi_tr['c'][selection2] = 0

        from dedalus.extras.plot_tools import plot_bot_2d
        figkw = {'figsize':(6,4), 'dpi':100}
        # Plot log magnitude of spectral coefficients
        psi_tr['c']
        # psi_tr.towards_coeff_space()
        log_mag = lambda xmesh, ymesh, data: (xmesh, ymesh, np.log10(np.abs(data)))
        plot_bot_2d(psi_tr, func=log_mag, clim=(-20, 0), cmap='viridis', title="log10(abs(f['c'])", figkw=figkw);
        plt.show()

        # psi_tr['g']
        # psi_tr.towards_grid_space()
        plot_psi = np.flipud(np.transpose(psi_tr['g'].real))
        print(plot_psi[0][0])
        plt.pcolormesh(t_tr[:], np.flip(z_tr[:]), plot_psi, cmap=cc.cm.bkr);
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
