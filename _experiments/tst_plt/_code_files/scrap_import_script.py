import h5py
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('dark_background')
import colorcet as cc
from dedalus import public as de
from dedalus.extras.plot_tools import quad_mesh, pad_limits

# Import SwitchBoard Parameters (sbp)
import switchboard as sbp


# Inspired by this post from Jeff https://groups.google.com/g/dedalus-users/c/mTZliZszO_M/m/-_I76EJDAQAJ

# In order to plot vs. time, I needed to fudge a `t_basis` which only works because
#   I force my simulations to have constant dt, and the Fourier basis is evenly spaced

with h5py.File("consolidated_analysis1.h5", mode='r') as file:
    # Get task
    psi_data = file['tasks']['psi']
    # Get dimensions
    t = psi_data.dims[0]['sim_time']
    z = psi_data.dims[1][0]
    # Bases and domain
    z_basis = de.Fourier('z', len(z), interval=(z[0], z[-1]), dealias=sbp.dealias)
    t_basis = de.Fourier('t', len(t), interval=(t[0], t[-1]), dealias=sbp.dealias)
    domain = de.Domain([t_basis,z_basis], grid_dtype=np.complex128)
    # set field
    psi = domain.new_field(name = 'psi')
    psi['g'] = psi_data
    # Plot data
    plt.figure(figsize=(7,4), dpi=100)
    plt.pcolormesh(t[:], z[:], np.transpose(psi['g'].real), shading='nearest', cmap=cc.cm.bkr)
    plt.colorbar()
    plt.xlabel('t')
    plt.ylabel('z')
    plt.title('Test plot of psi')
    plt.tight_layout()
    plt.show()
