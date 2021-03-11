import h5py
import numpy as np
import matplotlib.pyplot as plt
import colorcet as cc
from dedalus import public as de
from dedalus.extras.plot_tools import quad_mesh, pad_limits

# Import SwitchBoard Parameters (sbp)
import switchboard as sbp

# domain = sbp.domain
# slices = domain.dist.grid_layout.slices(scale=(1,1))

# see https://groups.google.com/g/dedalus-users/c/TWmorKHizZ4/m/J-kDsyKJBAAJ
def limits_from_grid(basis_type,grid):
    if basis_type == 'Fourier':
        # even grid spacing
        delta = grid[1] - grid[0]
        left_edge = grid[0]
        right_edge = grid[-1] + delta
    elif basis_type == 'SinCos':
        delta = grid[1] - grid[0]
        left_edge = grid[0] - delta/2.
        right_edge = grid[-1] + delta/2.
    elif basis_type == 'Chebyshev':
        # Use Gauss points
        raise NotImplementedError("Chebyshev limits not yet supported")
    else:
        raise ValueError("Unknown basis type:", basis_type)

    return left_edge, right_edge

def domain_from_file(filename,basis_types,dealias=3/2):
    with h5py.File(filename,'r') as dfile:
        testkey = next(iter(dfile['tasks']))
        testdata = dfile['tasks'][testkey]
        dims = testdata.shape[1:] # trim write dimension
        dim_names = [i.decode('utf-8') for i in testdata.attrs['DIMENSION_LABELS'][1:]]
        if not testdata.attrs['grid_space'][-1]:
            dims[-1] *= 2
        bases = []
        for n,b,d in zip(dim_names,basis_types,dims):
            if b == 'Fourier':
                limits = limits_from_grid('Fourier',dfile['/scales/'+n+'/1.0'])
                basis = de.Fourier(n, d, interval=limits, dealias=dealias)
            elif b == 'SinCos':
                limits = limits_from_grid('SinCos',dfile['/scales/'+n+'/1.0'])
                basis = de.SinCos(n, d, interval=limits, dealias=dealias)
            else:
                raise ValueError("Unknown basis type:", basis_type)
            bases.append(basis)
        d = de.Domain(bases,grid_dtype=np.float64) # hardcoded domain dtype for now
        return d

def fields_from_file(filename,basis_types,field,meta=None,index=-1):
    domain = domain_from_file(filename, basis_types)
    with h5py.File(filename,'r') as dfile:
        f = domain.new_field()
        dset = dfile['tasks'][field]
        if meta:
            for k,v in meta.items():
                f.meta[k] = v
        for layout in domain.dist.layouts:
            if np.allclose(layout.grid_space, dset.attrs['grid_space']):
                break
        else:
            raise ValueError("No matching layout")
        f[layout] = dset[(index,) + (slice(None))]
    return f

# psi = fields_from_file("consolidated_analysis1.h5", ['Fourier'], 'psi')
# print(psi)

# This one is by Jeff https://groups.google.com/g/dedalus-users/c/mTZliZszO_M/m/-_I76EJDAQAJ

with h5py.File("consolidated_analysis1.h5", mode='r') as file:
    psi_data = file['tasks']['psi']
    t = psi_data.dims[0]['sim_time']
    z = psi_data.dims[1][0]
    # Get bases
    # t  = np.array(file['scales']['sim_time'])
    print("len(t) = ",len(t))
    # z  = np.flip(np.array(file['scales']['z']['1.0']))
    print("len(z) = ",len(z))

    # Get domain
    # Bases and domain
    z_basis = de.Fourier('z', len(z), interval=(z[0], z[-1]), dealias=sbp.dealias)
    t_basis = de.Fourier('t', len(t), interval=(t[0], t[-1]), dealias=sbp.dealias)
    domain = de.Domain([t_basis,z_basis], grid_dtype=np.complex128)

    psi = domain.new_field(name = 'psi')
    # # Load datasets
    psi['g'] = psi_data
    # # psi['g'] = file.get('psi')[slices]
    print(psi['g'].shape)

    # Make meshgrid
    xmesh, ymesh = quad_mesh(x=t, y=z)
    # Plot data
    plt.figure(figsize=(6,7), dpi=100)
    # plt.imshow(psi['g'].real, cmap='RdBu_r')
    plt.pcolormesh(t[:], z[:], np.transpose(psi['g'].real), shading='nearest', cmap=cc.cm.bkr)
    # plt.colorbar()
    # plt.xlabel('t')
    # plt.ylabel('z')
    # plt.title('Test plot of psi')
    # plt.tight_layout()
    plt.show()
