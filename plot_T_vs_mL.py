"""

This script plots:
    Simulated and analytical values of the transmission coefficient

"""

import numpy as np
import matplotlib.pyplot as plt

###############################################################################
# Import SwitchBoard Parameters (sbp)
#   This import assumes the switchboard is in the same directory as the core code
import switchboard as sbp
# Add functions in helper file
import helper_functions as hf

# Physical parameters
nu          = sbp.nu            # [m^2/s] Viscosity (momentum diffusivity)
#kappa       = sbp.kappa         # [m^2/s] Thermal diffusivity
f_0         = sbp.f_0           # [s^-1]        Reference Coriolis parameter
g           = sbp.g             # [m/s^2] Acceleration due to gravity
# Problem parameters
theta       = sbp.theta         # [rad]         Propagation angle from vertical
omega       = sbp.omega         # [rad s^-1]    Wave frequency

###############################################################################
# Get depth and wavenumber axes
z           = sbp.z
ks          = sbp.ks

mL_array = [0, 1, 2, 3, 4, 5]
sim_data = [1.0087048549269324,
            0.4627298453648014,
            0.10604904571929356,
            0.024508789111336624,
            0.008528472600904731,
            0.00428457274527059]

def plot_T_vs_mL(mL_array, sim_data, theta, omega, title_str='Transmission Coefficient', filename='f_1D_T_vs_mL.png'):
    """
    Plots the amplitude of the downward propagating wave as a function of time
        for both the incident and transmission depths

    mL_array    1D array of mL values used in simulations
    sim_data    1D array of measured transmission coefficients for values in mL_array
    theta       Angle at which wave is incident on stratification structure
    omega       frequency of wave
    """
    # Set figure and axes for plot
    fig, axes = hf.set_fig_axes([1], [1])
    # Create array of analytical function
    mL_ana = np.linspace(mL_array[0], mL_array[-1], 100)
    ana_data = hf.SY_eq2_4(theta, mL_ana)
    # Plot line for analytical solution
    axes.plot(mL_ana, ana_data, color='black', label=r'$\mathbb{T}_{ana}$')
    # Plot points from measurements in simulations
    axes.plot(mL_array, sim_data, color='red', marker='o', label=r'$\mathbb{T}_{sim}$')
    axes.legend()
    # Add labels and titles
    axes.set_xlabel(r'$mL$')
    axes.set_ylabel(r'$\mathbb{T}$')
    plt.title(r'%s for $\theta=$%s' %(title_str, hf.rad_to_degs(theta)))
    #plt.show()
    plt.savefig(filename)

plot_T_vs_mL(mL_array, sim_data, theta, omega)
