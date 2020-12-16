"""

This script plots the simulated and analytical values of the transmission coefficient

Usage:
    plot_exp_data.py EXP SIMS

Options:
    EXP             # name of the experiment from -e
    SIMS            # number of simulations in the experiment from -s

"""

import numpy as np
import matplotlib.pyplot as plt

# Parse input parameters
from docopt import docopt
args = docopt(__doc__)
exp_name    = args['EXP']       # Experiment name, for finding csv files
num_sims    = args['SIMS']      # Number of simulations

###############################################################################
# Import SwitchBoard Parameters (sbp)
#   This import assumes the default simulation naming convention
import importlib
import switchboard as sbp
# switchboard_module = "000_" + exp_name + "." + switchboard
# sbp = importlib.import_module(switchboard_module)
# Add functions in helper file
import helper_functions as hf

# Physical parameters
# nu          = sbp.nu            # [m^2/s] Viscosity (momentum diffusivity)
# kappa       = sbp.kappa         # [m^2/s] Thermal diffusivity
# f_0         = sbp.f_0           # [s^-1]        Reference Coriolis parameter
# g           = sbp.g             # [m/s^2] Acceleration due to gravity
# Problem parameters
theta       = sbp.theta         # [rad]         Propagation angle from vertical
omega       = sbp.omega         # [rad s^-1]    Wave frequency

###############################################################################
# Read data from the csv files

import csv
data_arr = []
for i in range(num_sims):
    sim_csv = str(id) + "_" + str(exp_name) + "/sim_data.csv"
    with open(csv_file, 'r') as datafile:
        csvreader = csv.reader(datafile)
        np.append(data_arr, np.array(list(csvreader)), axis=1)

# csv read in as U32, convert data to float 64 before using
mL_array = data_arr[:,1].astype('float64')
sim_data = data_arr[:,3].astype('float64')

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
    axes.scatter(mL_array, sim_data, color='red', marker='o', label=r'$\mathbb{T}_{sim}$')
    axes.legend()
    # Add labels and titles
    axes.set_xlabel(r'$mL$')
    axes.set_ylabel(r'$\mathbb{T}$')
    plt.title(r'%s for $\theta=$%s' %(title_str, hf.rad_to_degs(theta)))
    #plt.show()
    plt.savefig(filename)

plot_T_vs_mL(mL_array, sim_data, theta, omega)
