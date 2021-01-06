"""

This script plots the simulated and analytical values of the transmission coefficient

Usage:
    plot_exp_data.py EXP SIMS

Options:
    EXP             # name of the experiment from -e
    SIMS            # number of simulations from -s

"""

import numpy as np
import matplotlib.pyplot as plt

# Parse input parameters
from docopt import docopt
args = docopt(__doc__)
exp_name    = args['EXP']       # Experiment name, needed to find all csv files
num_sims    = int(args['SIMS']) # Number of simulations in the experiment

###############################################################################
# Import relevant modules

# Import SwitchBoard Parameters (sbp)
#   This import assumes the default simulation naming convention
import importlib
helper_module = "_simulations.000_" + exp_name + ".helper_functions"
hf = importlib.import_module(helper_module)

# Problem parameters -> eventually should be taken from switchboard
# theta       = sbp.theta         # [rad]         Propagation angle from vertical
# omega       = sbp.omega         # [rad s^-1]    Wave frequency
theta   = np.arctan(1)
omega   = 1 * np.cos(theta)

###############################################################################
# Read data from the csv file

import csv
# Make blank list to hold csv data
data = [None]*num_sims
# Read in exp csv row by row
with open(csv_file, 'r') as datafile:
    csvreader = csv.reader(datafile)
    data = list(csvreader)[0]
# Format lists into a numpy array
data_arr = np.array(data)

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
