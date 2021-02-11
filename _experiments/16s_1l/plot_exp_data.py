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
helper_module = "_code_files.helper_functions"
# helper_module = "_simulations.000_" + exp_name + ".helper_functions"
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
with open("exp_data.csv", 'r') as datafile:
    csvreader = csv.reader(datafile)
    data = list(csvreader)
# Format lists into a numpy array
data_arr = np.array(data)

# csv read in as U32, convert data to float 64 before using
kL_array = data_arr[:,1].astype('float64')
th_array = data_arr[:,2].astype('float64')
sim_data = data_arr[:,3].astype('float64')

def plot_T_vs_kL(kL_array, sim_data, theta, omega, title_str='Transmission Coefficient', filename='f_1D_T_vs_kL.png'):
    """
    Plots the amplitude of the downward propagating wave as a function of time
        for both the incident and transmission depths

    kL_array    1D array of kL values used in simulations
    sim_data    1D array of measured transmission coefficients for values in kL_array
    theta       Angle at which wave is incident on stratification structure
    omega       frequency of wave
    """
    # Set figure and axes for plot
    fig, axes = hf.set_fig_axes([1], [1])
    # Create array of analytical function
    kL_ana = np.linspace(kL_array[0], kL_array[-1], 100)
    ana_data = hf.SY_eq2_4(theta, kL_ana)
    # Plot line for analytical solution
    axes.plot(kL_ana, ana_data, color='black', label=r'$\mathbb{T}_{ana}$')
    # Plot points from measurements in simulations
    axes.scatter(kL_array, sim_data, color='red', marker='o', label=r'$\mathbb{T}_{sim}$')
    axes.legend()
    # Add labels and titles
    axes.set_xlabel(r'$kL$')
    axes.set_ylabel(r'$\mathbb{T}$')
    plt.title(r'%s for $\theta=$%s' %(title_str, hf.rad_to_degs(theta)))
    #plt.show()
    plt.savefig(filename)

plot_T_vs_kL(kL_array, sim_data, theta, omega)
