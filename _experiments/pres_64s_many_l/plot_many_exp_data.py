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
csv_files = ["1l_exp_data.csv", "2l_exp_data.csv", "3l_exp_data.csv", "4l_exp_data.csv"]

def combine_exp_csv_files(csv_files, output_file):
    import csv
    # Find the number of experiments' data to combine
    n_exp = len(csv_files)
    # Make blank array to hold all the data
    #   n_cols = n_exp+1 to hold column of kL values
    combined_data_arr = np.array([[None]*(n_exp+1)]*num_sims)
    # print("shape of combined data array = ",combined_data_arr.shape)
    i = 0
    for csv_f in csv_files:
        # Make blank list to hold csv data
        data = [None]*num_sims
        # Read in exp csv row by row
        with open(csv_f, 'r') as datafile:
            csvreader = csv.reader(datafile)
            data = list(csvreader)
        # Format lists into a numpy array
        data_arr = np.array(data)

        # csv read in as U32, convert data to float 64 before using
        if i == 0:
            combined_data_arr[:,0] = data_arr[:,1].astype('float64')
        combined_data_arr[:,i+1] = data_arr[:,3].astype('float64')
        i += 1
    # Write data to new output csv file
    np.savetxt(output_file, combined_data_arr, delimiter=",")
    return combined_data_arr
all_exp_data = combine_exp_csv_files(csv_files, "all_exp_data.csv")

def plot_T_vs_kL(all_exp_data, theta, omega, title_str='Transmission Coefficient', filename='pres_64s_T_vs_kL.png'):
    """
    Plots the transmission coefficient as a function of kL

    all_exp_data    2D array, 1st column kL_arr, all other columns are experiments
    theta           Angle at which wave is incident on stratification structure
    omega           frequency of wave
    """
    # Set figure and axes for plot, fig_ratio=1 makes a square
    fig, axes = hf.set_fig_axes([1], [1], fig_ratio=0.8)
    # Create array of analytical function
    kL_ana = np.linspace(all_exp_data[0,0], all_exp_data[-1,0], 100)
    ana_data = hf.SY_eq2_4(theta, kL_ana)
    # Plot line for analytical solution
    axes.plot(kL_ana, ana_data, color='black', label=r'$\mathbb{T}_{ana}$')
    # Get colors
    plt_clrs = [hf.CSS4_COLORS['forestgreen'],
                hf.CSS4_COLORS['blueviolet'],
                hf.CSS4_COLORS['tomato'],
                hf.CSS4_COLORS['royalblue']]
    # Plot points from measurements in simulations
    for i in range(1,all_exp_data.shape[1]):
        axes.plot(all_exp_data[:,0], all_exp_data[:,i], color=plt_clrs[i-1], marker='.', label=rf'$T_{i}$')
    axes.legend()
    # Add labels and titles
    axes.set_xlabel(r'$kL$')
    axes.set_xlim(0,2.5)
    axes.set_ylabel(r'$\mathbb{T}$')
    plt.title(r'%s for $\theta=$%s' %(title_str, hf.rad_to_degs(theta)))
    #plt.show()
    plt.savefig(filename)

plot_T_vs_kL(all_exp_data, theta, omega)
