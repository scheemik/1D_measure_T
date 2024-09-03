"""

This script plots the simulated and analytical values of the transmission coefficient
The range in kL for the 2 layer case has been restricted to [0,2.5] from the 1 layer
    case of [0,5] in order for the plot to match the data shown in figure 3 a and b
    of Ghaemsaidi et al 2016

Usage:
    plot_1and2l_exp_data.py EXP SIMS

Options:
    EXP             # name of the experiment from -e
    SIMS            # number of simulations from -s

"""

import numpy as np
import matplotlib.pyplot as plt
import csv

# Parse input parameters
# from docopt import docopt
# args = docopt(__doc__)
exp_name    = 'repro_Ghaemsaridi_fig_3ab_MIT_alpha'#args['EXP']       # Experiment name, needed to find all csv files
num_sims    = 64#int(args['SIMS']) # Number of simulations in the experiment

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
# Set plot parameters
csv_files = ['1L_exp_data.csv',
             '2L_exp_data.csv',
             '1L_MIT_data.csv',
             '2L_MIT_data.csv']
plt_clrs = [hf.CSS4_COLORS['forestgreen'],
            hf.CSS4_COLORS['blueviolet'],
            # hf.my_clrs['warm-salty'], 
            # hf.my_clrs['cold-fresh']]
            'b', 
            'r']
labels   = [r'$\mathbb{T}_1$',
            r'$\mathbb{T}_2$',
            r'$\mathbb{T}_1$ G16',
            r'$\mathbb{T}_2$ G16']
markers  = ['o',
            'x',
            'o',
            'x']
linestyles=['-',
            '-',
            ':',
            ':']

def combine_exp_csv_files(csv_files, output_file, num_sims):
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
# all_exp_data = combine_exp_csv_files(csv_files, "all_exp_data.csv", 128)
# all_MIT_data = combine_exp_csv_files(MIT_files, "all_MIT_data.csv", 64)

def plot_T_vs_kL(csv_files, theta, omega, title_str='Transmission Coefficient', filename='pres_64s_1and2l_T_vs_kL.png', d_mode=False):
    """
    Plots the transmission coefficient as a function of kL

    csv_files       List of csv files with the values of kL and transmission coefficients
    theta           Angle at which wave is incident on stratification structure
    omega           frequency of wave
    """
    # Set figure and axes for plot, fig_ratio=1 makes a square
    fig, axes = hf.set_fig_axes([1], [1], fig_ratio=0.6)
    # Create array of analytical function
    kL_ana = np.linspace(0, 5, 100)
    ana_data = hf.SY_eq2_4(theta, kL_ana)
    # Plot line for analytical solution
    if d_mode:
        ana_clr = 'w'
    else:
        ana_clr = 'k'
    axes.plot(kL_ana, ana_data, color=ana_clr, label=r'$\mathbb{T}_{ana}$', zorder=5)
    # Plot points from MIT MATLAB code
    for i in range(len(csv_files)):
        # Read in exp csv row by row
        with open(csv_files[i], 'r') as datafile:
            csvreader = csv.reader(datafile)
            data = list(csvreader)
        # Format lists into a numpy array
        data_arr = np.array(data)
        # print(data_arr[:,1].astype('float64'))
        # exit(0)
        axes.plot(data_arr[:,1].astype('float64'), data_arr[:,3].astype('float64'), color=plt_clrs[i], marker=markers[i], linestyle=linestyles[i], label=labels[i], alpha=0.7)
    axes.legend()
    # Set x axis limits
    axes.set_xlim([0, 2.5])
    # Add labels and titles
    axes.set_xlabel(r'$kL$')
    axes.set_ylabel(r'$\mathbb{T}$')
    plt.title(r'%s for $\theta=$%s' %(title_str, hf.rad_to_degs(theta)))
    #plt.show()
    plt.savefig(filename, dpi=300, transparent=True)

file_name = exp_name + '.png'
plot_T_vs_kL(csv_files, theta, omega, filename=file_name)

# dark mode
file_name_d = exp_name + '_d.png'
plt.style.use('dark_background')
plot_T_vs_kL(csv_files, theta, omega, filename=file_name_d, d_mode=True)
