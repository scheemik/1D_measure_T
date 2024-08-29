"""

This script plots the simulated and analytical values of the transmission coefficient

"""

import numpy as np
import matplotlib.pyplot as plt

# Enable dark mode plotting
# plt.style.use('dark_background')

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
csv_files = ["1l_exp_data.csv", "2l_exp_data.csv", "3l_exp_data.csv", "4l_exp_data.csv", "5l_exp_data.csv"]
num_sims = 128

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

def plot_T_vs_kL(all_exp_data, theta, omega, title_str='Transmission Coefficient', filename='new_128s_T_vs_kL.png'):
    """
    Plots the transmission coefficient as a function of kL

    all_exp_data    2D array, 1st column kL_arr, all other columns are experiments
    theta           Angle at which wave is incident on stratification structure
    omega           frequency of wave
    """
    this_font_size = 16
    # Set figure and axes for plot, fig_ratio=1 makes a square
    fig, axes = hf.set_fig_axes([1], [1], fig_ratio=0.6)
    # Create array of analytical function
    kL_ana = np.linspace(all_exp_data[0,0], all_exp_data[-1,0], 100)
    ana_data = hf.SY_eq2_4(theta, kL_ana)
    # Plot line for analytical solution
    axes.plot(kL_ana, ana_data, color='black', label=r'$\mathbb{T}_{ana}$', zorder=10)
    # Get colors
    plt_clrs = [hf.CSS4_COLORS['forestgreen'],
                hf.CSS4_COLORS['blueviolet'],
                hf.CSS4_COLORS['tomato'],
                hf.CSS4_COLORS['royalblue'],
                hf.CSS4_COLORS['sienna']]
    labels   = [rf'$1$ layer',
                rf'$2$ layers',
                rf'$3$ layers',
                rf'$4$ layers',
                rf'$5$ layers']
    markers  = ['o',
                'x',
                '^',
                's',
                'p']
    linestyles=['-',
                '-',
                '-.',
                '--',
                ':']
    # Plot points from measurements in simulations
    for i in range(1,all_exp_data.shape[1]):
        axes.plot(all_exp_data[:,0], all_exp_data[:,i], color=plt_clrs[i-1], marker=markers[i-1], linestyle=linestyles[i-1], markersize=5, label=labels[i-1])#rf'$T_{i}$')
    axes.legend(fontsize=this_font_size-2)
    # Add labels and titles
    axes.set_xlabel(r'$k_xL$', fontsize=this_font_size)
    axes.set_xlim(0,2.5)
    axes.set_ylabel(r'$\mathbb{T}$', fontsize=this_font_size+6)
    plt.title(r'%s for $\theta=$%s' %(title_str, hf.rad_to_degs(theta)), fontsize=this_font_size)
    #plt.show()
    plt.savefig(filename, dpi=300)

plot_T_vs_kL(all_exp_data, theta, omega)
