# Helper functions for Dedalus experiment
"""
Description:
This contains helper functions for the Dedalus code so the same version of
    functions can be called by multiple scripts
"""

import numpy as np
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.colors as clrs
from matplotlib import ticker
from dedalus.extras.plot_tools import quad_mesh, pad_limits

###############################################################################
# The analytical solution to the tranmission coefficient in the linear,
#   non-rotating, non-viscous case. From Sutherland and Yewchuk 2004, equation 2.4
def SY_eq2_4(theta, kL):
    return 1 / (1 + (np.sinh(kL)/np.sin(2*theta))**2 )

###############################################################################
# Measuring the Transmission Coefficient

def find_nearest_index(array, value):
    """
    Assumes array is sorted, returns the integer of the index which contains
        the value closest to the value given. Assumes values are the same type
    array       a 1 D array of values
    value       the given value in question
    """
    class Error(Exception):
        pass
    # Check to make sure array is in the correct order
    if array[0] > array[-1]:
        raise Error("Array sorted incorrectly")
        return 0
    # Use a bisection search function from numpy
    idx = np.searchsorted(array, value, side="left")
    # Assumes there is a problem if the result is an endpoint
    if idx == 0:
        raise Error("Value to the left of index range")
        return None
    elif idx == len(array):
        raise Error("Value to the right of index range")
        return None
    # If closer to the right index, return index
    elif np.abs(value - array[idx]) < np.abs(value - array[idx-1]):
        return idx
    # If closer to the left index, return index-1
    else:
        return idx-1

def AAcc(data):
    """
    Returns the data multiplied by the complex conjugate of the data
    """
    return data * np.conj(data)

# depricated function
def max_amp_at_z(data, T_skip=None, T=None, t=None):
    """
    data        1D array of wave field for certain z, AAcc
    T_skip      oscillation periods to skip before measuring
    T           oscillation period in seconds
    t           array of time values
    """
    if T_skip != None and T_skip != 0.0:
        # find index of time after skip interval
        idt   = find_nearest_index(t, T_skip*T)
        arr_A = data[idt:]
        max_amp = max(arr_A)
    else:
        max_amp = max(data)
    # Multiply element wise by the complex conjugate, find maximum
    return max_amp

def avg_within_bounds(data, x_array=None, l_bound=None, r_bound=None):
    """
    Returns average value of an array between left and right bounds, if provided
    data        1D array of data, assumed to be real valued
    x_array     axis array for data
    l_bound     left bound, in x (not an index)
    r_bound     right bound, in x (not an index)
    """
    if x_array != None:
        # find index of left and right bounds
        idx_l    = find_nearest_index(x_array, l_bound)
        idx_r    = find_nearest_index(x_array, r_bound)
        data_avg = np.mean(data[idx_l:idx_r])
    else:
        data_avg = np.mean(data)
    return data_avg

def pull_depth_data(depth, data, z):
    """
    Returns a 1D array of the data at a specified depth

    depth       specified depth in meters
    data        2D array of wavefield data (t,z)
    z           array of z values
    """
    # Find the index of the z's closest to specified depth
    idx_depth = find_nearest_index(z, depth)
    # Pull relevant depths from the data, take absolute value
    data_at_depth = data[:][idx_depth]
    return data_at_depth

def measure_T(data, z, z_I, z_T, T_skip=0.0, T=None, t=None):
    """
    Measures the transmission coefficient for a filtered wave field
    Expects the data to be only one direction of propagation

    data        array of the psi wavefield, complex valued
    z           array of z values
    z_I         depth at which to measure incident wave
    z_T         depth at which to measure transmitted wave
    T_skip      oscillation periods to skip before measuring
    T           oscillation period in seconds
    t           array of time values
    """
    # Pull relevant depths from the data
    arr_I = pull_depth_data(z_I, data, z)
    arr_T = pull_depth_data(z_T, data, z)
    # Multiply element wise by the complex conjugate, assume imaginary parts go to 0
    AAcc_I = AAcc(arr_I).real
    AAcc_T = AAcc(arr_T).real
    # Find amplitude values for incident and transmitted waves
    use_max = False
    if use_max:
        I_ = max_amp_at_z(AAcc_I, T_skip, T, t)
        T_ = max_amp_at_z(AAcc_T, T_skip, T, t)
    else:
        I_ = avg_within_bounds(AAcc_I)
        T_ = avg_within_bounds(AAcc_T)
    # The ratio of T_/I_ is the Transmission Coefficient
    return I_, T_, AAcc_I, AAcc_T


###############################################################################
# Takes an exponential number and returns a string formatted nicely for latex
#   Expects numbers in the format 7.0E+2
def latex_exp(num):
    """
    Takes in a number, returns a string

    num         number to be formatted
    """
    if (isinstance(num, int)):
        # integer type, don't reformat
        return str(num)
    else:
        float_str = "{:.1E}".format(num)
        if "E" in float_str:
            base, exponent = float_str.split("E")
            exp = int(exponent)
            b   = float(base)
            str1 = '$'
            if (exp == -1):
                str1 = str1 + str(b/10.0)
            elif (exp == 0):
                str1 = str1 + str(base)
            elif (exp == 1):
                str1 = str1 + str(b*10.0)
            elif (exp == 2):
                str1 = str1 + str(b*100.0)
            else:
                str1 = str1 + str(base) + r'\cdot10^{' + str(exp) + '}'
            str1 = str1 + '$'
            return r"{0}".format(str1)
        else:
            return float_str

# takes in a value in radians and outputs the value in degrees
def rad_to_degs(num):
    """
    Takes in a number, returns a string

    num         number to be converted to degrees
    """
    deg = int(round(math.degrees(num)))
    deg_str = str(deg) + r'$^\circ$'
    return deg_str

###############################################################################

def extended_stop_time(sim_time_stop, dt):
    # Find number of time slices
    nt = sim_time_stop // dt
    i = 1
    # Find what power of 2 nt should round up to
    while (nt > 2**i):
        i += 1
    # Calculate new sim_time_stop
    new_sim_time_stop = dt * (2**i)
    return new_sim_time_stop, int(nt)

def trim_data(z_array, t_array, data, z0_dis=None, zf_dis=None, T_cutoff=None, T=None):
    """
    Trims data, removes time before T_cutoff and space outside display domain

    z_array     array of z values
    t_array     array of time values (in seconds)
    data        data to be trimmed (z, t)
    z0_dis      top of vertical structure extent
    zf_dis      bottom of vertical structure extent
    T_cutoff    integer, oscillations to cut off of beginning of simulation
    T           oscillation period (in seconds)
    """
    ## Trim time before T_cutoff
    if T_cutoff != None and T != None:
        # Find time in seconds of cutoff
        cutoff_time = T_cutoff * T
        # Find index of cutoff time
        idx_cutoff = find_nearest_index(t_array, cutoff_time)
        # Trim time
        t_data = data[:,idx_cutoff:]
        trimmed_t_array = t_array[idx_cutoff:]
    else:
        t_data = data
        trimmed_t_array = t_array
    ## Trim space to only encompass display bounds
    if z0_dis != None and zf_dis != None:
        # Find indices of display bounds
        idx_z0 = find_nearest_index(z_array, z0_dis)
        idx_zf = find_nearest_index(z_array, zf_dis)
        # Trim space ### careful of the order, zf is bottom, so should have lower index?
        trimmed_data = t_data[idx_zf:idx_z0,:]
        trimmed_z_array = z_array[idx_zf:idx_z0]
    else:
        trimmed_data = t_data
        trimmed_z_array = z_array
    return trimmed_z_array, trimmed_t_array, trimmed_data

def trim_data_z(z_array, data, z0_dis, zf_dis):
    """
    Trims data, removes space outside display domain

    z_array     1D array of z values
    data        1D or 2D data to be trimmed (z, t)
    z0_dis      top of vertical structure extent
    zf_dis      bottom of vertical structure extent
    """
    ## Trim space to only encompass display bounds
    #   Find indices of display bounds
    idx_z0 = find_nearest_index(z_array, z0_dis)
    idx_zf = find_nearest_index(z_array, zf_dis)
    # Trim space ### careful of the order, zf is bottom, so should have lower index
    #   Check whether data is 1D or 2D
    if len(data.shape) == 2:
        trimmed_data = data[idx_zf:idx_z0,:]
    elif len(data.shape) == 1:
        trimmed_data = data[idx_zf:idx_z0]
    else:
        print("trim_data_z only accepts 1D or 2D arrays")
    trimmed_z_array = z_array[idx_zf:idx_z0]
    return trimmed_z_array, trimmed_data

def trim_data_t(t_array, data, T_cutoff, T):
    """
    Trims data, removes time before T_cutoff

    t_array     1D array of time values (in seconds)
    data        2D data to be trimmed (z, t)
    T_cutoff    integer, oscillations to cut off of beginning of simulation
    T           oscillation period (in seconds)
    """
    ## Trim time before T_cutoff
    #   Find time in seconds of cutoff
    cutoff_time = T_cutoff * T
    # print("cutoff_time = ",cutoff_time)
    # Find index of cutoff time
    idx_cutoff = find_nearest_index(t_array, cutoff_time)
    # Trim time
    trimmed_data    = data[:,idx_cutoff:]
    trimmed_t_array = t_array[idx_cutoff:]
    return trimmed_t_array, trimmed_data

###############################################################################

def set_fig_axes(heights, widths, fig_ratio=0.5, share_x_axis=None, share_y_axis=None):
    """
    Creates fig and axes objects based on desired heights and widths of subplots
    Ex: if widths=[1,5], there will be 2 columns, the 1st 1/5 the width of the 2nd

    heights     array of integers for subplot height ratios, len=rows
    widths      array of integers for subplot width  ratios, len=cols
    fig_ratio   ratio of height to width of overall figure
    share_x_axis bool whether the subplots should share their x axes
    share_y_axis bool whether the subplots should share their y axes
    """
    # Set aspect ratio of overall figure
    w, h = mpl.figure.figaspect(fig_ratio)
    # Find rows and columns of subplots
    rows = len(heights)
    cols = len(widths)
    # This dictionary makes each subplot have the desired ratios
    # The length of heights will be nrows and likewise len(widths)=ncols
    plot_ratios = {'height_ratios': heights,
                   'width_ratios': widths}
    # Determine whether to share x or y axes
    if share_x_axis == None and share_y_axis == None:
        if rows == 1 and cols != 1: # if only one row, share y axis
            share_x_axis = False
            share_y_axis = True
        elif rows != 1 and cols == 1: # if only one column, share x axis
            share_x_axis = True
            share_y_axis = False
        else:                       # otherwise, forget about it
            share_x_axis = False
            share_y_axis = False
    # Set ratios by passing dictionary as 'gridspec_kw', and share y axis
    return plt.subplots(figsize=(w,h), nrows=rows, ncols=cols, gridspec_kw=plot_ratios, sharex=share_x_axis, sharey=share_y_axis)

# Create and format a colorbar
def set_colorbar(data, im, plt, axis):
    """
    Creates formats a colorbar for the provided axis
    Ideally, it would format the ticks with my Latex function, but I
        haven't figured that out yet

    data        data array which is put into the colormap
    im          image object of colormap
    plt         plot object
    axis        axis of plot object on which the colorbar will be placed
    """
    # Find max of absolute value for colorbar for limits symmetric around zero
    cmax = max(abs(data.flatten()))
    if cmax==0.0:
        cmax = 0.001 # to avoid the weird jump with the first frame
    # Set upper and lower limits on colorbar
    im.set_clim(-cmax, cmax)
    # Add colorbar to im
    cbar = plt.colorbar(im, ax=axis)#, format=ticker.FuncFormatter(latex_exp))
    cbar.ax.ticklabel_format(style='sci', scilimits=(-2,2), useMathText=True)

# A standard title for most of my plots
def add_plot_title(fig, title_str, kL, theta, omega):
    """
    Adds overall title (suptitle) to given figure
    fig         figure object on which to add the title
    title_str   string of title
    others      values of parameters to include in title
    """
    param_formated_str = latex_exp(kL)+', '+rad_to_degs(theta)+', '+latex_exp(omega)
    fig.suptitle(r'%s, $(kL,\theta,\omega)$=(%s)' %(title_str, param_formated_str))

# Background profile in N_0
def BP_n_layers(z, z0_str, n, L, R_i):
    """
    z           array of z values
    z0_str      top of vertical structure extent
    n           number of layers
    L           layer thickness
    R_i         ratio of interface to layer thickness
    """
    # create array of 1's the same length as z
    BP_array = z*0+1
    # start the first step's top at the top of the structure
    step_top = z0_str
    # Loop across each i step
    for i in range(n):
        # Set bottom of step
        step_bot = step_top - L
        # Set the vertical structure to 0 for that step
        for j in range(len(BP_array)):
            if z[j] < step_top and z[j] > step_bot:
                BP_array[j] = 0
        # Set step top for next step
        step_top = step_bot - R_i*L
    return BP_array

# # Mask to just keep display domain
# def make_DD_mask(z, z0_dis, zf_dis):
#     """
#     z           array of z values
#     z0_dis      bottom of display domain
#     zf_dis      top of display domain
#     """
#     # create blank array the same size as z
#     DD_array = z*0
#     # Set display domain range to have value 1
#     for i in range(len(DD_array)):
#         if z[i] < zf_dis and z[i] > z0_dis:
#             DD_array[i] = 1
#     return DD_array


def add_lines_to_ax(ax, z_I=None, z_T=None, z0_dis=None, zf_dis=None, T_cutoff=None):
    """
    Adds straight lines, either vertical or horizontal depending, to the provided axis

    ax          axis for plot
    z_I         depth at which to measure Incident wave
    z_T         depth at which to measure Transmitted wave
    z0_dis      top of vertical structure extent
    zf_dis      bottom of vertical structure extent
    T_cutoff    integer, oscillations to cut off of beginning of simulation
    """
    # Horizontal lines for incident and transmission depths
    if z_I != None: # separated to allow for plotting of I_ and T_ separately
        ax.axhline(y=z_I, color=my_clrs['incident'], linestyle='--')
    if z_T != None:
        ax.axhline(y=z_T, color=my_clrs['transmission'], linestyle='--')
    # Horizontal lines for display depths
    if z0_dis != None and zf_dis != None:
        line_color = my_clrs['black']
        ax.axhline(y=z0_dis, color=line_color, linestyle='--')
        ax.axhline(y=zf_dis, color=line_color, linestyle='--')
    # Vertical line for transient cutoff
    if T_cutoff != None:
        line_color = my_clrs['black']
        ax.axvline(x=T_cutoff, color=line_color, linestyle='--')

def set_plot_bounds(ax, plot_full_x, plot_full_y, z_array=None, z0_dis=None, zf_dis=None, t_array=None, T=None, T_cutoff=None):
    """
    Sets plot bounds, transient oscillations cutoff and display domain

    ax          axis for plot
    plot_full_x True or False, include and transient period?
    plot_full_y True or False, include depths outside display domain?
    z_array     array of z values
    t_array     array of time values (in seconds)
    T           oscillation period (in seconds)
    z0_dis      top of vertical structure extent
    zf_dis      bottom of vertical structure extent
    T_cutoff    integer, oscillations to cut off of beginning of simulation
    """
    if plot_full_x:
        if isinstance(t_array, np.ndarray) and T != None:
            ax.set_xlim([t_array[0]/T,t_array[-1]/T])
    else:
        if isinstance(t_array, np.ndarray) and T != None and T_cutoff != None:
            ax.set_xlim([T_cutoff,t_array[-1]/T])
    if plot_full_y:
        if isinstance(z_array, np.ndarray):
            ax.set_ylim([z_array[0],z_array[-1]])
    else:
        if z0_dis != None and zf_dis != None:
            ax.set_ylim([zf_dis, z0_dis]) # Careful of this order, zf is less than z0


# Plot background profile
def plot_BP(ax, BP, z, omega=None):
    """
    Plots the background profile as a function of z

    ax          axis object for plot
    BP          1D array of background profile in N_0
    z           1D array of z values
    omega       frequency of wave
    """
    ax.plot(BP, z, color=my_clrs['N_0'], label=r'$N_0$')
    ax.set_xlabel(r'$N_0$ (s$^{-1}$)')
    ax.set_ylabel(r'$z$ (m)')
    ax.set_title(r'Background Profile')
    ax.set_ylim([min(z),max(z)])
    if omega != None:
        ax.axvline(x=omega, color=my_clrs['omega'], linestyle='--', label=r'$\omega$')
        ax.legend()

###############################################################################

def plot_v_profiles(z_array, BP_array, bf_array, sp_array, kL=None, theta=None, omega=None, z_I=None, z_T=None, z0_dis=None, zf_dis=None, plot_full_x=True, plot_full_y=True, title_str='Forced 1D Wave', filename='f_1D_windows.png'):
    """
    Plots the vertical profiles: stratification, boundary forcing, sponge layer
        Note: all arrays must be imput as full-domain, no trimmed versions

    z_array     1D array of z values
    BP_array    1D array of the background profile values in z
    bf_array    1D array of boundary forcing window
    sp_array    1D array of sponge layer window
    kL          Non-dimensional number relating wavelength and layer thickness
    theta       Angle at which wave is incident on stratification structure
    omega       frequency of wave
    z_I         depth to measure incident wave
    z_T         depth to meausre transmitted wave
    z0_dis      top of vertical structure extent
    zf_dis      bottom of vertical structure extent
    plot_full...True or False, include depths outside display and transient period?
    """
    # Set figure and axes for plot
    fig, axes = set_fig_axes([1], [1,4], 0.75)
    # Plot the background profile of N_0
    plot_BP(axes[0], BP_array, z_array, omega)
    # Plot boudnary forcing and sponge layer windows
    axes[1].plot(bf_array, z_array, color=my_clrs['F_bf'], label='Boundary forcing')
    axes[1].plot(sp_array, z_array, color=my_clrs['F_sp'], label='Sponge layer')
    # Add horizontal lines
    add_lines_to_ax(axes[0], z_I=z_I, z_T=z_T, z0_dis=z0_dis, zf_dis=zf_dis)
    add_lines_to_ax(axes[1], z_I=z_I, z_T=z_T, z0_dis=z0_dis, zf_dis=zf_dis)
    # Set plot bounds to include full domain or not
    set_plot_bounds(axes[0], plot_full_x, plot_full_y, z_array=z_array, z0_dis=z0_dis, zf_dis=zf_dis)
    set_plot_bounds(axes[1], plot_full_x, plot_full_y, z_array=z_array, z0_dis=z0_dis, zf_dis=zf_dis)
    # Add labels
    axes[1].set_xlabel('Amplitude')
    axes[1].set_title(r'Windows')
    axes[1].legend()
    #
    fig.suptitle(r'%s' %(title_str))
    plt.savefig(filename)

###############################################################################
# Main plotting functions
###############################################################################

def plot_z_vs_t(z_array, t_array, T, data, BP_array, kL, theta, omega, z_I=None, z_T=None, z0_dis=None, zf_dis=None, plot_full_x=True, plot_full_y=True, T_cutoff=0.0, c_map='RdBu_r', title_str='Forced 1D Wave', filename='f_1D_wave.png'):
    """
    Plots the data as a colormap on z vs t with the vertical profile included to the left

    z_array     array of z values
    t_array     array of time values (in seconds)
    T           oscillation period (in seconds)
    data        data to be plotted on colormap (assmued real valued)
    BP_array    array of the background profile values in z
    kL          Non-dimensional number relating wavelength and layer thickness
    theta       Angle at which wave is incident on stratification structure
    omega       frequency of wave
    z_I         depth to measure incident wave
    z_T         depth to meausre transmitted wave
    z0_dis      top of vertical structure extent
    zf_dis      bottom of vertical structure extent
    plot_full...True or False, include depths outside display and transient period?
    T_cutoff    integer, oscillations to cut off of beginning of simulation
    """
    # Set figure and axes for plot
    fig, axes = set_fig_axes([1], [1,5])
    # Plot background profiles
    plot_BP(axes[0], BP_array, z_array, omega)
    # Create meshes and image of colormap
    xmesh, ymesh = quad_mesh(x=t_array/T, y=z_array)
    im = axes[1].pcolormesh(xmesh, ymesh, data, cmap=c_map)
    # Set colorbar
    set_colorbar(data, im, plt, axes[1])
    # Add straight lines to profile plot
    add_lines_to_ax(axes[0], z_I=z_I, z_T=z_T, z0_dis=z0_dis, zf_dis=zf_dis)
    # Add straight lines to wavefield plot
    add_lines_to_ax(axes[1], z_I=z_I, z_T=z_T, z0_dis=z0_dis, zf_dis=zf_dis, T_cutoff=T_cutoff)
    # Set plot bounds to include full domain or not
    set_plot_bounds(axes[0], plot_full_x, plot_full_y, z_array=z_array, z0_dis=z0_dis, zf_dis=zf_dis)
    set_plot_bounds(axes[1], plot_full_x, plot_full_y, z_array=z_array, z0_dis=z0_dis, zf_dis=zf_dis, t_array=t_array, T=T, T_cutoff=T_cutoff)
    # Add titles and labels
    axes[1].set_xlabel(r'$t/T$')
    axes[1].set_title(r'$\Psi$ (m$^2$/s)')
    add_plot_title(fig, title_str, kL, theta, omega)
    #plt.show()
    plt.savefig(filename)

###############################################################################

def plot_A_of_I_T(z_array, t_array, T, dn_array, kL, theta, omega, z_I, z_T, plot_full_x=True, plot_full_y=True, T_cutoff=0.0, title_str='Forced 1D Wave', filename='f_1D_A_of_I_T.png'):
    """
    Plots the amplitude of the downward propagating wave as a function of time
        for both the incident and transmission depths

    z_array     array of z values
    t_array     array of time values (in seconds)
    T           oscillation period (in seconds)
    dn_array    wavefield data for downward propagating wave (complex valued)
    z_I         depth to measure incident wave
    z_T         depth to meausre transmitted wave
    kL          Non-dimensional number relating wavelength and layer thickness
    theta       Angle at which wave is incident on stratification structure
    omega       frequency of wave
    z_I         depth to measure incident wave
    z_T         depth to meausre transmitted wave
    plot_full...True or False, include depths outside display and transient period?
    T_cutoff    integer, oscillations to cut off of beginning of simulation
    """
    # Set figure and axes for plot
    fig, axes = set_fig_axes([1,1], [1])
    # Find the AA^* for incident and transmitted wave depths
    I_, T_, AAcc_I, AAcc_T = measure_T(dn_array, z_array, z_I, z_T, T_skip=T_cutoff, T=T, t=t_array)
    # Plot lines of AAcc for incident and transmission depths
    axes[0].plot(t_array/T, AAcc_I, color=my_clrs['b'], label=r'$I$')
    axes[1].plot(t_array/T, AAcc_T, color=my_clrs['b'], label=r'$T$')
    # Add horizontal lines to show the overall value of each depth
    #   Add vertical lines to show where the transient ends
    add_lines_to_ax(axes[0], z_I=I_, T_cutoff=T_cutoff)
    add_lines_to_ax(axes[1], z_T=T_, T_cutoff=T_cutoff)
    # Set plot bounds to include full domain or not
    set_plot_bounds(axes[0], plot_full_x, plot_full_y, t_array=t_array, T=T, T_cutoff=T_cutoff)
    set_plot_bounds(axes[1], plot_full_x, plot_full_y, t_array=t_array, T=T, T_cutoff=T_cutoff)
    # Add labels and titles
    axes[1].set_xlabel(r'$t/T$')
    axes[0].set_ylabel(r'Amplitude')
    axes[1].set_ylabel(r'Amplitude')
    axes[0].set_title(r'$z_I=%s$' %(z_I))
    axes[1].set_title(r'$z_T=%s$' %(z_T))
    add_plot_title(fig, title_str, kL, theta, omega)
    #plt.show()
    plt.savefig(filename)

###############################################################################

def plot_AA_for_z(z_array, BP_array, dn_array, kL, theta, omega, z_I=None, z_T=None, z0_dis=None, zf_dis=None, plot_full_x=True, plot_full_y=True, title_str='Forced 1D Wave', filename='f_1D_AA_for_z.png'):
    """
    Plots the average of the downward-propagating wave's AA^* for all depths

    z_array     array of z values
    dn_array    wavefield data for downward propagating wave (AA^*, complex parts ~0)
    BP_array    array of the background profile values in z
    kL          Non-dimensional number relating wavelength and layer thickness
    theta       Angle at which wave is incident on stratification structure
    omega       frequency of wave
    z_I         depth to measure incident wave
    z_T         depth to meausre transmitted wave
    z0_dis      top of vertical structure extent
    zf_dis      bottom of vertical structure extent
    plot_full...True or False, include depths outside display and transient period?
    """
    # Set figure and axes for plot
    fig, axes = set_fig_axes([1], [1,4], 0.75)
    # Plot the background profile of N_0
    plot_BP(axes[0], BP_array, z_array, omega)
    # Find maximum value of the field times its complex conjugate for each z
    I_and_T_for_z = z_array*0.0
    for i in range(len(z_array)):
        # avg_within_bounds assumes real valued input
        I_and_T_for_z[i] = avg_within_bounds(dn_array[i].real)
    # Plot boudnary forcing and sponge layer windows
    axes[1].plot(I_and_T_for_z, z_array, color=my_clrs['incident'], label='Amp')
    # Add horizontal lines
    add_lines_to_ax(axes[0], z_I=z_I, z_T=z_T, z0_dis=z0_dis, zf_dis=zf_dis)
    add_lines_to_ax(axes[1], z_I=z_I, z_T=z_T, z0_dis=z0_dis, zf_dis=zf_dis)
    # Set plot bounds to include full domain or not
    set_plot_bounds(axes[0], plot_full_x, plot_full_y, z_array=z_array, z0_dis=z0_dis, zf_dis=zf_dis)
    set_plot_bounds(axes[1], plot_full_x, plot_full_y, z_array=z_array, z0_dis=z0_dis, zf_dis=zf_dis)
    # Add labels
    axes[1].set_xlabel('Amplitude')
    #axes[1].set_ylabel(r'$z$')
    axes[1].set_title(r'Windows')
    axes[1].legend()
    add_plot_title(fig, title_str, kL, theta, omega)
    plt.savefig(filename)

###############################################################################

###############################################################################

def sort_k_coeffs(arr, nz):
    if len(arr.shape) == 1:
        sorted_arr = np.zeros(nz-1)
        sorted_arr[0:(nz//2-1)] = arr[(nz//2):(nz-1)]
        sorted_arr[(nz//2-1):(nz-1)] = arr[0:(nz//2)]
        return sorted_arr
    else:
        for i in range(arr.shape[1]):
            arr[i][:] = sort_k_coeffs(arr[i][:], nz)
        return arr

def plot_spectral(k_array, f_array, r_array, i_array, kL, theta, omega, c_map='viridis', title_str='Forced 1D Wave', filename='f_1D_wave_spectra.png'):
    """
    k_array     1D array of wavenumber values
    f_array     1D array of frequency values
    r_array     2D real part of data
    i_array     2D imaginary part of data
    kL          Non-dimensional number relating wavelength and layer thickness
    theta       Angle at which wave is incident on stratification structure
    omega       frequency of wave
    """
    # Set figure and axes for plot
    fig, axes = set_fig_axes([1], [1,1])
    # Make grid meshes and plot colormaps on both plots
    xmesh, ymesh = quad_mesh(x=f_array, y=k_array)
    im_r = axes[0].pcolormesh(xmesh, ymesh, r_array, cmap=c_map)
    im_i = axes[1].pcolormesh(xmesh, ymesh, i_array, cmap=c_map)
    # Add colorbar to im
    cbar_r = plt.colorbar(im_r, ax=axes[0])#, format=ticker.FuncFormatter(latex_exp))
    cbar_r.ax.ticklabel_format(style='sci', scilimits=(-2,2), useMathText=True)
    cbar_i = plt.colorbar(im_i, ax=axes[1])#, format=ticker.FuncFormatter(latex_exp))
    cbar_i.ax.ticklabel_format(style='sci', scilimits=(-2,2), useMathText=True)
    axes[0].set_ylabel(r'$k$ (m$^{-1}$)')
    # Add labels
    axes[0].set_xlabel(r'$f$ (s$^{-1}$)')
    axes[1].set_xlabel(r'$f$ (s$^{-1}$)')
    axes[0].set_title(r'real')
    axes[1].set_title(r'imag')
    # Set limits
    axes[0].set_ylim([-100,100])
    axes[1].set_ylim([-100,100])
    add_plot_title(fig, title_str, kL, theta, omega)
    #plt.show()
    plt.savefig(filename)

def plot_k_f_spectra(z_array, dz, t_array, dt, T, k_array, f_array, k_data, f_data, kL, theta, omega, z_I, z_T, plot_full_x=True, plot_full_y=True, T_cutoff=0.0, title_str='Forced 1D Wave', filename='f_1D_k_and_f.png'):
    """
    Plots the spectral content of the data in wavenumber and frequency space for
        the T_cutoff time and z_I depth, respectively

    z_array     1D array of z values
    dz
    t_array     1D array of time values (in seconds)
    dt
    T           oscillation period (in seconds)
    k_array     1D array of wavenumber values
    f_array     1D array of frequency values
    k_data      2D data in wavenumber space (k, t)
    f_data      2D data in frequency  space (z, f)
    z_I         depth to measure incident wave
    z_T         depth to meausre transmitted wave
    kL          Non-dimensional number relating wavelength and layer thickness
    theta       Angle at which wave is incident on stratification structure
    omega       frequency of wave
    z_I         depth to measure incident wave
    z_T         depth to meausre transmitted wave
    plot_full...True or False, include depths outside display and transient period?
    T_cutoff    integer, oscillations to cut off of beginning of simulation
    """
    # Set figure and axes for plot (fig_ratio determines visibility of some lines, be careful)
    fig, axes = set_fig_axes([1], [1,1], fig_ratio=0.33, share_x_axis=False, share_y_axis=False)
    # 1D or 2D?
    plot_1D = False
    # Find wavenumber axis
    k_array = np.fft.fftshift(np.fft.fftfreq(len(z_array), dz))
    if plot_1D:
        # Take just a slice of the wavenumber data at T_cutoff time
        idx_t = find_nearest_index(t_array, T_cutoff*T)
        # Normalize by number of points squared to get power spectrum
        plot_k = np.fft.fftshift(k_data[:,idx_t]) / (len(k_array)**2) # allows comparison across datasets
        # Plot lines of wavenumber spectrum, real and imaginary
        axes[0].plot(k_array, plot_k.real, color=my_clrs['b'], label=r'real')
        axes[0].plot(k_array, plot_k.imag, color=my_clrs['w'], label=r'imag')
    else:
        c_map='viridis'
        # Make the mesh grid to be (t, k)
        xmesh, ymesh = quad_mesh(x=t_array, y=k_array)
        # Normalize by number of points squared to get power spectrum
        plot_k = (np.fft.fftshift(k_data) / (len(k_array)**2)).real
        # Create colormap of spectral data
        im_k = axes[0].pcolormesh(xmesh, ymesh, plot_k, cmap=c_map)
        set_colorbar(plot_k, im_k, plt, axes[0])
        k_dis_bound = 100
        axes[0].set_ylim([-k_dis_bound,k_dis_bound])
    #
    # Find frequency axis
    f_array = np.fft.fftshift(np.fft.fftfreq(len(t_array), dt))
    if plot_1D:
        # Take just a slice of the frequency data at z_I depth
        idx_z = find_nearest_index(z_array, z_I)
        # Normalize by number of points squared to get power spectrum
        plot_f = np.fft.fftshift(f_data[idx_z,:]) / (len(f_array)**2) # allows comparison across datasets
        # Plot lines of frequency spectrum, real and imaginary
        axes[1].plot(f_array, plot_f.real, color=my_clrs['b'], label=r'real')
        axes[1].plot(f_array, plot_f.imag, color=my_clrs['w'], label=r'imag')
    else:
        c_map='viridis'
        # Make the mesh grid to be (f, z)
        xmesh, ymesh = quad_mesh(x=f_array, y=z_array)
        # Normalize by number of points squared to get power spectrum
        plot_f = np.fft.fftshift(f_data) / (len(f_array)**2)
        # Create colormap of spectral data
        im_f = axes[1].pcolormesh(xmesh, ymesh, plot_f.real, cmap=c_map)
        set_colorbar(plot_f, im_f, plt, axes[1])
    #
    # Add labels and titles
    if plot_1D:
        axes[0].legend()
        axes[1].legend()
        axes[0].set_xlabel(r'$k$')
        axes[0].set_ylabel(r'Amplitude')
        axes[1].set_xlabel(r'$f$')
        # axes[1].set_ylabel(r'Amplitude')
        axes[0].set_title(r'$t_{1}=%s T$' %(T_cutoff))
        axes[1].set_title(r'$z_I=%s$' %(z_I))
    else:
        axes[0].set_xlabel(r'$t$')
        axes[0].set_ylabel(r'$k$')
        axes[1].set_xlabel(r'$f$')
        axes[1].set_ylabel(r'$z$')
        axes[0].set_title(r'Wavenumber spectrum')
        axes[1].set_title(r'Frequency spectrum')
    add_plot_title(fig, title_str, kL, theta, omega)
    #plt.show()
    plt.savefig(filename)

###############################################################################

# Make a plot for one time slice
def plot_task(ax, time_i, task_j, z_ax, dsets):
    # plot line of psi vs. z
    im = ax.plot(dsets[task_j][time_i][1], z_ax, color=my_clrs['b'])
    # Find max of absolute value for data to make symmetric around zero
    xmax = max(abs(max(dsets[task_j][time_i][1].flatten())), abs(min(dsets[task_j][time_i][1].flatten())))
    if xmax==0.0:
        xmax = 0.001 # to avoid the weird jump with the first frame
    # format range of plot extent
    ax.set_xlim(-xmax, xmax)
    ax.set_ylim(z_ax[0], z_ax[-1])

def find_tick_locations(x0, xf, n_ticks=3):
    locs = [0]*n_ticks
    locs[0]  = x0
    locs[-1] = xf
    if n_ticks <= 3 and x0 < 0 and xf>0:
        locs[1] = 0.0
    else:
        spacing = (xf-x0)/(n_ticks-1)
        for i in range(1,n_ticks):
            locs[i] = x0 + spacing*i
    return locs

def format_labels_and_ticks(ax, hori_label, n_ticks=3, tick_formatter=latex_exp):
    # add labels
    ax.set_xlabel(hori_label)
    # fix horizontal ticks
    x0, xf = ax.get_xlim()
    tick_locs = find_tick_locations(x0, xf, n_ticks)
    ax.xaxis.set_ticks(tick_locs)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(tick_formatter))

def format_labels_and_ticks_v(ax, vert_label):
    # add labels
    ax.set_ylabel(vert_label)
    # fix vertical ticks
    y0, yf = ax.get_ylim()
    ax.yaxis.set_ticks([y0, 0.0, yf])
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(latex_exp))

###############################################################################
# Plotting colors from style guide

CUSTOM_COLORS ={'lightcornflowerblue2': '#a4c2f4',
                'lightred3': '#f2c1c1'}

TAB_COLORS =   {'tab:blue': '#1f77b4',
                'tab:orange': '#ff7f0e',
                'tab:green': '#2ca02c',
                'tab:red': '#d62728',
                'tab:purple': '#ffffff',
                'tab:brown': '#ffffff',
                'tab:pink': '#ffffff',
                'tab:gray': '#ffffff',
                'tab:olive': '#ffffff',
                'tab:cyan': '#ffffff'
                }

CSS4_COLORS =  {'aliceblue': '#F0F8FF',
                'antiquewhite': '#FAEBD7',
                'aqua': '#00FFFF',
                'aquamarine': '#7FFFD4',
                'azure': '#F0FFFF',
                'beige': '#F5F5DC',
                'bisque': '#FFE4C4',
                'black': '#000000',
                'blanchedalmond': '#FFEBCD',
                'blue': '#0000FF',
                'blueviolet': '#8A2BE2',
                'brown': '#A52A2A',
                'burlywood': '#DEB887',
                'cadetblue': '#5F9EA0',
                'chartreuse': '#7FFF00',
                'chocolate': '#D2691E',
                'coral': '#FF7F50',
                'cornflowerblue': '#6495ED',
                'cornsilk': '#FFF8DC',
                'crimson': '#DC143C',
                'cyan': '#00FFFF',
                'darkblue': '#00008B',
                'darkcyan': '#008B8B',
                'darkgoldenrod': '#B8860B',
                'darkgray': '#A9A9A9',
                'darkgreen': '#006400',
                'darkgrey': '#A9A9A9',
                'darkkhaki': '#BDB76B',
                'darkmagenta': '#8B008B',
                'darkolivegreen': '#556B2F',
                'darkorange': '#FF8C00',
                'darkorchid': '#9932CC',
                'darkred': '#8B0000',
                'darksalmon': '#E9967A',
                'darkseagreen': '#8FBC8F',
                'darkslateblue': '#483D8B',
                'darkslategray': '#2F4F4F',
                'darkslategrey': '#2F4F4F',
                'darkturquoise': '#00CED1',
                'darkviolet': '#9400D3',
                'deeppink': '#FF1493',
                'deepskyblue': '#00BFFF',
                'dimgray': '#696969',
                'dimgrey': '#696969',
                'dodgerblue': '#1E90FF',
                'firebrick': '#B22222',
                'floralwhite': '#FFFAF0',
                'forestgreen': '#228B22',
                'fuchsia': '#FF00FF',
                'gainsboro': '#DCDCDC',
                'ghostwhite': '#F8F8FF',
                'gold': '#FFD700',
                'goldenrod': '#DAA520',
                'gray': '#808080',
                'green': '#008000',
                'greenyellow': '#ADFF2F',
                'grey': '#808080',
                'honeydew': '#F0FFF0',
                'hotpink': '#FF69B4',
                'indianred': '#CD5C5C',
                'indigo': '#4B0082',
                'ivory': '#FFFFF0',
                'khaki': '#F0E68C',
                'lavender': '#E6E6FA',
                'lavenderblush': '#FFF0F5',
                'lawngreen': '#7CFC00',
                'lemonchiffon': '#FFFACD',
                'lightblue': '#ADD8E6',
                'lightcoral': '#F08080',
                'lightcyan': '#E0FFFF',
                'lightgoldenrodyellow': '#FAFAD2',
                'lightgray': '#D3D3D3',
                'lightgreen': '#90EE90',
                'lightgrey': '#D3D3D3',
                'lightpink': '#FFB6C1',
                'lightsalmon': '#FFA07A',
                'lightseagreen': '#20B2AA',
                'lightskyblue': '#87CEFA',
                'lightslategray': '#778899',
                'lightslategrey': '#778899',
                'lightsteelblue': '#B0C4DE',
                'lightyellow': '#FFFFE0',
                'lime': '#00FF00',
                'limegreen': '#32CD32',
                'linen': '#FAF0E6',
                'magenta': '#FF00FF',
                'maroon': '#800000',
                'mediumaquamarine': '#66CDAA',
                'mediumblue': '#0000CD',
                'mediumorchid': '#BA55D3',
                'mediumpurple': '#9370DB',
                'mediumseagreen': '#3CB371',
                'mediumslateblue': '#7B68EE',
                'mediumspringgreen': '#00FA9A',
                'mediumturquoise': '#48D1CC',
                'mediumvioletred': '#C71585',
                'midnightblue': '#191970',
                'mintcream': '#F5FFFA',
                'mistyrose': '#FFE4E1',
                'moccasin': '#FFE4B5',
                'navajowhite': '#FFDEAD',
                'navy': '#000080',
                'oldlace': '#FDF5E6',
                'olive': '#808000',
                'olivedrab': '#6B8E23',
                'orange': '#FFA500',
                'orangered': '#FF4500',
                'orchid': '#DA70D6',
                'palegoldenrod': '#EEE8AA',
                'palegreen': '#98FB98',
                'paleturquoise': '#AFEEEE',
                'palevioletred': '#DB7093',
                'papayawhip': '#FFEFD5',
                'peachpuff': '#FFDAB9',
                'peru': '#CD853F',
                'pink': '#FFC0CB',
                'plum': '#DDA0DD',
                'powderblue': '#B0E0E6',
                'purple': '#800080',
                'rebeccapurple': '#663399',
                'red': '#FF0000',
                'rosybrown': '#BC8F8F',
                'royalblue': '#4169E1',
                'saddlebrown': '#8B4513',
                'salmon': '#FA8072',
                'sandybrown': '#F4A460',
                'seagreen': '#2E8B57',
                'seashell': '#FFF5EE',
                'sienna': '#A0522D',
                'silver': '#C0C0C0',
                'skyblue': '#87CEEB',
                'slateblue': '#6A5ACD',
                'slategray': '#708090',
                'slategrey': '#708090',
                'snow': '#FFFAFA',
                'springgreen': '#00FF7F',
                'steelblue': '#4682B4',
                'tan': '#D2B48C',
                'teal': '#008080',
                'thistle': '#D8BFD8',
                'tomato': '#FF6347',
                'turquoise': '#40E0D0',
                'violet': '#EE82EE',
                'wheat': '#F5DEB3',
                'white': '#FFFFFF',
                'whitesmoke': '#F5F5F5',
                'yellow': '#FFFF00',
                'yellowgreen': '#9ACD32'}


my_clrs       =  {'b': TAB_COLORS['tab:blue'],
                  'w': (1, 0, 0),               # - r
                  'u': (0, 0, 1),               # - b
                  'v'  : (0, 0.5, 0),           # - g
                  'p': CSS4_COLORS['plum'],
                  'diffusion': CSS4_COLORS['peru'],
                  'viscosity': CSS4_COLORS['peru'],
                  'N_0': TAB_COLORS['tab:blue'],
                  'rho': CSS4_COLORS['slateblue'],
                  'advection': CSS4_COLORS['indianred'],
                  'coriolis': CSS4_COLORS['teal'],
                  'omega': CSS4_COLORS['slategray'],
                  'F_bf': '#008080',            # - teal
                  'F_sp': '#CD853F',            # - peru
                  'temperature': '#B22222',     # - firebrick
                  'salinity': '#4682B4',        # - steelblue
                  'incident': '#8A2BE2',        # - blueviolet
                  'transmission': '#4169E1',    # - royalblue
                  'reflection': '#FF6347',      # - tomato
                  'linear': CSS4_COLORS['forestgreen'],
                  'nonlinear': CSS4_COLORS['indianred'],
                  'arctic': CSS4_COLORS['cornflowerblue'],
                  'cold-fresh': CUSTOM_COLORS['lightcornflowerblue2'],
                  'warm-salty': CUSTOM_COLORS['lightred3'],
                  'black': (0, 0, 0),
                  'white': (1, 1, 1)}
