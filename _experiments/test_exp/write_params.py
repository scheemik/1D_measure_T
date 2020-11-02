"""
Writes parameters for one simulation. Run with $ python3 post_process.py NAME PLOT_CHECKS snapshots/*.h5

Usage:
    write_params.py NAME ID N_SIMS

Options:
    NAME            # name of the experiment run from -n
    ID              # Simulation ID number
    N_SIMS          # Number of simulations

"""
###############################################################################
###############################################################################
# Import standard(ish) libraries and functions
import numpy as np

# Parse input parameters
from docopt import docopt
args = docopt(__doc__)
run_name    = args['NAME']      # Simulation name, used to route filepaths of plots
sim_id      = int(args['ID'])
n_sims      = int(args['N_SIMS'])

###############################################################################
# Make arrays of parameters

mL_start = 0
mL_stop  = 5

mL_arr = np.linspace(mL_start, mL_stop, n_sims)

###############################################################################
# Write out parameters based on sim ID

filename = run_name + "/params.py"
params_file = open(filename, "w")

params_file.write("mL = " + str(mL_arr[sim_id]))
params_file.write("\ntheta = 0.7")

params_file.close()
