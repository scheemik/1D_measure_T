"""
Author: Mikhail Schee

Writes parameters for one simulation.

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

# n_layers = 4
kL_start = 0
kL_stop  = 5 #/n_layers

kL_arr = np.linspace(kL_start, kL_stop, n_sims)

theta  = np.arctan(1.0) # 0.7853981634

###############################################################################
# Write out parameters based on sim ID

filename = "params.py"
params_file = open(filename, "w+")

# params_file.write("n_layers = " + str(n_layers))
params_file.write("\nkL = " + str(kL_arr[sim_id]))
# params_file.write("\nkL = " + str(1.55))
params_file.write("\ntheta = " + str(theta))

params_file.close()
