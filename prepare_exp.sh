#!/bin/bash

# Run this script to prepare an experiment directory
# Takes in optional arguments:
#	$ sh prepare_exp.sh -e <name of experiment> 				 Default: current datetime
#								  		-s <number of simulations>      Default: 5

echo ""
echo "Preparing experiment directory"

# Current datetime
DATETIME=`date +"%Y-%m-%d_%Hh%M"`

# Having a ":" after a flag means an option is required to invoke that flag
while getopts e:s: option
do
	case "${option}"
		in
		e) EXP=${OPTARG};;
		s) SIMS=${OPTARG};;
	esac
done

# check to see if arguments were passed
if [ -z "$EXP" ]
then
	EXP=$DATETIME
	echo "-e, No name specified, using EXP=$EXP"
else
  echo "-e, Name specified, using EXP=$EXP"
fi
if [ -z "$SIMS" ]
then
	SIMS=5
	echo "-s, No number of simulations specified, using SIMS=$SIMS"
else
  echo "-s, Number of simulations specified, using SIMS=$SIMS"
fi
if [ -z "$CORES" ]
then
	CORES=1
	echo "-c, No number of cores specified, using CORES=$CORES"
else
  echo "-c, Number of cores specified, using CORES=$CORES"
fi
if [ "$RUN" = r ]
then
	ARGS=true
	echo "-r, Run the simulation"
fi
if [ "$MER" = m ]
then
	ARGS=true
	echo "-m, Merge the data"
fi
if [ "$PRO" = o ]
then
	ARGS=true
	echo "-o, Post-process data"
fi
if [ "$PLT" = p ]
then
	ARGS=true
	echo "-p, Plot the simulation"
fi
if [ "$GIF" = g ]
then
	ARGS=true
	echo "-g, Create gif of simulation"
fi
if [ "$VID" = v ]
then
	ARGS=true
	echo "-v, Create video (mp4) of the simulation"
fi

###############################################################################
# Miscellaneous variables

# Name of experiment plotting file
plot_exp="plot_exp_data.py"

###############################################################################
# Code file variables needed for each simulation

# Name of the run script
run_file="run.sh"
# Name of the script to write the unique parameters file
params_script='write_params.py'
# Name of the main code file
code_file='main.py'
# Name of switchboard file
switch_file="switchboard.py"
# Name of merging file
merge_file="merge.py"
# Helper code files
helper_funcs="helper_functions.py"
helper_funcs_CD="helper_functions_CD.py"
# Name of post processing file
post_process="post_process.py"
# Name of slice plotting file
plot_file="plot_slices.py"
# Name of gif creation file
gif_cre_file="create_gif.py"
# Group all the code files for ease of calling
CODE_FILES="$run_file $params_script $code_file $switch_file $merge_file $helper_funcs $helper_funcs_CD $post_process $plot_file $gif_cre_file"

###############################################################################

# Check if directory for experiment exists already
if [ -d _experiments/${EXP} ]
then
	# If it does exist already, ask if it should be overwritten
	# Begin loop waiting for user to confirm
	while true
	do
		read -r -p "${EXP} exists already. Overwrite [y/n]? " input
		case $input in
			[yY][eE][sS]|[yY])
		echo "Overwriting ${EXP}"
		rm -rf _experiments/${EXP}
		break
		;;
			[nN][oO]|[nN])
		echo "Not overwriting ${EXP}. Aborting script"
		exit 1
		break
				;;
			*)
		echo "Invalid input"
		;;
		esac
	done
fi

# Make the experiment directory
mkdir _experiments/${EXP}

# Copy experiment plotting file to experiment directory
cp $plot_exp _experiments/${EXP}

###############################################################################
# Get each individual simulation's directory ready to run

# Loop over each simulation ID
for (( id=0; id<$SIMS; id++ ))
do
	# Pad ID with zeros so all ID's have 3 digits
	printf -v NID "%03d" $id
	# Generate simulation name
	SIM_NAME="${NID}_${EXP}"
	# Make directory for that simulation
	mkdir _experiments/${EXP}/${SIM_NAME}
	# Copy code files to simulation directory
	cp $CODE_FILES _experiments/${EXP}/${SIM_NAME}
done
