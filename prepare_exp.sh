#!/bin/bash

# Run this script to prepare an experiment directory
# Takes in optional arguments:
#	$ sh prepare_exp.sh -e <name of experiment> 				 Default: current datetime
#								  		-s <number of simulations>       Default: 5

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

###############################################################################
# Variables

# Name of experiment plotting file
plot_exp="plot_exp_data.py"
# Name of folder that contains simulation code files
CODE_FILES="_code_files"

###############################################################################
# Check if experiments folder exists
if [ -e _experiments ]
then
	echo 'Experiment folder exists'
else
	echo 'Experiment folder not found. Aborting prepare_exp script'
	exit 1
fi

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
		echo "Not overwriting ${EXP}. Aborting prepare_exp script"
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

# Make a general copy of all the code files
cp -r $CODE_FILES _experiments/${EXP}

###############################################################################
# Get each individual simulation's directory ready to run

# Make directory for all the simulations
mkdir _experiments/${EXP}/_simulations

# Loop over each simulation ID
for (( id=0; id<$SIMS; id++ ))
do
	# Pad ID with zeros so all ID's have 3 digits
	printf -v NID "%03d" $id
	# Generate simulation name
	SIM_NAME="${NID}_${EXP}"
	# Copy code files to experiment directory
	cp -r $CODE_FILES _experiments/${EXP}/_simulations
	# Rename code file folder for specific simulation
	mv _experiments/${EXP}/_simulations/${CODE_FILES} _experiments/${EXP}/_simulations/${SIM_NAME}
done
echo ''
