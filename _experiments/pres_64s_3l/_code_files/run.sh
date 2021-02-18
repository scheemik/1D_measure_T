#!/bin/bash

# A bash script to run the Dedalus python code
#		Assumes it is run inside a simulation directory
# Takes in optional arguments:
#	$ sh run.sh -e <name of experiment> 				Default: test_exp
#             -i <sim ID>          						Default: 0
#							-n <name of simulation>         Default: ID + name of experiment
#							-s <number of simulations>      Default: 1
#							-c <number of cores>            Default: 1
#							-a <ask before deleting>				Default: False
#							-r <run simulation>             Default: False
#             -m <merge h5 files>             Default: False
#							-o <post-process data>          Default: False
#							-p <plot simulation>            Default: False
#             -g <create gif>                 Default: False
#             -v <create video, mp4>          Default: False

# Current datetime
DATETIME=`date +"%Y-%m-%d_%Hh%M"`

# Having a ":" after a flag means an option is required to invoke that flag
while getopts e:i:n:s:c:armopgv option
do
	case "${option}"
		in
		e) EXP=${OPTARG};;
		i) ID=${OPTARG};;
		n) NAME=${OPTARG};;
		s) SIMS=${OPTARG};;
		c) CORES=${OPTARG};;
		a) ASK=true;;
		r) RUN=true;;
    m) MER=true;;
    o) PRO=true;;
    p) PLT=true;;
    g) GIF=true;;
    v) VID=true;;
	esac
done

# check to see what arguments were passed
if [ -z "$EXP" ]
then
	EXP="test_exp"
	echo "-e, No name specified, using EXP=$EXP"
else
  echo "-e, Name specified, using EXP=$EXP"
fi
if [ -z "$ID" ]
then
	ID=0
	echo "-i, No simulation ID specified, using ID=$ID"
else
  echo "-i, Simulation ID specified, using ID=$ID"
fi
if [ -z "$NAME" ]
then
	# Pad ID with zeros
	printf -v NID "%03d" $ID
	NAME="${NID}_${EXP}"
	echo "-n, No simulation name specified, using NAME=$NAME"
else
  echo "-n, Simulation name specified, using NAME=$NAME"
fi
if [ -z "$SIMS" ]
then
	SIMS=1
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
if [ "$ASK" = true ]
then
	echo "-a, Will ask before overwriting files"
else
	ASK=false
	echo "-a, Will not ask before overwriting files"
fi
if [ "$RUN" = true ]
then
	echo "-r, Run the simulation"
fi
if [ "$MER" = true ]
then
	echo "-m, Merge the data"
fi
if [ "$PRO" = true ]
then
	echo "-o, Post-process data"
fi
if [ "$PLT" = true ]
then
	echo "-p, Plot the simulation"
else
	PLT=false
fi
if [ "$GIF" = true ]
then
	echo "-g, Create gif of simulation"
fi
if [ "$VID" = true ]
then
	echo "-v, Create video (mp4) of the simulation"
fi

###############################################################################
# The command and arguments for running scripts with mpi
# mpiexec_command="mpiexec"
# The version of python to use
python_command="python3"
# Path to snapshot files
snapshot_name="snapshots"
snapshot_path="${snapshot_name}"
# Name of output directory
output_dir="outputs"
# Path to frames
frames_path='frames'
# Clean up the snapshots after merging
cleanup_snapshots="True"

# Name of the script to write the parameters
params_script='write_params.py'
# Name of parameter file
params_file='params.py'
# Name of the main code file
code_file='main.py'
# Name of switchboard file
switchboard="switchboard"
switch_file="${switchboard}.py"
# Name of merging file
merge_file="merge.py"
# Name of post processing file
post_process="post_process.py"
# Name of gif creation file
gif_cre_file="create_gif.py"

###############################################################################
# parameters file
echo ''
echo "${ID}-Writing parameters file-"
# Check if parameter file already exists
if [ -e $params_file ]
then
	echo "Found $params_file"
else
	if [ -e $params_script ]
	then
		${python_command} $params_script $NAME $ID $SIMS
	else
		echo "Cannot find $params_script, aborting run script"
		exit 1
	fi
fi

###############################################################################
# run simulation
if [ "$RUN" = true ]
then
	echo ''
	echo "${ID}-Running script-"
	# Check if snapshots already exist. If so, remove them
	if [ -e $snapshot_path ]
	then
		echo "Removing old snapshots"
		rm -rf "${snapshot_path}"		# Careful, easy to accidentally remove $NAME dir
	fi
	echo "Running Dedalus script"
	${python_command} $code_file $NAME $ID $PLT
  echo ""
	echo 'Done running script'
fi

###############################################################################
# merge h5 files (snapshots)
if [ "$MER" = true ]
then
	echo ''
	echo "${ID}-Merging snapshots-"
	# Check to make sure snapshots folder exists
	echo "Checking for snapshots in directory: $snapshot_path"
	if [ -e $snapshot_path ]
	then
		echo "Found snapshots"
	else
		echo "Cannot find snapshots. Aborting run script"
		exit 1
	fi
	# Check if snapshots have already been merged
	if [ -e $snapshot_path/consolidated_analysis.h5 ] || [ -e $snapshot_path/snapshots_s1.h5 ] || [ -e $snapshot_path/snapshots_s01.h5 ]
	then
		echo "Snapshots already merged"
	else
		echo "Merging snapshots"
		${python_command} $merge_file $snapshot_path --cleanup=$cleanup_snapshots
	fi
  echo 'Done merging snapshots'

	# Reformat snapshot file names if necessary
	if [ -e $snapshot_path/snapshots_s10.h5 ]
	then
		echo 'Reformatting snapshot file names'
		#for i in {1..9..1}
		for ((i=1;i<=9;++i))
		do
			old_name=${snapshot_path}/snapshots_s${i}.h5
			if [ -e $old_name ]
			then
				new_name=${snapshot_path}/snapshots_s0${i}.h5
				echo "Reformatting: ${old_name} to ${new_name}"
				mv $old_name $new_name
			fi
		done
	fi
fi

###############################################################################
# post-process data, make plots if requested through PLT
if [ "$PRO" = true ]
then
	echo ''
	echo "${ID}-Post processing-"
	# Check to make sure snapshots exists
	echo "Checking for snapshots in directory: $snapshot_path"
	if [ -e $snapshot_path/consolidated_analysis.h5 ] || [ -e $snapshot_path/snapshots_s1.h5 ] || [ -e $snapshot_path/snapshots_s01.h5 ]
	then
		echo "Snapshots found"
	else
		echo "Cannot find snapshots. Aborting run script"
		exit 1
	fi
	echo 'Running post processing script'
	${python_command} $post_process $NAME $PLT $snapshot_path/*.h5
	echo 'Done post processing'
fi

###############################################################################
# create gif
if [ "$GIF" = true ]
then
	echo ''
	echo "${ID}-Creating gif-"
	gif_name="${DATETIME}.gif"
	# Check if output directory exists
	if [ ! -e $output_dir ]
	then
		echo "Creating $output_dir directory"
		mkdir $output_dir
	fi
	# Check if gis already exists
	if [ -e $output_dir/$gif_name ]
	then
		echo "Overwriting $gif_name"
		rm $output_dir/$gif_name
	fi
	files=/$frames_path/*
	if [ -e $frames_path ] && [ ${#files[@]} -gt 0 ]
	then
		echo "Executing gif script"
		${python_command} $gif_cre_file $NAME $output_dir/$gif_name $frames_path
	else
		echo "No frames found"
	fi
	echo 'Done with gif creation'
fi

echo ''
echo "${ID}-Done running simulation-"
echo ''
